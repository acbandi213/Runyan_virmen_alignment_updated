"""Training loops with logging for CTRNN models.

Handles session-level training with dynamic feedback: the model's actual
predictions determine the feedback signal for the next trial, rather than
assuming correct behavior. Hidden state persists across trials within a session.
"""

from typing import Dict, List, Optional, Any

import numpy as np
import torch
import torch.nn as nn

from src.utils import save_checkpoint


def inject_feedback(trial_inputs: torch.Tensor,
                    feedback: List[float]) -> torch.Tensor:
    """Replace feedback channels [6:8] in a trial's input tensor.

    Args:
        trial_inputs: Input tensor [n_timesteps, 8].
        feedback: Two-element list [fb1, fb2].

    Returns:
        Modified copy of the input tensor.
    """
    modified = trial_inputs.clone()
    modified[:, 6:8] = torch.tensor(feedback, dtype=torch.float32)
    return modified


def compute_dynamic_feedback(model_output: torch.Tensor,
                             target: int,
                             timesteps_resp: int = 5) -> List[float]:
    """Compute feedback from model's prediction on a trial.

    Args:
        model_output: Output logits [n_timesteps, output_size].
        target: True target class (0 or 1).
        timesteps_resp: Number of response timesteps to average over.

    Returns:
        [1.0, 0.0] if correct, [0.0, 1.0] if incorrect.
    """
    resp_logits = model_output[-timesteps_resp:].mean(dim=0)
    predicted = resp_logits.argmax().item()
    return [1.0, 0.0] if predicted == target else [0.0, 1.0]


def run_session(model: nn.Module, session, device: torch.device,
                tbptt_trials: Optional[int] = None) -> Dict[str, Any]:
    """Process a session trial-by-trial with dynamic feedback.

    Hidden state persists across trials. Feedback for trial i+1 is based
    on the model's prediction at trial i.

    Args:
        model: CTRNN model.
        session: Task1Session or Task2Session instance.
        device: torch device.
        tbptt_trials: Detach hidden state every N trials. None = full BPTT.

    Returns:
        Dict with all_outputs, all_hidden, all_targets, all_metadata,
        trial_predictions, trial_correct.
    """
    n_trials = len(session)
    h = model.init_hidden(1).to(device)

    all_outputs = []
    all_hidden = []
    all_targets = []
    all_metadata = []
    trial_predictions = []
    trial_correct = []

    feedback = [0.0, 0.0]  # first trial: no feedback

    for i in range(n_trials):
        trial = session[i]
        trial_input = inject_feedback(trial['inputs'], feedback)
        trial_input = trial_input.unsqueeze(0).to(device)
        target = trial['target']

        outputs, hidden_states, h = model(trial_input, h_init=h)

        all_outputs.append(outputs.squeeze(0))
        all_hidden.append(hidden_states.squeeze(0))
        all_targets.append(target)
        all_metadata.append(trial['metadata'])

        with torch.no_grad():
            fb = compute_dynamic_feedback(outputs.squeeze(0), target)
            predicted = outputs.squeeze(0)[-5:].mean(dim=0).argmax().item()
        feedback = fb
        trial_predictions.append(predicted)
        trial_correct.append(predicted == target)

        if tbptt_trials is not None and (i + 1) % tbptt_trials == 0:
            h = h.detach()

    return {
        'all_outputs': all_outputs,
        'all_hidden': all_hidden,
        'all_targets': all_targets,
        'all_metadata': all_metadata,
        'trial_predictions': trial_predictions,
        'trial_correct': trial_correct,
    }


def compute_loss(session_result: Dict[str, Any], model: nn.Module,
                 timesteps_resp: int = 5,
                 l1_lambda: float = 1e-4) -> tuple:
    """Compute cross-entropy on response epoch + L1 on firing rates.

    Args:
        session_result: Output from run_session.
        model: CTRNN model (for activation function).
        timesteps_resp: Number of response timesteps.
        l1_lambda: L1 regularization coefficient.

    Returns:
        (total_loss, metrics_dict) where metrics_dict has
        ce_loss, l1_loss, total_loss, accuracy.
    """
    ce_fn = nn.CrossEntropyLoss()

    resp_outputs = []
    resp_targets = []

    for outputs, target in zip(session_result['all_outputs'],
                               session_result['all_targets']):
        resp_outputs.append(outputs[-timesteps_resp:])
        resp_targets.extend([target] * timesteps_resp)

    resp_outputs = torch.cat(resp_outputs, dim=0)
    resp_targets = torch.tensor(resp_targets, dtype=torch.long,
                                device=resp_outputs.device)

    ce_loss = ce_fn(resp_outputs, resp_targets)

    all_hidden = torch.cat(session_result['all_hidden'], dim=0)
    firing_rates = model.activation(all_hidden)
    l1_loss = l1_lambda * firing_rates.mean()

    total_loss = ce_loss + l1_loss
    accuracy = sum(session_result['trial_correct']) / len(session_result['trial_correct'])

    metrics = {
        'ce_loss': ce_loss.item(),
        'l1_loss': l1_loss.item(),
        'total_loss': total_loss.item(),
        'accuracy': accuracy,
    }
    return total_loss, metrics


def train_epoch(model: nn.Module, task_class, task_kwargs: dict,
                optimizer: torch.optim.Optimizer, device: torch.device,
                n_sessions: int = 8, grad_clip: float = 1.0,
                tbptt_trials: Optional[int] = None,
                l1_lambda: float = 1e-4,
                seed: int = 42) -> Dict[str, float]:
    """Train for one epoch over multiple sessions.

    Args:
        model: CTRNN model.
        task_class: Task1Session or Task2Session.
        task_kwargs: Task kwargs (excluding seed).
        optimizer: Optimizer.
        device: torch device.
        n_sessions: Sessions per epoch.
        grad_clip: Max gradient norm.
        tbptt_trials: Truncated BPTT window. None = full BPTT.
        l1_lambda: L1 coefficient.
        seed: Base seed (each session gets seed + i).

    Returns:
        Dict of averaged metrics.
    """
    model.train()
    epoch_metrics = {'ce_loss': 0, 'l1_loss': 0, 'total_loss': 0, 'accuracy': 0}

    for i in range(n_sessions):
        session = task_class(**task_kwargs, seed=seed + i)

        optimizer.zero_grad()
        result = run_session(model, session, device, tbptt_trials=tbptt_trials)
        loss, metrics = compute_loss(result, model, l1_lambda=l1_lambda)

        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), grad_clip)
        optimizer.step()

        for k in epoch_metrics:
            epoch_metrics[k] += metrics[k]

    for k in epoch_metrics:
        epoch_metrics[k] /= n_sessions

    return epoch_metrics


def evaluate(model: nn.Module, task_class, task_kwargs: dict,
             device: torch.device, n_sessions: int = 4,
             seed: int = 9999) -> Dict[str, Any]:
    """Evaluate model on fresh sessions.

    Args:
        model: CTRNN model.
        task_class: Task1Session or Task2Session.
        task_kwargs: Task kwargs (excluding seed).
        device: torch device.
        n_sessions: Number of eval sessions.
        seed: Base seed.

    Returns:
        Dict with accuracy, ce_loss, per_trial_accuracy, per_trial_metadata,
        and (Task 2) congruent_accuracy, incongruent_accuracy.
    """
    model.eval()
    all_correct = []
    all_metadata = []
    total_ce = 0.0

    with torch.no_grad():
        for i in range(n_sessions):
            session = task_class(**task_kwargs, seed=seed + i)
            result = run_session(model, session, device)
            _, metrics = compute_loss(result, model)

            all_correct.extend(result['trial_correct'])
            all_metadata.extend(result['all_metadata'])
            total_ce += metrics['ce_loss']

    eval_result = {
        'accuracy': sum(all_correct) / len(all_correct),
        'ce_loss': total_ce / n_sessions,
        'per_trial_accuracy': all_correct,
        'per_trial_metadata': all_metadata,
    }

    # Task 2 specific metrics
    if any('is_congruent' in m for m in all_metadata):
        cong = [c for c, m in zip(all_correct, all_metadata)
                if m.get('is_congruent')]
        incong = [c for c, m in zip(all_correct, all_metadata)
                  if not m.get('is_congruent')]
        if cong:
            eval_result['congruent_accuracy'] = sum(cong) / len(cong)
        if incong:
            eval_result['incongruent_accuracy'] = sum(incong) / len(incong)

    return eval_result


def train_model(model: nn.Module, task_class, task_kwargs: dict,
                device: torch.device, n_epochs: int = 100,
                lr: float = 1e-3, grad_clip: float = 1.0,
                n_sessions_train: int = 8, n_sessions_eval: int = 4,
                eval_every: int = 5, tbptt_trials: Optional[int] = None,
                l1_lambda: float = 1e-4, checkpoint_dir: str = 'checkpoints',
                checkpoint_prefix: str = 'model', seed: int = 42,
                verbose: bool = True) -> Dict[str, list]:
    """Full training loop with periodic evaluation and checkpointing.

    Args:
        model: CTRNN model.
        task_class: Task1Session or Task2Session.
        task_kwargs: Task kwargs (excluding seed).
        device: torch device.
        n_epochs: Number of training epochs.
        lr: Learning rate for Adam.
        grad_clip: Max gradient norm.
        n_sessions_train: Sessions per training epoch.
        n_sessions_eval: Sessions per evaluation.
        eval_every: Evaluate every N epochs.
        tbptt_trials: Truncated BPTT window. None = full BPTT.
        l1_lambda: L1 regularization coefficient.
        checkpoint_dir: Directory for checkpoints.
        checkpoint_prefix: Filename prefix for checkpoints.
        seed: Random seed.
        verbose: Print progress.

    Returns:
        History dict with per-epoch metrics lists.
    """
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)

    history = {
        'train_loss': [], 'train_accuracy': [],
        'eval_loss': [], 'eval_accuracy': [],
        'eval_congruent_acc': [], 'eval_incongruent_acc': [],
        'epoch': [],
    }

    for epoch in range(1, n_epochs + 1):
        epoch_seed = seed + epoch * 1000

        train_metrics = train_epoch(
            model, task_class, task_kwargs, optimizer, device,
            n_sessions=n_sessions_train, grad_clip=grad_clip,
            tbptt_trials=tbptt_trials, l1_lambda=l1_lambda,
            seed=epoch_seed
        )

        history['train_loss'].append(train_metrics['total_loss'])
        history['train_accuracy'].append(train_metrics['accuracy'])

        if epoch % eval_every == 0 or epoch == 1:
            eval_result = evaluate(
                model, task_class, task_kwargs, device,
                n_sessions=n_sessions_eval, seed=9999
            )

            history['eval_loss'].append(eval_result['ce_loss'])
            history['eval_accuracy'].append(eval_result['accuracy'])
            history['eval_congruent_acc'].append(
                eval_result.get('congruent_accuracy'))
            history['eval_incongruent_acc'].append(
                eval_result.get('incongruent_accuracy'))
            history['epoch'].append(epoch)

            if verbose:
                msg = (f"Epoch {epoch:3d} | "
                       f"Train Loss: {train_metrics['total_loss']:.4f} | "
                       f"Train Acc: {train_metrics['accuracy']:.3f} | "
                       f"Eval Acc: {eval_result['accuracy']:.3f}")
                if 'congruent_accuracy' in eval_result:
                    msg += (f" | Cong: {eval_result['congruent_accuracy']:.3f}"
                            f" | Incong: {eval_result['incongruent_accuracy']:.3f}")
                print(msg)

        if epoch % eval_every == 0:
            path = f"{checkpoint_dir}/{checkpoint_prefix}_epoch{epoch}.pt"
            save_checkpoint(model, optimizer, epoch,
                            {'history': history}, path)

    path = f"{checkpoint_dir}/{checkpoint_prefix}_final.pt"
    save_checkpoint(model, optimizer, n_epochs, {'history': history}, path)

    return history
