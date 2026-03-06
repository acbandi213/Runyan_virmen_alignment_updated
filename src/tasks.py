"""Task 1 and Task 2 trial generators for the T-maze decision task.

Task 1 (Unisensory block-switching):
    Mice decide left/right based on either visual or auditory cues. Context
    (which modality is relevant) switches in blocks and is explicitly cued.
    The irrelevant modality carries random, uninformative stimuli.

Task 2 (Multisensory implicit context):
    Both modalities carry informative stimuli on every trial. Context switches
    are implicit — the network must infer from reward feedback. Trials are
    congruent (both cues agree) or incongruent (cues conflict).

Input encoding (8 dimensions):
    [0:2] Visual stimulus: [cos(θ), sin(θ)], θ ∈ {0°, 90°}
    [2:4] Auditory stimulus: [-1, 0] left / [1, 0] right
    [4:6] Context cue: [1,0] visual / [0,1] auditory / [0,0] none
    [6:8] Previous trial feedback: [1,0] correct / [0,1] incorrect / [0,0] first

Stimulus-to-choice mapping:
    Visual: 0° → left (target=0), 90° → right (target=1)
    Auditory: left speaker → left (target=0), right speaker → right (target=1)
"""

import math
from typing import List, Dict, Any

import numpy as np
import torch


def _encode_visual(direction: str) -> List[float]:
    """Encode visual stimulus as [cos(θ), sin(θ)].

    Args:
        direction: 'left' (0°) or 'right' (90°).

    Returns:
        Two-element list [cos(θ), sin(θ)].
    """
    if direction == 'left':
        theta = 0.0
    else:
        theta = math.pi / 2.0
    return [math.cos(theta), math.sin(theta)]


def _encode_auditory(direction: str) -> List[float]:
    """Encode auditory stimulus as speaker location.

    Args:
        direction: 'left' or 'right'.

    Returns:
        Two-element list: [-1, 0] for left, [1, 0] for right.
    """
    if direction == 'left':
        return [-1.0, 0.0]
    else:
        return [1.0, 0.0]


def _direction_to_target(direction: str) -> int:
    """Convert direction string to target integer.

    Args:
        direction: 'left' or 'right'.

    Returns:
        0 for left, 1 for right.
    """
    return 0 if direction == 'left' else 1


def _build_context_schedule(n_trials: int, block_size: int,
                            starting_context: str,
                            n_switches: int = None) -> List[Dict[str, Any]]:
    """Build the context schedule for a session.

    Generates block assignments for each trial. Context alternates every
    block_size trials. By default, switching continues for the entire
    session. If n_switches is set, stops switching after that many times.

    Args:
        n_trials: Total number of trials.
        block_size: Trials per block before a switch.
        starting_context: 'visual' or 'auditory'.
        n_switches: Maximum number of context switches. None = unlimited.

    Returns:
        List of dicts with 'context', 'block_number', 'trial_in_block' per trial.
    """
    schedule = []
    current_context = starting_context
    block_number = 0
    trial_in_block = 0
    switches_done = 0

    for _ in range(n_trials):
        schedule.append({
            'context': current_context,
            'block_number': block_number,
            'trial_in_block': trial_in_block,
        })
        trial_in_block += 1

        can_switch = n_switches is None or switches_done < n_switches
        if trial_in_block >= block_size and can_switch:
            current_context = 'auditory' if current_context == 'visual' else 'visual'
            block_number += 1
            trial_in_block = 0
            switches_done += 1

    return schedule


def _build_input_tensor(visual_enc: List[float], auditory_enc: List[float],
                        context_cue: List[float], feedback: List[float],
                        timesteps_fix: int, timesteps_stim: int,
                        timesteps_resp: int) -> torch.Tensor:
    """Construct the input tensor for a single trial.

    Temporal structure:
        Fixation: stimulus channels zero, context cue + feedback active.
        Stimulus: all channels active.
        Response: all channels active (stimuli persist).

    Args:
        visual_enc: 2-element visual encoding.
        auditory_enc: 2-element auditory encoding.
        context_cue: 2-element context cue.
        feedback: 2-element previous trial feedback.
        timesteps_fix: Fixation epoch length.
        timesteps_stim: Stimulus epoch length.
        timesteps_resp: Response epoch length.

    Returns:
        Tensor of shape [total_timesteps, 8].
    """
    total = timesteps_fix + timesteps_stim + timesteps_resp
    inputs = torch.zeros(total, 8)

    vis_t = torch.tensor(visual_enc, dtype=torch.float32)
    aud_t = torch.tensor(auditory_enc, dtype=torch.float32)
    ctx_t = torch.tensor(context_cue, dtype=torch.float32)
    fb_t = torch.tensor(feedback, dtype=torch.float32)

    # Fixation: context cue + feedback only (no stimuli)
    inputs[:timesteps_fix, 4:6] = ctx_t
    inputs[:timesteps_fix, 6:8] = fb_t

    # Stimulus epoch: all channels
    stim_start = timesteps_fix
    stim_end = timesteps_fix + timesteps_stim
    inputs[stim_start:stim_end, 0:2] = vis_t
    inputs[stim_start:stim_end, 2:4] = aud_t
    inputs[stim_start:stim_end, 4:6] = ctx_t
    inputs[stim_start:stim_end, 6:8] = fb_t

    # Response epoch: stimuli persist
    inputs[stim_end:, 0:2] = vis_t
    inputs[stim_end:, 2:4] = aud_t
    inputs[stim_end:, 4:6] = ctx_t
    inputs[stim_end:, 6:8] = fb_t

    return inputs


class Task1Session:
    """Generates a session of unisensory block-switching trials (Task 1).

    In each block, one modality is relevant (informative) and the other is
    irrelevant (random). The context is explicitly cued. Previous trial
    feedback assumes correct behavior during session generation.

    Args:
        n_trials: Total number of trials in the session.
        block_size: Number of trials per context block before a switch.
        n_switches: Maximum context switches. None (default) = unlimited alternating.
        timesteps_fix: Fixation epoch length in timesteps.
        timesteps_stim: Stimulus epoch length in timesteps.
        timesteps_resp: Response epoch length in timesteps.
        seed: Random seed for reproducibility.
    """

    def __init__(self, n_trials: int = 300, block_size: int = 50,
                 n_switches: int = None, timesteps_fix: int = 5,
                 timesteps_stim: int = 10, timesteps_resp: int = 5,
                 seed: int = 42):
        self.n_trials = n_trials
        self.block_size = block_size
        self.n_switches = n_switches
        self.timesteps_fix = timesteps_fix
        self.timesteps_stim = timesteps_stim
        self.timesteps_resp = timesteps_resp
        self.total_timesteps = timesteps_fix + timesteps_stim + timesteps_resp
        self.seed = seed

        self.trials = self._generate_session()

    def _generate_session(self) -> List[Dict[str, Any]]:
        """Generate all trials for the session.

        Returns:
            List of trial dicts with 'inputs', 'target', and 'metadata'.
        """
        rng = np.random.RandomState(self.seed)

        starting_context = rng.choice(['visual', 'auditory'])
        schedule = _build_context_schedule(
            self.n_trials, self.block_size, starting_context,
            n_switches=self.n_switches
        )

        trials = []
        for i in range(self.n_trials):
            ctx = schedule[i]['context']
            block_num = schedule[i]['block_number']
            tib = schedule[i]['trial_in_block']

            # Choose correct direction randomly
            correct_dir = rng.choice(['left', 'right'])

            # Set relevant modality to correct direction, irrelevant is random
            if ctx == 'visual':
                visual_dir = correct_dir
                auditory_dir = rng.choice(['left', 'right'])
            else:
                auditory_dir = correct_dir
                visual_dir = rng.choice(['left', 'right'])

            # Encode stimuli
            vis_enc = _encode_visual(visual_dir)
            aud_enc = _encode_auditory(auditory_dir)

            # Context cue
            ctx_cue = [1.0, 0.0] if ctx == 'visual' else [0.0, 1.0]

            # Previous trial feedback (assume correct, [0,0] for first trial)
            feedback = [0.0, 0.0] if i == 0 else [1.0, 0.0]

            inputs = _build_input_tensor(
                vis_enc, aud_enc, ctx_cue, feedback,
                self.timesteps_fix, self.timesteps_stim, self.timesteps_resp
            )

            target = _direction_to_target(correct_dir)

            trials.append({
                'inputs': inputs,
                'target': target,
                'metadata': {
                    'context': ctx,
                    'block_number': block_num,
                    'trial_in_block': tib,
                    'correct_direction': correct_dir,
                },
            })

        return trials

    def __len__(self) -> int:
        return len(self.trials)

    def __getitem__(self, idx: int) -> Dict[str, Any]:
        return self.trials[idx]


class Task2Session:
    """Generates a session of multisensory implicit-context trials (Task 2).

    Both modalities carry informative stimuli on every trial. Trials are
    congruent (both cues indicate the same direction) or incongruent (cues
    conflict). Context is NOT explicitly cued — the network must infer it
    from reward feedback.

    Args:
        n_trials: Total number of trials in the session.
        block_size: Number of trials per context block before a switch.
        n_switches: Maximum context switches. None (default) = unlimited alternating.
        congruent_ratio: Fraction of trials that are congruent.
        timesteps_fix: Fixation epoch length in timesteps.
        timesteps_stim: Stimulus epoch length in timesteps.
        timesteps_resp: Response epoch length in timesteps.
        seed: Random seed for reproducibility.
    """

    def __init__(self, n_trials: int = 300, block_size: int = 50,
                 n_switches: int = None, congruent_ratio: float = 0.5,
                 timesteps_fix: int = 5, timesteps_stim: int = 10,
                 timesteps_resp: int = 5, seed: int = 42):
        self.n_trials = n_trials
        self.block_size = block_size
        self.n_switches = n_switches
        self.congruent_ratio = congruent_ratio
        self.timesteps_fix = timesteps_fix
        self.timesteps_stim = timesteps_stim
        self.timesteps_resp = timesteps_resp
        self.total_timesteps = timesteps_fix + timesteps_stim + timesteps_resp
        self.seed = seed

        self.trials = self._generate_session()

    def _generate_session(self) -> List[Dict[str, Any]]:
        """Generate all trials for the session.

        Trial generation:
            1. Decide congruent/incongruent based on congruent_ratio.
            2. Sample the correct direction (relevant modality direction).
            3. Set the irrelevant modality to match (congruent) or oppose
               (incongruent) the relevant modality.

        Returns:
            List of trial dicts with 'inputs', 'target', and 'metadata'.
        """
        rng = np.random.RandomState(self.seed)

        starting_context = rng.choice(['visual', 'auditory'])
        schedule = _build_context_schedule(
            self.n_trials, self.block_size, starting_context,
            n_switches=self.n_switches
        )

        trials = []
        for i in range(self.n_trials):
            ctx = schedule[i]['context']
            block_num = schedule[i]['block_number']
            tib = schedule[i]['trial_in_block']

            # Decide congruent vs incongruent
            is_congruent = rng.random() < self.congruent_ratio

            # Choose the correct direction (what the relevant modality says)
            correct_dir = rng.choice(['left', 'right'])
            opposite_dir = 'right' if correct_dir == 'left' else 'left'

            # Assign directions to each modality
            if ctx == 'visual':
                visual_dir = correct_dir
                auditory_dir = correct_dir if is_congruent else opposite_dir
            else:
                auditory_dir = correct_dir
                visual_dir = correct_dir if is_congruent else opposite_dir

            # Encode stimuli
            vis_enc = _encode_visual(visual_dir)
            aud_enc = _encode_auditory(auditory_dir)

            # No context cue in Task 2
            ctx_cue = [0.0, 0.0]

            # Previous trial feedback (assume correct, [0,0] for first trial)
            feedback = [0.0, 0.0] if i == 0 else [1.0, 0.0]

            inputs = _build_input_tensor(
                vis_enc, aud_enc, ctx_cue, feedback,
                self.timesteps_fix, self.timesteps_stim, self.timesteps_resp
            )

            target = _direction_to_target(correct_dir)

            trials.append({
                'inputs': inputs,
                'target': target,
                'metadata': {
                    'context': ctx,
                    'block_number': block_num,
                    'trial_in_block': tib,
                    'is_congruent': bool(is_congruent),
                    'visual_direction': visual_dir,
                    'auditory_direction': auditory_dir,
                    'correct_direction': correct_dir,
                },
            })

        return trials

    def __len__(self) -> int:
        return len(self.trials)

    def __getitem__(self, idx: int) -> Dict[str, Any]:
        return self.trials[idx]
