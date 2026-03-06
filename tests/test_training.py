"""Unit tests for training infrastructure."""

import torch
import pytest

from src.models import CTRNN
from src.tasks import Task1Session, Task2Session
from src.training import (inject_feedback, compute_dynamic_feedback,
                          run_session, compute_loss)


class TestFeedback:

    def test_inject_feedback_replaces_channels(self):
        """inject_feedback should set channels [6:8] at all timesteps."""
        inputs = torch.zeros(20, 8)
        modified = inject_feedback(inputs, [1.0, 0.0])
        assert (modified[:, 6] == 1.0).all()
        assert (modified[:, 7] == 0.0).all()
        assert (modified[:, :6] == 0).all()

    def test_inject_feedback_does_not_modify_original(self):
        """inject_feedback should return a clone, not modify in-place."""
        inputs = torch.zeros(20, 8)
        _ = inject_feedback(inputs, [1.0, 0.0])
        assert (inputs[:, 6:8] == 0).all()

    def test_dynamic_feedback_correct(self):
        """Should return [1, 0] when prediction matches target."""
        output = torch.zeros(20, 2)
        output[-5:, 1] = 10.0  # predict class 1
        fb = compute_dynamic_feedback(output, target=1)
        assert fb == [1.0, 0.0]

    def test_dynamic_feedback_incorrect(self):
        """Should return [0, 1] when prediction doesn't match target."""
        output = torch.zeros(20, 2)
        output[-5:, 0] = 10.0  # predict class 0
        fb = compute_dynamic_feedback(output, target=1)
        assert fb == [0.0, 1.0]


class TestRunSession:

    def test_output_keys(self):
        """run_session should return expected keys."""
        model = CTRNN(input_size=8, hidden_size=32, seed=42)
        session = Task1Session(n_trials=10, seed=42)
        result = run_session(model, session, torch.device('cpu'))
        expected = {'all_outputs', 'all_hidden', 'all_targets',
                    'all_metadata', 'trial_predictions', 'trial_correct'}
        assert set(result.keys()) == expected

    def test_output_lengths(self):
        """All lists should have n_trials elements."""
        model = CTRNN(input_size=8, hidden_size=32, seed=42)
        n = 15
        session = Task1Session(n_trials=n, seed=42)
        result = run_session(model, session, torch.device('cpu'))
        assert len(result['all_outputs']) == n
        assert len(result['all_hidden']) == n
        assert len(result['all_targets']) == n
        assert len(result['trial_predictions']) == n
        assert len(result['trial_correct']) == n

    def test_output_shapes(self):
        """Individual trial outputs should have correct shapes."""
        model = CTRNN(input_size=8, hidden_size=32, seed=42)
        session = Task1Session(n_trials=5, seed=42)
        result = run_session(model, session, torch.device('cpu'))
        assert result['all_outputs'][0].shape == (20, 2)
        assert result['all_hidden'][0].shape == (20, 32)


class TestComputeLoss:

    def test_loss_is_scalar(self):
        """compute_loss should return a scalar loss tensor."""
        model = CTRNN(input_size=8, hidden_size=32, seed=42)
        session = Task1Session(n_trials=5, seed=42)
        result = run_session(model, session, torch.device('cpu'))
        loss, metrics = compute_loss(result, model)
        assert loss.dim() == 0
        assert loss.requires_grad

    def test_metrics_keys(self):
        """Metrics dict should have expected keys."""
        model = CTRNN(input_size=8, hidden_size=32, seed=42)
        session = Task1Session(n_trials=5, seed=42)
        result = run_session(model, session, torch.device('cpu'))
        _, metrics = compute_loss(result, model)
        assert set(metrics.keys()) == {'ce_loss', 'l1_loss', 'total_loss', 'accuracy'}

    def test_backward_pass(self):
        """Loss should support backward pass without error."""
        model = CTRNN(input_size=8, hidden_size=32, seed=42)
        session = Task1Session(n_trials=5, seed=42)
        result = run_session(model, session, torch.device('cpu'))
        loss, _ = compute_loss(result, model)
        loss.backward()
        # Check gradients exist
        assert model.W_rec.weight.grad is not None
