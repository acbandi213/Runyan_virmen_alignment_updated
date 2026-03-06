"""Unit tests for the CTRNN model."""

import pytest
import torch

from src.models import CTRNN


class TestCTRNN:

    def test_output_shape(self):
        """Forward pass should return correct shapes."""
        model = CTRNN(input_size=8, hidden_size=256, output_size=2, seed=42)
        x = torch.randn(4, 20, 8)
        outputs, hidden, h_final = model(x)
        assert outputs.shape == (4, 20, 2)
        assert hidden.shape == (4, 20, 256)
        assert h_final.shape == (4, 256)

    def test_hidden_state_persistence(self):
        """Different h_init should produce different outputs."""
        model = CTRNN(input_size=8, hidden_size=64, seed=42)
        model.eval()
        x = torch.randn(1, 20, 8)
        _, _, h1 = model(x, h_init=torch.zeros(1, 64))
        _, _, h2 = model(x, h_init=torch.randn(1, 64))
        assert not torch.allclose(h1, h2)

    def test_eval_deterministic(self):
        """Eval mode should be deterministic (no noise)."""
        model = CTRNN(input_size=8, hidden_size=64, sigma_rec=0.1, seed=42)
        model.eval()
        x = torch.randn(1, 20, 8)
        out1, _, _ = model(x)
        out2, _, _ = model(x)
        assert torch.allclose(out1, out2)

    def test_train_stochastic(self):
        """Train mode should be stochastic (noise added)."""
        model = CTRNN(input_size=8, hidden_size=64, sigma_rec=0.1, seed=42)
        model.train()
        x = torch.randn(1, 20, 8)
        out1, _, _ = model(x)
        out2, _, _ = model(x)
        assert not torch.allclose(out1, out2)

    def test_gamma_value(self):
        """gamma should be dt/tau = 20/100 = 0.2."""
        model = CTRNN()
        assert model.gamma == pytest.approx(0.2)

    def test_init_hidden_zeros(self):
        """init_hidden should return zeros."""
        model = CTRNN(hidden_size=128)
        h = model.init_hidden(batch_size=3)
        assert h.shape == (3, 128)
        assert (h == 0).all()

    def test_none_h_init(self):
        """Passing h_init=None should work (defaults to zeros)."""
        model = CTRNN(input_size=8, hidden_size=64, seed=42)
        model.eval()
        x = torch.randn(2, 10, 8)
        outputs, _, _ = model(x, h_init=None)
        assert outputs.shape == (2, 10, 2)
