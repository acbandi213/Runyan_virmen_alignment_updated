"""CTRNN model following Yang et al. 2019 (Nature Neuroscience).

Dynamics:
    h(t+1) = (1 - γ) * h(t) + γ * f(W_rec @ h(t) + W_in @ u(t) + b) + noise

where γ = Δt/τ, f = Softplus, noise ~ N(0, σ² * 2γ) during training only.
"""

import math

import torch
import torch.nn as nn
import torch.nn.functional as F


class CTRNN(nn.Module):
    """Continuous-time recurrent neural network.

    Args:
        input_size: Dimension of input (default: 8).
        hidden_size: Number of recurrent units (default: 256).
        output_size: Number of output classes (default: 2, left/right).
        dt: Integration timestep in ms (default: 20).
        tau: Membrane time constant in ms (default: 100).
        sigma_rec: Recurrent noise standard deviation (default: 0.05).
        seed: Random seed for weight initialization.
    """

    def __init__(self, input_size: int = 8, hidden_size: int = 256,
                 output_size: int = 2, dt: float = 20.0, tau: float = 100.0,
                 sigma_rec: float = 0.05, seed: int = 42):
        super().__init__()
        self.input_size = input_size
        self.hidden_size = hidden_size
        self.output_size = output_size
        self.dt = dt
        self.tau = tau
        self.gamma = dt / tau  # 0.2
        self.sigma_rec = sigma_rec
        self.noise_scale = sigma_rec * math.sqrt(2.0 * self.gamma)

        torch.manual_seed(seed)

        self.W_in = nn.Linear(input_size, hidden_size, bias=False)
        self.W_rec = nn.Linear(hidden_size, hidden_size, bias=True)
        self.W_out = nn.Linear(hidden_size, output_size, bias=True)
        self.activation = nn.Softplus()

        self._initialize_weights()

    def _initialize_weights(self):
        """Initialize weights.

        W_in: Xavier uniform.
        W_rec: Normal with std 1/sqrt(hidden_size) (near edge of chaos).
        W_out: Xavier uniform.
        Biases: zeros.
        """
        nn.init.xavier_uniform_(self.W_in.weight)
        nn.init.normal_(self.W_rec.weight, mean=0.0,
                        std=1.0 / math.sqrt(self.hidden_size))
        nn.init.xavier_uniform_(self.W_out.weight)
        nn.init.zeros_(self.W_rec.bias)
        nn.init.zeros_(self.W_out.bias)

    def init_hidden(self, batch_size: int) -> torch.Tensor:
        """Initialize hidden state to zeros.

        Args:
            batch_size: Number of parallel sessions.

        Returns:
            Zeros tensor of shape [batch_size, hidden_size].
        """
        return torch.zeros(batch_size, self.hidden_size,
                           device=self.W_rec.weight.device)

    def recurrence(self, u_t: torch.Tensor, h_t: torch.Tensor) -> torch.Tensor:
        """Single-step CTRNN dynamics.

        Args:
            u_t: Input at time t, shape [batch, input_size].
            h_t: Hidden state at time t, shape [batch, hidden_size].

        Returns:
            Hidden state at time t+1, shape [batch, hidden_size].
        """
        pre_act = self.W_rec(h_t) + self.W_in(u_t)
        r_t = self.activation(pre_act)
        h_next = (1.0 - self.gamma) * h_t + self.gamma * r_t

        if self.training and self.sigma_rec > 0:
            noise = torch.randn_like(h_next) * self.noise_scale
            h_next = h_next + noise

        return h_next

    def forward(self, inputs: torch.Tensor,
                h_init: torch.Tensor = None) -> tuple:
        """Forward pass through a sequence of timesteps.

        Args:
            inputs: Shape [batch, n_timesteps, input_size].
            h_init: Initial hidden state [batch, hidden_size]. Zeros if None.

        Returns:
            outputs: Logits [batch, n_timesteps, output_size].
            hidden_states: All states [batch, n_timesteps, hidden_size].
            h_final: Final hidden state [batch, hidden_size].
        """
        batch_size, n_timesteps, _ = inputs.shape

        if h_init is None:
            h_init = self.init_hidden(batch_size)

        h_t = h_init
        hidden_list = []

        for t in range(n_timesteps):
            h_t = self.recurrence(inputs[:, t, :], h_t)
            hidden_list.append(h_t)

        hidden_states = torch.stack(hidden_list, dim=1)
        outputs = self.W_out(hidden_states)

        return outputs, hidden_states, h_t
