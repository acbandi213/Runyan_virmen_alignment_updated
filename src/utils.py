"""Utility functions for seeding, checkpointing, and plotting."""

import random
import numpy as np
import torch


def set_seed(seed: int = 42) -> None:
    """Set random seeds for reproducibility across torch, numpy, and random.

    Args:
        seed: Random seed value.
    """
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def save_checkpoint(model, optimizer, epoch: int, metrics: dict, path: str) -> None:
    """Save model checkpoint to disk.

    Args:
        model: PyTorch model.
        optimizer: PyTorch optimizer.
        epoch: Current epoch number.
        metrics: Dictionary of training metrics.
        path: File path to save checkpoint.
    """
    torch.save({
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(),
        'epoch': epoch,
        'metrics': metrics,
    }, path)


def load_checkpoint(path: str) -> dict:
    """Load a checkpoint from disk.

    Args:
        path: File path to load checkpoint from.

    Returns:
        Dictionary with keys: model_state_dict, optimizer_state_dict, epoch, metrics.
    """
    checkpoint = torch.load(path, map_location='cpu', weights_only=False)
    return checkpoint
