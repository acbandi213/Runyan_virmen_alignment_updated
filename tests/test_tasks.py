"""Unit tests for Task 1 and Task 2 trial generators."""

import numpy as np
import pytest

from src.tasks import Task1Session, Task2Session


class TestTask1Session:
    """Tests for Task1Session trial generator."""

    def test_input_shape(self):
        """Each trial input tensor should be [20, 8]."""
        session = Task1Session(n_trials=50, seed=42)
        for trial in session.trials:
            assert trial['inputs'].shape == (20, 8), (
                f"Expected shape (20, 8), got {trial['inputs'].shape}"
            )

    def test_irrelevant_modality_uncorrelated(self):
        """Irrelevant modality should not predict target (|r| < 0.15)."""
        session = Task1Session(n_trials=500, block_size=50, n_switches=5, seed=42)

        # Collect irrelevant modality values and targets per context
        irrel_values = []
        targets = []

        for trial in session.trials:
            ctx = trial['metadata']['context']
            target = trial['target']
            # Extract stimulus encoding from the stimulus epoch (timestep 5)
            stim = trial['inputs'][5]

            if ctx == 'visual':
                # Auditory is irrelevant: channel [2] encodes direction
                # [-1, 0] for left, [1, 0] for right
                irrel_values.append(stim[2].item())
            else:
                # Visual is irrelevant: channel [1] encodes direction
                # sin(0°)=0 for left, sin(90°)=1 for right
                irrel_values.append(stim[1].item())

            targets.append(target)

        irrel_values = np.array(irrel_values)
        targets = np.array(targets, dtype=float)
        correlation = np.abs(np.corrcoef(irrel_values, targets)[0, 1])
        assert correlation < 0.15, (
            f"Irrelevant modality correlation with target = {correlation:.3f}, "
            f"expected |r| < 0.15"
        )

    def test_relevant_modality_predicts_target(self):
        """Relevant modality should perfectly predict target."""
        session = Task1Session(n_trials=500, block_size=50, n_switches=5, seed=42)
        correct = 0

        for trial in session.trials:
            ctx = trial['metadata']['context']
            target = trial['target']
            stim = trial['inputs'][5]  # stimulus epoch

            if ctx == 'visual':
                # Visual: sin(θ) > 0.5 → right (1), else left (0)
                predicted = 1 if stim[1].item() > 0.5 else 0
            else:
                # Auditory: channel[2] > 0 → right (1), else left (0)
                predicted = 1 if stim[2].item() > 0 else 0

            if predicted == target:
                correct += 1

        accuracy = correct / len(session.trials)
        assert accuracy == 1.0, (
            f"Relevant modality prediction accuracy = {accuracy:.3f}, expected 1.0"
        )

    def test_context_cue_matches_modality(self):
        """Context cue should correctly reflect the active modality."""
        session = Task1Session(n_trials=300, seed=42)

        for trial in session.trials:
            ctx = trial['metadata']['context']
            # Context cue is on channels [4:6], check during stimulus epoch
            cue = trial['inputs'][5, 4:6].tolist()

            if ctx == 'visual':
                assert cue == [1.0, 0.0], (
                    f"Visual context should have cue [1,0], got {cue}"
                )
            else:
                assert cue == [0.0, 1.0], (
                    f"Auditory context should have cue [0,1], got {cue}"
                )

    def test_number_of_switches_explicit(self):
        """When n_switches is set explicitly, should stop after that many."""
        for n_sw in [1, 2, 3, 5]:
            session = Task1Session(
                n_trials=500, block_size=50, n_switches=n_sw, seed=42
            )
            contexts = [t['metadata']['context'] for t in session.trials]
            switch_count = sum(
                1 for i in range(1, len(contexts))
                if contexts[i] != contexts[i - 1]
            )
            assert switch_count == n_sw, (
                f"Expected {n_sw} switches, got {switch_count}"
            )

    def test_unlimited_switching_default(self):
        """Default n_switches=None should alternate every block_size trials."""
        session = Task1Session(n_trials=300, block_size=50, seed=42)
        contexts = [t['metadata']['context'] for t in session.trials]
        switch_count = sum(
            1 for i in range(1, len(contexts))
            if contexts[i] != contexts[i - 1]
        )
        expected = (300 // 50) - 1  # 5 switches for 6 blocks
        assert switch_count == expected, (
            f"Expected {expected} switches with unlimited, got {switch_count}"
        )

    def test_fixation_epoch_no_stimulus(self):
        """During fixation, stimulus channels should be zero."""
        session = Task1Session(n_trials=10, seed=42)
        for trial in session.trials:
            fix_inputs = trial['inputs'][:5, :4]  # first 5 timesteps, stim channels
            assert (fix_inputs == 0).all(), "Stimulus channels should be zero during fixation"

    def test_feedback_first_trial(self):
        """First trial should have feedback [0, 0]."""
        session = Task1Session(n_trials=10, seed=42)
        fb = session.trials[0]['inputs'][0, 6:8].tolist()
        assert fb == [0.0, 0.0], f"First trial feedback should be [0,0], got {fb}"

    def test_feedback_subsequent_trials(self):
        """Subsequent trials should have feedback [1, 0] (assumed correct)."""
        session = Task1Session(n_trials=10, seed=42)
        for trial in session.trials[1:]:
            fb = trial['inputs'][0, 6:8].tolist()
            assert fb == [1.0, 0.0], (
                f"Subsequent trial feedback should be [1,0], got {fb}"
            )


class TestTask2Session:
    """Tests for Task2Session trial generator."""

    def test_input_shape(self):
        """Each trial input tensor should be [20, 8]."""
        session = Task2Session(n_trials=50, seed=42)
        for trial in session.trials:
            assert trial['inputs'].shape == (20, 8), (
                f"Expected shape (20, 8), got {trial['inputs'].shape}"
            )

    def test_congruent_trials_matching_directions(self):
        """Congruent trials must have matching visual and auditory directions."""
        session = Task2Session(n_trials=300, seed=42)
        for trial in session.trials:
            meta = trial['metadata']
            if meta['is_congruent']:
                assert meta['visual_direction'] == meta['auditory_direction'], (
                    f"Congruent trial has mismatched directions: "
                    f"visual={meta['visual_direction']}, "
                    f"auditory={meta['auditory_direction']}"
                )

    def test_incongruent_correct_matches_relevant(self):
        """On incongruent trials, correct answer matches the relevant modality."""
        session = Task2Session(n_trials=500, n_switches=5, seed=42)
        for trial in session.trials:
            meta = trial['metadata']
            if not meta['is_congruent']:
                # Correct direction should match the relevant modality
                if meta['context'] == 'visual':
                    assert meta['correct_direction'] == meta['visual_direction'], (
                        f"In visual context, correct should match visual direction"
                    )
                else:
                    assert meta['correct_direction'] == meta['auditory_direction'], (
                        f"In auditory context, correct should match auditory direction"
                    )

    def test_congruent_proportion(self):
        """Proportion of congruent trials should approximate congruent_ratio.

        With congruent_ratio=0.5 (default), roughly half the trials should
        be congruent. With a large enough sample, this is reliable.
        """
        session = Task2Session(n_trials=1000, seed=42)
        n_congruent = sum(
            1 for t in session.trials if t['metadata']['is_congruent']
        )
        proportion = n_congruent / len(session.trials)
        assert abs(proportion - session.congruent_ratio) < 0.06, (
            f"Congruent proportion = {proportion:.3f}, expected ~{session.congruent_ratio}"
        )

    def test_context_cue_always_zero(self):
        """Task 2 context cue channels should always be [0, 0]."""
        session = Task2Session(n_trials=300, seed=42)
        for trial in session.trials:
            # Check all timesteps
            ctx_cue = trial['inputs'][:, 4:6]
            assert (ctx_cue == 0).all(), (
                "Task 2 context cue should be [0, 0] at all timesteps"
            )

    def test_both_modalities_informative(self):
        """Both modalities should carry non-zero stimuli during stimulus epoch."""
        session = Task2Session(n_trials=100, seed=42)
        for trial in session.trials:
            stim = trial['inputs'][5]  # stimulus epoch
            vis = stim[0:2]
            aud = stim[2:4]
            assert vis.abs().sum() > 0, "Visual stimulus should be non-zero"
            assert aud.abs().sum() > 0, "Auditory stimulus should be non-zero"

    def test_incongruent_directions_differ(self):
        """Incongruent trials must have different visual and auditory directions."""
        session = Task2Session(n_trials=300, seed=42)
        for trial in session.trials:
            meta = trial['metadata']
            if not meta['is_congruent']:
                assert meta['visual_direction'] != meta['auditory_direction'], (
                    f"Incongruent trial should have different directions"
                )
