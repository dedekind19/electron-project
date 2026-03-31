"""
Tests for simulation functions in plasma_sim.simulation.
"""

import numpy as np
import pytest
from plasma_sim.simulation import sample_gamma_powerlaw


def test_sample_gamma_powerlaw_correct_number_of_samples():
    """Test that the function returns the correct number of samples.

    GIVEN: a power law distribution with standard parameters
    WHEN: n_electrons samples are requested
    THEN: the output array has exactly n_electrons elements
    """
    rng = np.random.default_rng(42)
    samples = sample_gamma_powerlaw(
        gamma_min=100, gamma_max=1e6, spectral_index=2.0,
        n_electrons=1000, rng=rng
    )
    assert len(samples) == 1000


def test_sample_gamma_powerlaw_within_bounds():
    """Test that all sampled gamma values are within the specified range.

    GIVEN: a power law distribution with gamma_min=100 and gamma_max=1e6
    WHEN: samples are drawn
    THEN: all values must lie within [gamma_min, gamma_max]
    """
    rng = np.random.default_rng(42)
    samples = sample_gamma_powerlaw(
        gamma_min=100, gamma_max=1e6, spectral_index=2.0,
        n_electrons=1000, rng=rng
    )
    assert np.all(samples >= 100)
    assert np.all(samples <= 1e6)


def test_sample_gamma_powerlaw_reproducible_with_seed():
    """Test that the same seed always produces the same samples.

    GIVEN: two RNG instances initialised with the same seed
    WHEN: samples are drawn from both
    THEN: the two sample arrays must be identical
    """
    rng1 = np.random.default_rng(42)
    rng2 = np.random.default_rng(42)
    samples1 = sample_gamma_powerlaw(
        gamma_min=100, gamma_max=1e6, spectral_index=2.0,
        n_electrons=1000, rng=rng1
    )
    samples2 = sample_gamma_powerlaw(
        gamma_min=100, gamma_max=1e6, spectral_index=2.0,
        n_electrons=1000, rng=rng2
    )
    assert np.allclose(samples1, samples2)


def test_sample_gamma_powerlaw_different_seeds_differ():
    """Test that different seeds produce different samples.

    GIVEN: two RNG instances initialised with different seeds
    WHEN: samples are drawn from both
    THEN: the two sample arrays must not be identical
    """
    rng1 = np.random.default_rng(42)
    rng2 = np.random.default_rng(99)
    samples1 = sample_gamma_powerlaw(
        gamma_min=100, gamma_max=1e6, spectral_index=2.0,
        n_electrons=1000, rng=rng1
    )
    samples2 = sample_gamma_powerlaw(
        gamma_min=100, gamma_max=1e6, spectral_index=2.0,
        n_electrons=1000, rng=rng2
    )
    assert not np.allclose(samples1, samples2)