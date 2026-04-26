"""
Tests for simulation functions in plasma_sim.simulation.
"""
import numpy as np
from plasma_sim.simulation import sample_gamma_powerlaw, compute_interaction_probability, sample_new_B


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


   

def test_interaction_probability_between_zero_and_one():
    """Test that the interaction probability is always a valid probability.

    GIVEN: typical plasma density, Lorentz factor and timestep
    WHEN: the interaction probability is computed
    THEN: the result must be between 0 and 1
    """
    p = compute_interaction_probability(gamma=1e4, n_plasma=1e3, dt=1e11)
    assert 0 <= p <= 1


def test_interaction_probability_zero_for_zero_density():
    """Test that probability vanishes when there is no plasma.

    GIVEN: zero plasma density
    WHEN: the interaction probability is computed
    THEN: the result must be zero since there are no ions to collide with
    """
    p = compute_interaction_probability(gamma=1e4, n_plasma=0.0, dt=1e11)
    assert p == 0.0


def test_interaction_probability_increases_with_density():
    """Test that probability increases with plasma density.

    GIVEN: two identical electrons in plasmas with different densities
    WHEN: the interaction probability is computed for both
    THEN: higher density must give higher probability
    """
    p1 = compute_interaction_probability(gamma=1e4, n_plasma=1e3, dt=1e11)
    p2 = compute_interaction_probability(gamma=1e4, n_plasma=2e3, dt=1e11)
    assert p2 > p1


def test_interaction_probability_increases_with_dt():
    """Test that probability increases with timestep.

    GIVEN: two identical electrons with different timesteps
    WHEN: the interaction probability is computed for both
    THEN: longer timestep must give higher probability
    """
    p1 = compute_interaction_probability(gamma=1e4, n_plasma=1e3, dt=1e11)
    p2 = compute_interaction_probability(gamma=1e4, n_plasma=1e3, dt=2e11)
    assert p2 > p1


def test_interaction_probability_decreases_with_gamma():
    """Test that probability decreases with Lorentz factor.

    GIVEN: two electrons with different Lorentz factors
    WHEN: the interaction probability is computed for both
    THEN: higher gamma must give lower probability since σ_coulomb ∝ 1/γ
    """
    p1 = compute_interaction_probability(gamma=1e4, n_plasma=1e3, dt=1e11)
    p2 = compute_interaction_probability(gamma=2e4, n_plasma=1e3, dt=1e11)
    assert p2 < p1

##########################
#MAGNETIC FIELD VARIATIONS
##########################

def test_sample_new_B_zero_sigma_returns_initial():
    """Test that zero sigma always returns the initial B value.

    GIVEN: an initial B value and sigma_B = 0
    WHEN: a new B is sampled
    THEN: the result must equal B_initial exactly
    """
    rng = np.random.default_rng(42)
    B = sample_new_B(B_initial=1e-9, sigma_B=0.0, rng=rng)
    assert B == 1e-9


def test_sample_new_B_is_positive():
    """Test that sampled B values are always positive.

    GIVEN: a typical initial B and sigma_B
    WHEN: many B values are sampled
    THEN: all results must be positive since negative B is unphysical
    """
    rng = np.random.default_rng(42)
    samples = [sample_new_B(B_initial=1e-9, sigma_B=1e-10, rng=rng)
               for _ in range(1000)]
    assert all(b > 0 for b in samples)


def test_sample_new_B_reproducible_with_seed():
    """Test that the same seed always produces the same B value.

    GIVEN: two RNG instances initialised with the same seed
    WHEN: a new B is sampled from both
    THEN: the two results must be identical
    """
    rng1 = np.random.default_rng(42)
    rng2 = np.random.default_rng(42)
    B1 = sample_new_B(B_initial=1e-9, sigma_B=1e-10, rng=rng1)
    B2 = sample_new_B(B_initial=1e-9, sigma_B=1e-10, rng=rng2)
    assert B1 == B2


def test_sample_new_B_mean_close_to_initial():
    """Test that the mean of many samples is close to B_initial.

    GIVEN: a Gaussian centred on B_initial
    WHEN: many B values are sampled
    THEN: the sample mean must be close to B_initial
    """
    rng = np.random.default_rng(42)
    samples = [sample_new_B(B_initial=1e-9, sigma_B=1e-10, rng=rng)
               for _ in range(10000)]
    assert np.isclose(np.mean(samples), 1e-9, rtol=0.05)