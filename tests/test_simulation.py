"""
Tests for simulation functions in plasma_sim.simulation.
"""
import numpy as np
from plasma_sim.simulation import sample_gamma_powerlaw, compute_interaction_probability, sample_new_B, run_simulation


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


###########################
#MAIN SIMULATION TESTS
###########################

def test_run_simulation_time_increases_monotonically():
    """Test that simulation time increases by dt at each step.

    GIVEN: a valid configuration
    WHEN: the simulation is run
    THEN: time must increase monotonically by dt at each step
    """
    config = {
        "physical": {
            "B_field": 1e-9, "sigma_B": 0.0, "n_plasma": 1e3,
            "gamma_min_init": 1e3, "gamma_max_init": 1e5,
            "spectral_index": 2.0, "gamma_min": 10.0,
            "epsilon_stop": 0.01, "t_snap": 1e15, "redshift": 0.0
        },
        "numerical": {
            "dt": 1e13, "n_electrons": 10, "random_seed": 42,
            "B_update_steps": 100, "n_bins": 50
        }
    }
    results = run_simulation(config)
    time = results["time_series"]["time"]
    diffs = [time[i+1] - time[i] for i in range(len(time)-1)]
    assert all(d > 0 for d in diffs)


def test_run_simulation_n_alive_only_decreases():
    """Test that the number of living electrons never increases.

    GIVEN: a valid configuration
    WHEN: the simulation is run
    THEN: n_alive must be non-increasing at every step
    """
    config = {
        "physical": {
            "B_field": 1e-9, "sigma_B": 0.0, "n_plasma": 1e3,
            "gamma_min_init": 1e3, "gamma_max_init": 1e5,
            "spectral_index": 2.0, "gamma_min": 10.0,
            "epsilon_stop": 0.01, "t_snap": 1e15, "redshift": 0.0
        },
        "numerical": {
            "dt": 1e13, "n_electrons": 10, "random_seed": 42,
            "B_update_steps": 100, "n_bins": 50
        }
    }
    results = run_simulation(config)
    n_alive = results["time_series"]["n_alive"]
    assert all(n_alive[i+1] <= n_alive[i] for i in range(len(n_alive)-1))


def test_run_simulation_gamma_mean_only_decreases():
    """Test that the mean Lorentz factor never increases.

    GIVEN: a valid configuration
    WHEN: the simulation is run
    THEN: gamma_mean must be non-increasing at every step
    """
    config = {
        "physical": {
            "B_field": 1e-9, "sigma_B": 0.0, "n_plasma": 1e3,
            "gamma_min_init": 1e3, "gamma_max_init": 1e5,
            "spectral_index": 2.0, "gamma_min": 10.0,
            "epsilon_stop": 0.01, "t_snap": 1e15, "redshift": 0.0
        },
        "numerical": {
            "dt": 1e13, "n_electrons": 10, "random_seed": 42,
            "B_update_steps": 100, "n_bins": 50
        }
    }
    results = run_simulation(config)
    gamma_mean = results["time_series"]["gamma_mean"]
    assert all(gamma_mean[i+1] <= gamma_mean[i] for i in range(len(gamma_mean)-1))


def test_run_simulation_loss_rates_are_negative():
    """Test that all energy loss rates are negative at every step.

    GIVEN: a valid configuration
    WHEN: the simulation is run
    THEN: all loss rates must be negative since electrons only lose energy
    """
    config = {
        "physical": {
            "B_field": 1e-9, "sigma_B": 0.0, "n_plasma": 1e3,
            "gamma_min_init": 1e3, "gamma_max_init": 1e5,
            "spectral_index": 2.0, "gamma_min": 10.0,
            "epsilon_stop": 0.01, "t_snap": 1e15, "redshift": 0.0
        },
        "numerical": {
            "dt": 1e13, "n_electrons": 10, "random_seed": 42,
            "B_update_steps": 100, "n_bins": 50
        }
    }
    results = run_simulation(config)
    assert all(v <= 0 for v in results["time_series"]["dEdt_sync_mean"])
    assert all(v <= 0 for v in results["time_series"]["dEdt_IC_mean"])
    assert all(v <= 0 for v in results["time_series"]["dEdt_brems_mean"])
    assert all(v <= 0 for v in results["time_series"]["dEdt_coulomb_mean"])