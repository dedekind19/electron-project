"""
Monte Carlo simulation of relativistic electron energy losses
in a magnetized plasma.
"""

from time import time

import numpy as np
from plasma_sim.constants import SIGMA_T, ALPHA_F, C

from plasma_sim.losses import (
    synchrotron_loss,
    inverse_compton_loss,
    bremsstrahlung_loss,
    coulomb_loss,
)   


def sample_gamma_powerlaw(
    gamma_min: float,
    gamma_max: float,
    spectral_index: float,
    n_electrons: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Sample initial Lorentz factors from a power law distribution.

    Uses the inverse transform method to sample from N(γ) ∝ γ^(-p).

    Parameters
    ----------
    gamma_min : float
        Minimum Lorentz factor (dimensionless)
    gamma_max : float
        Maximum Lorentz factor (dimensionless)
    spectral_index : float
        Power law index p, must be > 1
    n_electrons : int
        Number of electrons to sample
    rng : np.random.Generator
        Random number generator instance for reproducibility

    Returns
    -------
    np.ndarray
        Array of n_electrons Lorentz factors sampled from the power law
    """
    r = rng.uniform(0, 1, n_electrons)
    exponent = 1 - spectral_index
    gamma = (
        r * (gamma_max**exponent - gamma_min**exponent) + gamma_min**exponent
    ) ** (1 / exponent)
    return gamma

def compute_interaction_probability(
    gamma: float,
    n_plasma: float,
    dt: float,
) -> float:
    """Compute the probability of a Coulomb collision in a timestep dt.

    Uses the Poisson process formula p = 1 - exp(-dt / t_collision)
    where t_collision = 1 / (n_plasma * σ_coulomb(γ) * c)
    and σ_coulomb(γ) = σ_T * ln(Λ) / (8 * α_f * γ)

    Parameters
    ----------
    gamma : float
        Lorentz factor of the electron (dimensionless)
    n_plasma : float
        Plasma electron number density (m^-3)
    dt : float
        Timestep (s)

    Returns
    -------
    float
        Probability of a Coulomb collision in dt, between 0 and 1
    """
    if n_plasma == 0.0:
        return 0.0
    COULOMB_LOG = 30.0
    sigma_coulomb = SIGMA_T * COULOMB_LOG / (8 * ALPHA_F * gamma)
    t_collision = 1 / (n_plasma * sigma_coulomb * C)
    return 1 - np.exp(-dt / t_collision)

###############################
#MAGNETIC FIELD VARIATION
###############################

def sample_new_B(
    B_initial: float,
    sigma_B: float,
    rng: np.random.Generator,
) -> float:
    """ Models turbulent fluctuations in the lobe magnetic field by sampling
    from a Gaussian centred on B_initial. If sigma_B is zero, returns
    B_initial exactly. Negative values are rejected by resampling.

    Parameters
    ----------
    B_initial : float
        Central value of the magnetic field (T)
    sigma_B : float
        Standard deviation of the Gaussian (T)
    rng : np.random.Generator
        Random number generator instance for reproducibility

    Returns
    -------
    float
        New magnetic field value (T), always positive
    """
    if sigma_B == 0.0:
        return B_initial
    B = rng.normal(B_initial, sigma_B)
    while B <= 0:
        B = rng.normal(B_initial, sigma_B)
    return B

###################################
#MAIN SIMULATION FUNCTION
###################################

def run_simulation(config: dict) -> dict:
    """Run the Monte Carlo simulation of electron energy losses.

    Evolves a population of electrons through synchrotron, inverse Compton,
    bremsstrahlung and Coulomb losses until all electrons have cooled below
    gamma_min or the synchrotron power drops below epsilon_stop * P_initial.

    Parameters
    ----------
    config : dict
        Configuration dictionary with 'physical' and 'numerical' sections
        as described in the README.

    Returns
    -------
    dict
        Time series of population-averaged quantities with keys:
        time, n_alive, gamma_mean, dEdt_sync_mean, dEdt_IC_mean,
        dEdt_brems_mean, dEdt_coulomb_mean
    """
    phys = config["physical"]
    num = config["numerical"]

    B_field = phys["B_field"]
    sigma_B = phys["sigma_B"]
    n_plasma = phys["n_plasma"]
    gamma_min = phys["gamma_min"]
    epsilon_stop = phys["epsilon_stop"]
    redshift = phys["redshift"]
    dt = num["dt"]
    n_electrons = num["n_electrons"]
    B_update_steps = num["B_update_steps"]

    rng = np.random.default_rng(num["random_seed"])

    # sample initial gamma values from power law
    gammas = sample_gamma_powerlaw(
        gamma_min=phys["gamma_min_init"],
        gamma_max=phys["gamma_max_init"],
        spectral_index=phys["spectral_index"],
        n_electrons=n_electrons,
        rng=rng
    )

    # compute initial synchrotron power for stopping condition
    P_initial = np.abs(synchrotron_loss(gammas, B_field))

    # track which electrons are still alive
    alive = np.ones(n_electrons, dtype=bool)

    # storage for time series
    time_series = []
    n_alive_series = []
    gamma_mean_series = []
    dEdt_sync_series = []
    dEdt_IC_series = []
    dEdt_brems_series = []
    dEdt_coulomb_series = []

    B = B_field
    step = 0

    while np.any(alive):
        import time
        start = time.time()
        # update B every B_update_steps
        if step % B_update_steps == 0:
            B = sample_new_B(B_field, sigma_B, rng)

        # compute losses for all electrons at once (vectorized)
        dEdt_sync = synchrotron_loss(gammas, B)
        dEdt_IC = inverse_compton_loss(gammas, redshift)
        dEdt_brems = bremsstrahlung_loss(gammas, n_plasma)

        # Coulomb - stochastic, vectorized
        dEdt_coulomb = np.zeros(n_electrons)
        p = compute_interaction_probability(gammas, n_plasma, dt)
        collisions = rng.uniform(size=n_electrons) < p
        dEdt_coulomb[collisions & alive] = coulomb_loss(
            gammas[collisions & alive], n_plasma
        )

        # update gammas for living electrons only
        gammas[alive] += (
            dEdt_sync[alive] + dEdt_IC[alive] +
            dEdt_brems[alive] + dEdt_coulomb[alive]
        ) * dt

        # check stopping conditions vectorized
        P_sync = np.abs(synchrotron_loss(gammas, B))
        alive &= (gammas >= gamma_min) & (P_sync >= epsilon_stop * P_initial)

        # record time series for living electrons
        n_alive = int(np.sum(alive))
        if n_alive > 0:
            time_series.append(step * dt)
            n_alive_series.append(n_alive)
            gamma_mean_series.append(float(np.mean(gammas[alive])))
            dEdt_sync_series.append(float(np.mean(dEdt_sync[alive])))
            dEdt_IC_series.append(float(np.mean(dEdt_IC[alive])))
            dEdt_brems_series.append(float(np.mean(dEdt_brems[alive])))
            dEdt_coulomb_series.append(float(np.mean(dEdt_coulomb[alive])))

        step += 1

    return {
        "time": time_series,
        "n_alive": n_alive_series,
        "gamma_mean": gamma_mean_series,
        "dEdt_sync_mean": dEdt_sync_series,
        "dEdt_IC_mean": dEdt_IC_series,
        "dEdt_brems_mean": dEdt_brems_series,
        "dEdt_coulomb_mean": dEdt_coulomb_series,
    }
    print(f"Steps: {step}, Time: {time.time()-start:.2f}s")