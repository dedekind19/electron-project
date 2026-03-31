"""
Monte Carlo simulation of relativistic electron energy losses
in a magnetized plasma.
"""

import numpy as np
from plasma_sim.constants import SIGMA_T, ALPHA_F, C


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