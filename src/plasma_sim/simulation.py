"""
Monte Carlo simulation of relativistic electron energy losses
in a magnetized plasma.
"""

import numpy as np


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