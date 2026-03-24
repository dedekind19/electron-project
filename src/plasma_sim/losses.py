"""
Energy loss functions for relativistic electrons in a magnetized plasma.

Each function returns dγ/dt for a single loss process, given the current
Lorentz factor and the relevant physical parameters.
"""

from plasma_sim.constants import SIGMA_T, C, M_E, MU_0


def synchrotron_loss(gamma: float, B: float) -> float:
    """Compute the synchrotron energy loss rate dγ/dt.

    dγ/dt = -(4/3) * (σ_T * c) / (m_e * c²) * γ² * U_B
    where U_B = B² / (2 * μ₀)

    Parameters
    ----------
    gamma : float
        Lorentz factor of the electron (dimensionless)
    B : float
        Magnetic field strength (T)

    Returns
    -------
    float
        dγ/dt due to synchrotron radiation (s^-1), always <= 0
    """
    U_B = B**2 / (2 * MU_0)
    return -(4 / 3) * (SIGMA_T * C) / (M_E * C**2) * gamma**2 * U_B