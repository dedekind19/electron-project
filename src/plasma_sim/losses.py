"""
Energy loss functions for relativistic electrons in a magnetized plasma.

Each function returns dγ/dt for a single loss process, given the current
Lorentz factor and the relevant physical parameters.
"""

from plasma_sim.constants import SIGMA_T, C, M_E, MU_0, U_RAD_0, ALPHA_F


"SYNCHROTRON"

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


"INVERSE COMPTON"

def inverse_compton_loss(gamma: float, redshift: float) -> float:
    """Compute the inverse Compton energy loss rate dγ/dt.

    dγ/dt = -(4/3) * (σ_T * c) / (m_e * c²) * γ² * U_rad
    where U_rad = U_RAD_0 * (1 + z)^4

    Parameters
    ----------
    gamma : float
        Lorentz factor of the electron (dimensionless)
    redshift : float
        Redshift of the source (dimensionless), scales the CMB energy density

    Returns
    -------
    float
        dγ/dt due to inverse Compton scattering (s^-1), always <= 0
    """
    U_rad = U_RAD_0 * (1 + redshift)**4
    return -(4 / 3) * (SIGMA_T * C) / (M_E * C**2) * gamma**2 * U_rad

"BREMMSTRAHLUNG"
def bremsstrahlung_loss(gamma: float, n_plasma: float) -> float:
    """Compute the bremsstrahlung energy loss rate dγ/dt.

    dγ/dt = -n_e * c * σ_T * α_f * γ * 14.3
    where 14.3 ≈ ln(183) + 1/18 is the numerical factor for Z=1.

    Parameters
    ----------
    gamma : float
        Lorentz factor of the electron (dimensionless)
    n_plasma : float
        Plasma electron number density (m^-3)

    Returns
    -------
    float
        dγ/dt due to bremsstrahlung (s^-1), always <= 0
    """
    return -n_plasma * C * SIGMA_T * ALPHA_F * gamma * 14.3