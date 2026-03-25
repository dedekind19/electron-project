"""
Physical constants used in the plasma simulation.

All values are in SI units and taken from NIST CODATA 2018:
https://physics.nist.gov/cuu/Constants/index.html
"""

# Speed of light (m/s)
C = 2.99792458e8

# Electron mass (kg)
M_E = 9.1093837015e-31

# Thomson cross-section (m^2)
SIGMA_T = 6.6524587158e-29

# Vacuum permeability (H/m)
MU_0 = 1.25663706212e-6

# Fine structure constant (dimensionless)
ALPHA_F = 7.2973525693e-3

# Radiation constant a_rad (J m^-3 K^-4)
A_RAD = 7.5657e-16

# CMB temperature at z=0 (K)
T_CMB = 2.7255

# CMB energy density at z=0 (J/m^3)
# Computed as A_RAD * T_CMB^4
U_RAD_0 = A_RAD * T_CMB**4