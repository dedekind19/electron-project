# Electron Monte Carlo Simulation in Magnetized Plasma

Monte Carlo simulation of relativistic electron energy losses in a
magnetized plasma.

## Table of Contents
- [Physical Model](#physical-model)
- [Installation](#installation)
- [Usage](#usage)
- [Configuration](#configuration)
- [Output](#output)
- [Testing](#testing)
- [References](#references)

---

## Physical Model

### Assumptions

- Electrons are **ultra-relativistic**, described by Lorentz factor γ >> 1
- Magnetic field B is initialised by the user and updated every
  `B_update_steps` timesteps by sampling from a Gaussian distribution
  centred on the initial value with s.d. `sigma_B`. This
  accounts for turbulent fluctuations.
- **Boundary conditions are neglected**: the simulation tracks energy
  evolution in time only, not spatial motion. This is justified when
  the cooling length is much smaller than the size of the system
- Inverse Compton scattering is treated in the **Thomson regime**, valid
  for γ << 10⁸. 
- The plasma is fully ionized hydrogen (Z=1)
- The initial Lorentz factors of the electrons are sampled from a power
  law distribution N(γ) ∝ γ^(-p) with user-defined γ_min and γ_max.
  For example, in radio galaxy lobes, typical values are γ_min ~ 100 and
  γ_max ~ 10⁶.

### Physical context and parameter guidance

The simulation is general. Users should set parameters appropriate to their physical context. As a reference:

**Radio galaxy lobes:**
- B_field ~ 1e-10 to 3e-10 T (1-3 μG)
- n_plasma ~ 1e-1 m⁻³ (10⁻⁴ cm⁻³)
- In this regime, Bremsstrahlung and Coulomb losses are negligible
  compared to synchrotron and inverse Compton

**Dense laboratory plasma:**
- B_field ~ 1e-2 T
- n_plasma ~ 1e18 m⁻³
- In this regime, Coulomb losses dominate at low γ

### Energy Loss Processes

At each timestep dt, the Lorentz factor γ is updated by summing
contributions from all four processes:
dγ/dt = (dγ/dt)_sync + (dγ/dt)_IC + (dγ/dt)_brems + (dγ/dt)_coulomb

#### 1. Synchrotron Radiation

An electron gyrating in magnetic field B radiates energy continuously:
(dγ/dt)_sync = - (4/3) * (σ_T * c) / (m_e * c²) * γ² * U_B
U_B = B² / (2 * μ₀)

| Constant | Value | Description |

| σ_T | 6.6524 × 10⁻²⁹ m² | Thomson cross-section |
| c | 2.9979 × 10⁸ m/s | Speed of light |
| m_e | 9.1094 × 10⁻³¹ kg | Electron mass |
| μ₀ | 1.2566 × 10⁻⁶ H/m | Vacuum permeability |

#### 2. Inverse Compton Scattering

Electrons scatter off CMB photons. In the Thomson regime the loss rate
has the same form as synchrotron with U_B replaced by U_rad:
(dγ/dt)_IC = - (4/3) * (σ_T * c) / (m_e * c²) * γ² * U_rad
U_rad = a_rad * T_CMB⁴ * (1 + z)⁴

| Constant | Value | Description |

| a_rad | 7.5657 × 10⁻¹⁶ J m⁻³ K⁻⁴ | Radiation constant |
| T_CMB | 2.7255 K | CMB temperature at z=0 |
| U_rad(z=0) | 4.1719 × 10⁻¹⁴ J/m³ | CMB energy density at z=0 |

The (1+z)⁴ factor accounts for the redshift dependence of U_rad.

#### 3. Bremsstrahlung

Electrons decelerate in the Coulomb field of plasma ions. For a fully
ionized hydrogen plasma (Z=1) in the relativistic regime:
(dγ/dt)_brems = - n_e * c * σ_T * α_f * γ * 14.3

where 14.3 ≈ ln(183) + 1/18 is the numerical factor for Z=1.

| Constant | Value | Description |

| α_f | 7.2974 × 10⁻³ (≈ 1/137) | Fine structure constant |
| n_e | user parameter (m⁻³) | Plasma electron number density |

#### 4. Coulomb Collisions

Inelastic collisions with thermal plasma particles, dominant at low γ:
(dγ/dt)_coulomb = - n_e * c * σ_T / (8 * α_f) * (1/γ) * ln(Λ)

where ln(Λ) ≈ 30 is the Coulomb logarithm

### Derived Quantities

#### Cooling Time
t_cool = γ / |dγ/dt_total|

This is the time for an electron to lose a factor of 1/e of its energy
to all processes combined. It is a function of γ and is tracked as the
population evolves.

### Monte Carlo Method

Synchrotron, inverse Compton and Bremsstrahlung are treated as
continuous losses and applied at every timestep. Coulomb collisions
are treated as discrete stochastic events: at each timestep, a random
number is drawn and compared to the interaction probability:
p = 1 - exp(-dt / t_collision)
t_collision = 1 / (n_plasma * σ_coulomb(γ) * c)
σ_coulomb(γ) = σ_T * ln(Λ) / (8 * α_f * γ)

If the random number is less than p, the collision occurs and the
Coulomb energy loss is applied. The cross section is recomputed at
every step from the current γ.

The initial Lorentz factors are sampled from a power law distribution
using the inverse transform method. The magnetic field is updated every
`B_update_steps` timesteps by sampling from a Gaussian centred on the
initial value.

Each electron evolves with an adaptive timestep computed from its
current cooling time:
dt = 0.01 * γ / |dγ/dt_total|


---

## Installation

Requires Python >= 3.10.

```bash
git clone https://github.com/dedekind19/electron-project.git
cd electron-project
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -e ".[dev]"
```

---

## Usage

### Run a simulation

```bash
python -m plasma_sim run --config configs/default.json --output results/
```

---

## Configuration

All parameters are split into two sections: `physical` for the
plasma system and `numerical` for the simulation itself. The default parameters are the ones found in radio galxy lobes in condition of equipartition (U_B=U_rad). They can be cahnged in the following json:


Example `configs/default.json`:

```json
{
  "physical": {
    "B_field":         3e-10,
    "sigma_B":         3e-11,
    "n_plasma":        1e2,
    "gamma_min_init":  100,
    "gamma_max_init":  1e6,
    "spectral_index":  2.0,
    "gamma_min":       10.0,
    "epsilon_stop":    0.01,
    "redshift":        0.0
  },
  "numerical": {
    "n_electrons":     1000,
    "n_bins":          50,
    "random_seed":     42,
    "B_update_steps":  100
  }
}
```



### Physical parameters

| Parameter | Description | Default | Units | Valid range |

| B_field | Magnetic field strength | 3e-10 | T | > 0 |
| sigma_B | Standard deviation of B fluctuations | 3e-11 | T | ≥ 0 |
| n_plasma | Plasma electron number density | 1e2 | m⁻³ | > 0 |
| gamma_min_init | Minimum γ for initial power law sampling | 100 | — | ≥ 1 |
| gamma_max_init | Maximum γ for initial power law sampling | 1e6 | — | > gamma_min_init |
| spectral_index | Power law index p for initial γ distribution | 2.0 | — | > 1 |
| gamma_min | Minimum γ: electron stops below this | 10.0 | — | ≥ 1 |
| epsilon_stop | Synchrotron power cutoff fraction | 0.01 | — | 0 < ε < 1 |
| redshift | Source redshift (scales U_rad) | 0.0 | — | ≥ 0 |

### Numerical parameters

| Parameter | Description | Default | Units | Valid range |

| n_electrons | Number of simulated electrons | 1000 | — | ≥ 1 |
| n_bins | Number of bins for all histograms | 50 | — | ≥ 10 |
| random_seed | Random seed for reproducibility | 42 | — | any integer |
| B_update_steps | Number of steps between B updates | 100 | — | ≥ 1 |

### Validation and warnings

| Condition | Type | Reason |

| B_field ≤ 0 or n_plasma ≤ 0 | Error | Unphysical |
| sigma_B < 0 | Error | Unphysical |
| gamma_min_init ≥ gamma_max_init | Error | Empty sampling range |
| spectral_index ≤ 1 | Error | Non-normalisable power law |
| gamma_min_init < gamma_min | Warning | Electrons may stop immediately |
| gamma_max_init > 1e7 | Warning | Thomson approximation becoming inaccurate |
| redshift < 0 | Error | Unphysical |
| epsilon_stop ≤ 0 or ≥ 1 | Error | Unphysical stopping condition |

---

## Output

The `run` command produces a single `results.json` in the output
directory with the following structure:

```json
{
  "config": { ... },

  "time_series": {
    "time":              [...],
    "n_alive":           [...],
    "gamma_mean":        [...],
    "dEdt_sync_mean":    [...],
    "dEdt_IC_mean":      [...],
    "dEdt_brems_mean":   [...],
    "dEdt_coulomb_mean": [...]
  },

  "losses_vs_gamma": {
    "gamma_bins":        [...],
    "dEdgamma_sync":     [...],
    "dEdgamma_IC":       [...],
    "dEdgamma_brems":    [...],
    "dEdgamma_coulomb":  [...]
  },

  "summary": {
    "t_cool_mean":       ...,
    "t_cool_std":        ...,
    "gamma_final_mean":  ...,
    "gamma_final_std":   ...,
    "frac_sync_mean":    ...,
    "frac_IC_mean":      ...,
    "frac_brems_mean":   ...,
    "frac_coulomb_mean": ...
  }
}
```
## Plotting
The following command is used to generate some plot derived from the quantities in the Json. These are just exemplary, the code that generates them is in 'plot.py' The file with all the data is stored in the same folder:
```bash
python -m plasma_sim plot --input results/results.json --output results/



---

## Testing

```bash
pytest tests/ --cov=src/ --cov-report=term-missing
```

---

## References

- Rybicki & Lightman, *Radiative Processes in Astrophysics*, Wiley (1979)
- Blumenthal & Gould, Rev. Mod. Phys. **42**, 237 (1970)
- Longair, M.S., *High Energy Astrophysics*