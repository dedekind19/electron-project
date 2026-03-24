# Electron Monte Carlo Simulation in Magnetized Plasma

Monte Carlo simulation of relativistic electron energy losses in a
magnetized plasma, motivated by the physical conditions found in radio
galaxy lobes where synchrotron radiation is produced.

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
- Magnetic field **B** is uniform and constant, justified for large-scale
  lobes where B ~ 1-100 μG over scales of tens of kpc
- **Boundary conditions are neglected**: the cooling length λ_cool for all
  processes is much smaller than the typical lobe size (~10-100 kpc), so
  electrons lose all their energy before reaching the boundary
- **The simulation tracks energy evolution in time, not spatial motion**
- Inverse Compton scattering is treated in the **Thomson regime**, valid
  for γ << 10⁸. For typical radio galaxy lobe electrons (γ ~ 10³-10⁶)
  this condition is satisfied by several orders of magnitude, making the
  Klein-Nishina correction negligible
- The plasma is fully ionized hydrogen (Z=1)

### Energy Loss Processes

At each timestep dt, the Lorentz factor γ is updated by summing
At each timestep dt, the Lorentz factor γ is updated according to:
```
dγ/dt = (dγ/dt)_sync + (dγ/dt)_IC + (dγ/dt)_brems + (dγ/dt)_coulomb
```
dγ/dt = (dγ/dt)_sync + (dγ/dt)_IC + (dγ/dt)_brems + (dγ/dt)_coulomb

#### 1. Synchrotron Radiation

An electron gyrating in magnetic field B radiates energy continuously:
```
(dγ/dt)_sync = - (4/3) * (σ_T * c) / (m_e * c²) * γ² * U_B

U_B = B² / (2 * μ₀)
```

| Constant | Value | Description |
|---|---|---|
| σ_T | 6.6524 × 10⁻²⁹ m² | Thomson cross-section |
| c | 2.9979 × 10⁸ m/s | Speed of light |
| m_e | 9.1094 × 10⁻³¹ kg | Electron mass |
| μ₀ | 1.2566 × 10⁻⁶ H/m | Vacuum permeability |


#### 2. Inverse Compton Scattering

Electrons scatter off CMB photons. In the Thomson regime the loss rate
has the same form as synchrotron with U_B replaced by U_rad:
```
(dγ/dt)_IC = - (4/3) * (σ_T * c) / (m_e * c²) * γ² * U_rad

U_rad = a_rad * T_CMB⁴ * (1 + z)⁴
```

| Constant | Value | Description |
|---|---|---|
| a_rad | 7.5657 × 10⁻¹⁶ J m⁻³ K⁻⁴ | Radiation constant |
| T_CMB | 2.7255 K | CMB temperature at z=0 |
| U_rad(z=0) | 4.1719 × 10⁻¹⁴ J/m³ | CMB energy density at z=0 |

The (1+z)⁴ factor accounts for the redshift dependence of U_rad.


#### 3. Bremsstrahlung

Electrons decelerate in the Coulomb field of plasma ions. For a fully
ionized hydrogen plasma (Z=1) in the relativistic regime:
```
(dγ/dt)_brems = - n_e * c * σ_T * α_f * γ * 14.3
```

where 14.3 ≈ ln(183) + 1/18 is the numerical factor for Z=1.

| Constant | Value | Description |
|---|---|---|
| α_f | 7.2974 × 10⁻³ (≈ 1/137) | Fine structure constant |
| n_e | user parameter (m⁻³) | Plasma electron number density |


#### 4. Coulomb Collisions

Inelastic collisions with thermal plasma particles, dominant at low γ:
```
(dγ/dt)_coulomb = - n_e * c * σ_T / (8 * α_f) * (1/γ) * ln(Λ)
```

where ln(Λ) ≈ 30 is the Coulomb logarithm, weakly dependent on plasma
parameters and typical for radio galaxy lobe conditions.

### Derived Quantities

#### Cooling Length and Cooling Time
```
λ_cool = c * γ / |dγ/dt_total|
t_cool = γ / |dγ/dt_total|
```

These are the distance and time for an electron to lose a factor of
1/e of its energy to all processes combined.

#### Synchrotron Power Spectrum

The synchrotron power emitted per unit frequency by the electron
population, averaged over the γ distribution:
```
P(ν) = ∫ N(γ) * p(ν, γ) dγ

p(ν, γ) = (√3 * e³ * B) / (m_e * c²) * F(ν / ν_c)

ν_c = (3 * e * B * γ²) / (4π * m_e * c)
```

where F(x) = x ∫_x^∞ K_{5/3}(t) dt is the standard synchrotron
function, evaluated numerically via scipy.special.

#### Energy Equipartition
```
U_B         = B² / (2 * μ₀)
U_rad       = a_rad * T_CMB⁴ * (1 + z)⁴
U_particles = ∫ N(γ) * γ * m_e * c² dγ
```

Comparing these three quantities is a standard diagnostic for whether
a radio galaxy lobe is in equipartition.

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

This runs the simulation and saves a single `results.json` file in the
output directory. Once this file exists, any plot can be regenerated
without re-running the simulation.

### Generate plots
```bash
# γ distribution of the electron population at t_snap
python -m plasma_sim plot --spectrum --input results/results.json

# dE/dt for each process averaged over all electrons, as a function of time
python -m plasma_sim plot --losses-time --input results/results.json

# dE/dγ for each process binned from simulation output
python -m plasma_sim plot --losses-gamma --input results/results.json

# Cooling length and cooling time as a function of γ
python -m plasma_sim plot --cooling --input results/results.json

# Synchrotron power spectrum P(ν)
python -m plasma_sim plot --spectrum-sync --input results/results.json

# Energy equipartition: U_B vs U_rad vs U_particles
python -m plasma_sim plot --equipartition --input results/results.json

# Generate all plots at once
python -m plasma_sim plot --all --input results/results.json
```

---

## Configuration

All parameters are set in a JSON file with two sections: `physical`
for the astrophysical system and `numerical` for the simulation itself.

Example `configs/default.json`:
```json
{
  "physical": {
    "B_field":        1e-9,
    "n_plasma":       1e3,
    "gamma_initial":  1e4,
    "gamma_min":      10.0,
    "epsilon_stop":   0.01,
    "t_snap":         1e12,
    "redshift":       0.0
  },
  "numerical": {
    "dt":             1e9,
    "n_electrons":    1000,
    "n_bins":         50,
    "random_seed":    42
  }
}
```

### Physical parameters

| Parameter | Description | Default | Units | Valid range |
|---|---|---|---|---|
| B_field | Magnetic field strength | 1e-9 | T | > 0 |
| n_plasma | Plasma electron number density | 1e3 | m⁻³ | > 0 |
| gamma_initial | Initial Lorentz factor | 1e4 | — | > gamma_min |
| gamma_min | Minimum γ: electron stops below this | 10.0 | — | ≥ 1 |
| epsilon_stop | Synchrotron power cutoff fraction | 0.01 | — | 0 < ε < 1 |
| t_snap | Time at which to snapshot γ distribution | 1e12 | s | > 0 |
| redshift | Source redshift (scales U_rad) | 0.0 | — | ≥ 0 |

### Numerical parameters

| Parameter | Description | Default | Units | Valid range |
|---|---|---|---|---|
| dt | Timestep | 1e9 | s | > 0, see note |
| n_electrons | Number of simulated electrons | 1000 | — | ≥ 1 |
| n_bins | Number of bins for all histograms | 50 | — | ≥ 10 |
| random_seed | Random seed for reproducibility | 42 | — | any integer |

> **Note on dt**: the timestep must be small relative to the cooling
> time t_cool = γ / |dγ/dt| to ensure accurate results. The simulation
> will raise a warning if dt > 0.01 * t_cool at the start of the run.

### Validation and warnings

The simulation checks all parameters before running and will either
raise an error (simulation cannot proceed) or a warning (simulation
proceeds but results may be unreliable):

| Condition | Type | Reason |
|---|---|---|
| B_field ≤ 0 or n_plasma ≤ 0 | Error | Unphysical |
| gamma_initial ≤ gamma_min | Error | Simulation stops immediately |
| gamma_initial > 1e7 | Warning | Thomson approximation becoming inaccurate |
| dt > 0.01 * t_cool | Warning | Timestep too large, results may be inaccurate |
| t_snap > t_cool | Warning | Most electrons will have cooled before snapshot |
| redshift < 0 | Error | Unphysical |
| epsilon_stop ≤ 0 or ≥ 1 | Error | Unphysical stopping condition |

---

## Output

The `run` command produces a single `results.json` in the output
directory with the following structure:
```json
{
  "config": { ... },

  "energy_spectrum": {
    "gamma_bins": [...],
    "counts":     [...]
  },

  "losses_vs_time": {
    "time":         [...],
    "dEdt_sync":    [...],
    "dEdt_IC":      [...],
    "dEdt_brems":   [...],
    "dEdt_coulomb": [...]
  },

  "losses_vs_gamma": {
    "gamma_bins":        [...],
    "dEdgamma_sync":     [...],
    "dEdgamma_IC":       [...],
    "dEdgamma_brems":    [...],
    "dEdgamma_coulomb":  [...]
  },

  "cooling": {
    "gamma_bins":  [...],
    "lambda_cool": [...],
    "t_cool":      [...]
  },

  "synchrotron_spectrum": {
    "nu_bins": [...],
    "power":   [...]
  },

  "equipartition": {
    "U_B":         ...,
    "U_rad":       ...,
    "U_particles": ...
  },

  "summary_statistics": {
    "gamma_mean":        ...,
    "gamma_std":         ...,
    "frac_sync_mean":    ...,
    "frac_sync_std":     ...,
    "frac_IC_mean":      ...,
    "frac_IC_std":       ...,
    "frac_brems_mean":   ...,
    "frac_brems_std":    ...,
    "frac_coulomb_mean": ...,
    "frac_coulomb_std":  ...,
    "lambda_cool_mean":  ...,
    "lambda_cool_std":   ...,
    "t_cool_mean":       ...,
    "t_cool_std":        ...
  }
}
```

---

## Testing
```bash
pytest tests/ --cov=src/ --cov-report=term-missing
```

---

## References

- Rybicki & Lightman, *Radiative Processes in Astrophysics*, Wiley (1979)
- Blumenthal & Gould, Rev. Mod. Phys. **42**, 237 (1970)
- Sarazin, Rev. Mod. Phys. **58**, 1 (1986)