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

    Each electron evolves independently with an adaptive timestep computed
    from its current cooling time. Results are interpolated onto shared
    output times for population averaging.

    Parameters
    ----------
    config : dict
        Configuration dictionary with 'physical' and 'numerical' sections
        as described in the README.

    Returns
    -------
    dict
        Contains time series of population-averaged quantities and
        energy loss distributions binned by gamma.
    """
    phys = config["physical"]
    num = config["numerical"]

    B_field = phys["B_field"]
    sigma_B = phys["sigma_B"]
    n_plasma = phys["n_plasma"]
    gamma_min = phys["gamma_min"]
    epsilon_stop = phys["epsilon_stop"]
    redshift = phys["redshift"]
    n_electrons = num["n_electrons"]
    B_update_steps = num["B_update_steps"]
    n_bins = num["n_bins"]

    rng = np.random.default_rng(num["random_seed"])

    # sample initial gamma values from power law
    gammas_initial = sample_gamma_powerlaw(
        gamma_min=phys["gamma_min_init"],
        gamma_max=phys["gamma_max_init"],
        spectral_index=phys["spectral_index"],
        n_electrons=n_electrons,
        rng=rng
    )

    # storage for all electron histories
    all_times = []
    all_gammas = []
    all_dEdt_sync = []
    all_dEdt_IC = []
    all_dEdt_brems = []
    all_dEdt_coulomb = []
    t_cool_list = []
    gamma_final_list = []
    energy_sync_list = []
    energy_IC_list = []
    energy_brems_list = []
    energy_coulomb_list = []

    for i in range(n_electrons):
        gamma = gammas_initial[i]
        B = B_field
        t = 0.0
        step = 0

        P_initial = abs(synchrotron_loss(gamma, B))

        # per-electron history
        times = [t]
        gammas = [gamma]
        dEdt_sync_hist = []
        dEdt_IC_hist = []
        dEdt_brems_hist = []
        dEdt_coulomb_hist = []

        # cumulative energy lost to each process
        E_sync = 0.0
        E_IC = 0.0
        E_brems = 0.0
        E_coulomb = 0.0

        while True:
            # update B every B_update_steps
            if step % B_update_steps == 0:
                B = sample_new_B(B_field, sigma_B, rng)

            # compute losses
            dg_sync = synchrotron_loss(gamma, B)
            dg_IC = inverse_compton_loss(gamma, redshift)
            dg_brems = bremsstrahlung_loss(gamma, n_plasma)
            dg_coulomb = 0.0

            # stochastic Coulomb — probability computed below with adaptive dt
            dg_total = dg_sync + dg_IC + dg_brems

            # adaptive timestep — 1% of cooling time
            dt = 0.01 * gamma / abs(dg_total)

            # now check Coulomb with correct dt
            p = compute_interaction_probability(gamma, n_plasma, dt)
            if rng.uniform() < p:
                dg_coulomb = coulomb_loss(gamma, n_plasma)
                dg_total += dg_coulomb

            # update gamma
            gamma += dg_total * dt
            t += dt
            step += 1

            # accumulate energy losses
            E_sync += abs(dg_sync) * dt
            E_IC += abs(dg_IC) * dt
            E_brems += abs(dg_brems) * dt
            E_coulomb += abs(dg_coulomb) * dt

            # record history
            times.append(t)
            gammas.append(gamma)
            dEdt_sync_hist.append(dg_sync)
            dEdt_IC_hist.append(dg_IC)
            dEdt_brems_hist.append(dg_brems)
            dEdt_coulomb_hist.append(dg_coulomb)

            # check stopping conditions
            P_sync = abs(synchrotron_loss(gamma, B))
            if gamma < gamma_min or P_sync < epsilon_stop * P_initial:
                break

        # store electron history
        all_times.append(times)
        all_gammas.append(gammas)
        all_dEdt_sync.append(dEdt_sync_hist)
        all_dEdt_IC.append(dEdt_IC_hist)
        all_dEdt_brems.append(dEdt_brems_hist)
        all_dEdt_coulomb.append(dEdt_coulomb_hist)
        t_cool_list.append(t)
        gamma_final_list.append(gamma)
        energy_sync_list.append(E_sync)
        energy_IC_list.append(E_IC)
        energy_brems_list.append(E_brems)
        energy_coulomb_list.append(E_coulomb)

    # build shared time axis (log spaced from min to max cooling time)
    t_min = min(t[1] for t in all_times)
    t_max = max(t[-1] for t in all_times)
    shared_times = np.logspace(np.log10(t_min), np.log10(t_max), 200)

    # interpolate each electron onto shared times
    gamma_interp = np.zeros((n_electrons, len(shared_times)))
    dEdt_sync_interp = np.zeros((n_electrons, len(shared_times)))
    dEdt_IC_interp = np.zeros((n_electrons, len(shared_times)))
    dEdt_brems_interp = np.zeros((n_electrons, len(shared_times)))
    dEdt_coulomb_interp = np.zeros((n_electrons, len(shared_times)))
    alive_mask = np.zeros((n_electrons, len(shared_times)), dtype=bool)

    for i in range(n_electrons):
        t_arr = np.array(all_times[i])
        g_arr = np.array(all_gammas[i])
        sync_arr = np.array([all_dEdt_sync[i][0]] + all_dEdt_sync[i])
        IC_arr = np.array([all_dEdt_IC[i][0]] + all_dEdt_IC[i])
        brems_arr = np.array([all_dEdt_brems[i][0]] + all_dEdt_brems[i])
        coulomb_arr = np.array([all_dEdt_coulomb[i][0]] + all_dEdt_coulomb[i])

        # only interpolate within this electron's lifetime
        mask = (shared_times >= t_arr[0]) & (shared_times <= t_arr[-1])
        alive_mask[i, mask] = True

        if np.sum(mask) > 0:
            gamma_interp[i, mask] = np.interp(
                shared_times[mask], t_arr, g_arr
            )
            dEdt_sync_interp[i, mask] = np.interp(
                shared_times[mask], t_arr, sync_arr
            )
            dEdt_IC_interp[i, mask] = np.interp(
                shared_times[mask], t_arr, IC_arr
            )
            dEdt_brems_interp[i, mask] = np.interp(
                shared_times[mask], t_arr, brems_arr
            )
            dEdt_coulomb_interp[i, mask] = np.interp(
                shared_times[mask], t_arr, coulomb_arr
            )

    # compute population averages at each shared time
    n_alive = np.sum(alive_mask, axis=0)
    valid = n_alive > 0

    gamma_mean = np.zeros(len(shared_times))
    dEdt_sync_mean = np.zeros(len(shared_times))
    dEdt_IC_mean = np.zeros(len(shared_times))
    dEdt_brems_mean = np.zeros(len(shared_times))
    dEdt_coulomb_mean = np.zeros(len(shared_times))

    for t_idx in range(len(shared_times)):
        if n_alive[t_idx] > 0:
            mask = alive_mask[:, t_idx]
            gamma_mean[t_idx] = np.mean(gamma_interp[mask, t_idx])
            dEdt_sync_mean[t_idx] = np.mean(dEdt_sync_interp[mask, t_idx])
            dEdt_IC_mean[t_idx] = np.mean(dEdt_IC_interp[mask, t_idx])
            dEdt_brems_mean[t_idx] = np.mean(dEdt_brems_interp[mask, t_idx])
            dEdt_coulomb_mean[t_idx] = np.mean(
                dEdt_coulomb_interp[mask, t_idx]
            )

    # build losses vs gamma by binning all history data
    all_g = np.concatenate([all_gammas[i][:-1] for i in range(n_electrons)])
    all_s = np.concatenate([all_dEdt_sync[i] for i in range(n_electrons)])
    all_ic = np.concatenate([all_dEdt_IC[i] for i in range(n_electrons)])
    all_b = np.concatenate([all_dEdt_brems[i] for i in range(n_electrons)])
    all_c = np.concatenate([all_dEdt_coulomb[i] for i in range(n_electrons)])

    gamma_bins = np.logspace(
        np.log10(gamma_min),
        np.log10(phys["gamma_max_init"]),
        n_bins + 1
    )
    bin_indices = np.digitize(all_g, gamma_bins)

    dEdgamma_sync = np.zeros(n_bins)
    dEdgamma_IC = np.zeros(n_bins)
    dEdgamma_brems = np.zeros(n_bins)
    dEdgamma_coulomb = np.zeros(n_bins)

    for b in range(1, n_bins + 1):
        mask = bin_indices == b
        if np.sum(mask) > 0:
            dEdgamma_sync[b-1] = np.mean(all_s[mask])
            dEdgamma_IC[b-1] = np.mean(all_ic[mask])
            dEdgamma_brems[b-1] = np.mean(all_b[mask])
            dEdgamma_coulomb[b-1] = np.mean(all_c[mask])

    gamma_bin_centers = 0.5 * (gamma_bins[:-1] + gamma_bins[1:])

    # summary statistics
    t_cool_arr = np.array(t_cool_list)
    E_total = np.array([
        energy_sync_list[i] + energy_IC_list[i] +
        energy_brems_list[i] + energy_coulomb_list[i]
        for i in range(n_electrons)
    ])

    return {
        "time_series": {
            "time": shared_times[valid].tolist(),
            "n_alive": n_alive[valid].tolist(),
            "gamma_mean": gamma_mean[valid].tolist(),
            "dEdt_sync_mean": dEdt_sync_mean[valid].tolist(),
            "dEdt_IC_mean": dEdt_IC_mean[valid].tolist(),
            "dEdt_brems_mean": dEdt_brems_mean[valid].tolist(),
            "dEdt_coulomb_mean": dEdt_coulomb_mean[valid].tolist(),
        },
        "losses_vs_gamma": {
            "gamma_bins": gamma_bin_centers.tolist(),
            "dEdgamma_sync": dEdgamma_sync.tolist(),
            "dEdgamma_IC": dEdgamma_IC.tolist(),
            "dEdgamma_brems": dEdgamma_brems.tolist(),
            "dEdgamma_coulomb": dEdgamma_coulomb.tolist(),
        },
        "summary": {
            "t_cool_mean": float(np.mean(t_cool_arr)),
            "t_cool_std": float(np.std(t_cool_arr)),
            "gamma_final_mean": float(np.mean(gamma_final_list)),
            "gamma_final_std": float(np.std(gamma_final_list)),
            "frac_sync_mean": float(np.mean(
                np.array(energy_sync_list) / E_total
            )),
            "frac_IC_mean": float(np.mean(
                np.array(energy_IC_list) / E_total
            )),
            "frac_brems_mean": float(np.mean(
                np.array(energy_brems_list) / E_total
            )),
            "frac_coulomb_mean": float(np.mean(
                np.array(energy_coulomb_list) / E_total
            )),
        }
    }