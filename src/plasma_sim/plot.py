"""
Plotting functions for the plasma simulation results.
"""

import json
import numpy as np
import matplotlib.pyplot as plt


def plot_gamma_evolution(results: dict, output_path: str = None):
    """Plot the mean Lorentz factor evolution over time.

    Parameters
    ----------
    results : dict
        Results dictionary from run_simulation
    output_path : str, optional
        Path to save the figure. If None, displays interactively.
    """
    time = np.array(results["time_series"]["time"])
    gamma_mean = np.array(results["time_series"]["gamma_mean"])
    n_alive = np.array(results["time_series"]["n_alive"])

    time_years = time / 3.15e7

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    ax1.loglog(time_years, gamma_mean, color="steelblue", linewidth=2)
    ax1.set_ylabel("Mean Lorentz factor γ")
    ax1.set_title("Electron population evolution")
    ax1.grid(True, alpha=0.3)

    ax2.semilogx(time_years, n_alive, color="coral", linewidth=2)
    ax2.set_xlabel("Time (years)")
    ax2.set_ylabel("Number of living electrons")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150)
        print(f"Saved to {output_path}")
    else:
        plt.show()


def plot_loss_rates_vs_time(results: dict, output_path: str = None):
    """Plot the energy loss rates for each process over time.

    Parameters
    ----------
    results : dict
        Results dictionary from run_simulation
    output_path : str, optional
        Path to save the figure. If None, displays interactively.
    """
    time = np.array(results["time_series"]["time"])
    time_years = time / 3.15e7

    sync = np.abs(results["time_series"]["dEdt_sync_mean"])
    IC = np.abs(results["time_series"]["dEdt_IC_mean"])
    brems = np.abs(results["time_series"]["dEdt_brems_mean"])
    coulomb = np.abs(results["time_series"]["dEdt_coulomb_mean"])

    plt.figure(figsize=(10, 6))
    plt.loglog(time_years, sync, label="Synchrotron", linewidth=2)
    plt.loglog(time_years, IC, label="Inverse Compton", linewidth=2)
    plt.loglog(time_years, brems, label="Bremsstrahlung", linewidth=2)
    plt.loglog(time_years, coulomb, label="Coulomb", linewidth=2)

    plt.xlabel("Time (years)")
    plt.ylabel("|dγ/dt| (s⁻¹)")
    plt.title("Energy loss rates vs time")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150)
        print(f"Saved to {output_path}")
    else:
        plt.show()


def plot_loss_rates_vs_gamma(results: dict, output_path: str = None):
    """Plot the energy loss rates for each process as a function of gamma.

    Parameters
    ----------
    results : dict
        Results dictionary from run_simulation
    output_path : str, optional
        Path to save the figure. If none, displays interactively.
    """
    gamma_bins = np.array(results["losses_vs_gamma"]["gamma_bins"])
    sync = np.abs(results["losses_vs_gamma"]["dEdgamma_sync"])
    IC = np.abs(results["losses_vs_gamma"]["dEdgamma_IC"])
    brems = np.abs(results["losses_vs_gamma"]["dEdgamma_brems"])
    coulomb = np.abs(results["losses_vs_gamma"]["dEdgamma_coulomb"])

    plt.figure(figsize=(10, 6))
    plt.loglog(gamma_bins, sync, label="Synchrotron", linewidth=2)
    plt.loglog(gamma_bins, IC, label="Inverse Compton", linewidth=2)
    plt.loglog(gamma_bins, brems, label="Bremsstrahlung", linewidth=2)
    plt.loglog(gamma_bins, coulomb, label="Coulomb", linewidth=2)

    plt.xlabel("Lorentz factor γ")
    plt.ylabel("|dγ/dt| (s⁻¹)")
    plt.title("Energy loss rates vs Lorentz factor")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150)
        print(f"Saved to {output_path}")
    else:
        plt.show()


def plot_energy_fractions_vs_time(results: dict, output_path: str = None):
    """Plot the fractional contribution of each process to total energy loss vs time.

    Parameters
    ----------
    results : dict
        Results dictionary from run_simulation
    output_path : str, optional
        Path to save the figure. If none, displays interactively.
    """
    time = np.array(results["time_series"]["time"])
    time_years = time / 3.15e7

    sync = np.abs(results["time_series"]["dEdt_sync_mean"])
    IC = np.abs(results["time_series"]["dEdt_IC_mean"])
    brems = np.abs(results["time_series"]["dEdt_brems_mean"])
    coulomb = np.abs(results["time_series"]["dEdt_coulomb_mean"])

    total = np.array(sync) + np.array(IC) + np.array(brems) + np.array(coulomb)
    total = np.where(total == 0, 1, total)  # avoid division by zero

    frac_sync = sync / total
    frac_IC = IC / total
    frac_brems = brems / total
    frac_coulomb = coulomb / total

    plt.figure(figsize=(10, 6))
    plt.stackplot(time_years,
                  frac_sync, frac_IC, frac_brems, frac_coulomb,
                  labels=["Synchrotron", "Inverse Compton",
                          "Bremsstrahlung", "Coulomb"],
                  alpha=0.8)
    plt.xscale("log")
    plt.xlabel("Time (years)")
    plt.ylabel("Fractional energy loss")
    plt.title("Energy loss fractions vs time")
    plt.legend(loc="upper left")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150)
        print(f"Saved to {output_path}")
    else:
        plt.show()


def plot_energy_budget(results: dict, output_path: str = None):
    """Plot a pie chart of total energy lost to each process.

    Parameters
    ----------
    results : dict
        Results dictionary from run_simulation
    output_path : str, optional
        Path to save the figure. If None, displays interactively.
    """
    fracs = [
        results["summary"]["frac_sync_mean"],
        results["summary"]["frac_IC_mean"],
        results["summary"]["frac_brems_mean"],
        results["summary"]["frac_coulomb_mean"],
    ]
    labels = ["Synchrotron", "Inverse Compton", "Bremsstrahlung", "Coulomb"]
    colors = ["steelblue", "orange", "green", "red"]

    plt.figure(figsize=(8, 8))
    plt.pie(fracs, labels=labels, colors=colors, autopct="%1.2f%%",
            startangle=90, pctdistance=0.85)
    plt.title("Total energy budget")
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150)
        print(f"Saved to {output_path}")
    else:
        plt.show()


def plot_IC_sync_ratio(results: dict, output_path: str = None):
    """Plot the ratio of IC to synchrotron loss rate vs time.

    Both processes scale as gamma^2 so the ratio should be roughly constant
    and equal to U_rad / U_B. Deviations reveal the effect of B fluctuations.

    Parameters
    ----------
    results : dict
        Results dictionary from run_simulation
    output_path : str, optional
        Path to save the figure. If None, displays interactively.
    """
    time = np.array(results["time_series"]["time"])
    time_years = time / 3.15e7

    sync = np.abs(np.array(results["time_series"]["dEdt_sync_mean"]))
    IC = np.abs(np.array(results["time_series"]["dEdt_IC_mean"]))

    ratio = IC / np.where(sync == 0, 1, sync)

    plt.figure(figsize=(10, 6))
    plt.semilogx(time_years, ratio, color="purple", linewidth=2)
    plt.xlabel("Time (years)")
    plt.ylabel("IC / Synchrotron loss rate")
    plt.title("Inverse Compton to Synchrotron ratio vs time")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150)
        print(f"Saved to {output_path}")
    else:
        plt.show()


def plot_all(results: dict, output_dir: str = "results/"):
    """Generate all plots and save them to the output directory.

    Parameters
    ----------
    results : dict
        Results dictionary from run_simulation
    output_dir : str
        Directory to save the figures
    """
    plot_gamma_evolution(results, output_path=f"{output_dir}gamma_evolution.png")
    plot_loss_rates_vs_time(results, output_path=f"{output_dir}loss_rates_time.png")
    plot_loss_rates_vs_gamma(results, output_path=f"{output_dir}loss_rates_gamma.png")
    plot_energy_fractions_vs_time(results, output_path=f"{output_dir}energy_fractions.png")
    plot_energy_budget(results, output_path=f"{output_dir}energy_budget.png")
    plot_IC_sync_ratio(results, output_path=f"{output_dir}IC_sync_ratio.png")
    print("All plots saved!")