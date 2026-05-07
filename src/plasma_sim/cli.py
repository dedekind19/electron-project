"""
Command line interface for the plasma simulation.

Usage:
    python -m plasma_sim run --config configs/default.json --output results/
    python -m plasma_sim plot --input results/results.json --output results/
"""

import argparse
import json
import os
import sys
import warnings

from plasma_sim.simulation import run_simulation
from plasma_sim.plot import plot_all


def validate_config(config: dict) -> None:
    """Validate the configuration dictionary and raise errors or warnings.

    Parameters
    ----------
    config : dict
        Configuration dictionary with 'physical' and 'numerical' sections.

    Raises
    ------
    ValueError
        If any parameter is unphysical or would cause the simulation to fail.
    """
    phys = config["physical"]
    num = config["numerical"]

    # errors — simulation cannot proceed
    if phys["B_field"] <= 0:
        raise ValueError("B_field must be positive.")
    if phys["n_plasma"] <= 0:
        raise ValueError("n_plasma must be positive.")
    if phys.get("sigma_B", 0) < 0:
        raise ValueError("sigma_B must be non-negative.")
    if phys["redshift"] < 0:
        raise ValueError("redshift must be non-negative.")
    if phys["epsilon_stop"] <= 0 or phys["epsilon_stop"] >= 1:
        raise ValueError("epsilon_stop must be between 0 and 1.")
    if phys["gamma_min_init"] >= phys["gamma_max_init"]:
        raise ValueError("gamma_min_init must be less than gamma_max_init.")
    if phys["spectral_index"] <= 1:
        raise ValueError("spectral_index must be greater than 1.")
    if phys["gamma_min"] < 1:
        raise ValueError("gamma_min must be >= 1.")
    if num["n_electrons"] < 1:
        raise ValueError("n_electrons must be at least 1.")
    if num["n_bins"] < 10:
        raise ValueError("n_bins must be at least 10.")
    if num["B_update_steps"] < 1:
        raise ValueError("B_update_steps must be at least 1.")

    # warnings — simulation can proceed but results may be unreliable
    if phys["gamma_max_init"] > 1e7:
        warnings.warn(
            "gamma_max_init > 1e7: Thomson approximation becoming inaccurate.",
            UserWarning
        )
    if phys["gamma_min_init"] < phys["gamma_min"]:
        warnings.warn(
            "gamma_min_init < gamma_min: some electrons may stop immediately.",
            UserWarning
        )


def cmd_run(args):
    """Run the simulation with the given config file."""

    # load config
    if not os.path.exists(args.config):
        print(f"Error: config file '{args.config}' not found.")
        sys.exit(1)

    with open(args.config) as f:
        config = json.load(f)

    # validate
    try:
        validate_config(config)
    except ValueError as e:
        print(f"Error in config: {e}")
        sys.exit(1)

    # create output directory if it does not exist
    os.makedirs(args.output, exist_ok=True)

    # run simulation
    print(f"Running simulation with {config['numerical']['n_electrons']} electrons...")
    results = run_simulation(config)
    results["config"] = config

    # save results
    output_file = os.path.join(args.output, "results.json")
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"Results saved to {output_file}")


def cmd_plot(args):
    """Generate plots from an existing results file."""

    if not os.path.exists(args.input):
        print(f"Error: results file '{args.input}' not found.")
        sys.exit(1)

    with open(args.input) as f:
        results = json.load(f)

    os.makedirs(args.output, exist_ok=True)
    plot_all(results, output_dir=args.output + "/")
    print(f"Plots saved to {args.output}/")


def main():
    """Main entry point for the command line interface."""

    parser = argparse.ArgumentParser(
        description="Monte Carlo simulation of relativistic electron energy losses."
    )
    subparsers = parser.add_subparsers(dest="command")

    # run command
    run_parser = subparsers.add_parser("run", help="Run the simulation")
    run_parser.add_argument(
        "--config", required=True,
        help="Path to the JSON configuration file"
    )
    run_parser.add_argument(
        "--output", default="results/",
        help="Directory to save results (default: results/)"
    )

    # plot command
    plot_parser = subparsers.add_parser("plot", help="Generate plots from results")
    plot_parser.add_argument(
        "--input", required=True,
        help="Path to the results JSON file"
    )
    plot_parser.add_argument(
        "--output", default="results/",
        help="Directory to save plots (default: results/)"
    )

    args = parser.parse_args()

    if args.command == "run":
        cmd_run(args)
    elif args.command == "plot":
        cmd_plot(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()