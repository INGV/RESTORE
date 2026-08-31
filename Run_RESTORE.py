"""
Run RESTORE for one or more stochastic realizations.
Input parameters are defined in "input_file.txt".

EXAMPLE
-------
# Run with a specific input file
restore --config example/amatrice/input_file.txt

# Show the installed RESTORE version
restore --version
"""

import argparse
import hashlib
from itertools import product
import json
import os
import platform
import shutil
import traceback
import numpy as np
import pandas as pd
import config
import RESTORE

CATALOG_COLUMNS = [
    "longitude",
    "latitude",
    "year",
    "month",
    "day",
    "magnitude",
    "depth_km",
    "hour",
    "minute",
    "second",
    "flag",
]

REPRESENTATIVE_METRIC_WEIGHTS = {
    # Number of eqs reconstructed
    "n_events_reconstructed": 1.0,
    # Magnitude distribution: total weight = 1
    "_reconstructed_magnitude_median": 1.0,
    # Spatial distribution: total weight = 1
    "_reconstructed_latitude_median": 0.5,
    "_reconstructed_longitude_median": 0.5,
    # Depth distribution: total weight = 1
    "_reconstructed_depth_km_median": 1.0,
    # STAI characteristics: total weight = 1
    "total_stai_duration_days": 0.5,
    "n_stai_gaps": 0.5,
}


def make_folder(folder):
    os.makedirs(folder, exist_ok=True)


def _file_sha256(path):
    """
    Provide a unique identifier for the supplied file, useful for reproducibility.

    RESTORE uses it for the config and the catalog files, and stores it in the realization metadata
    --> the outputs of each realization can always be associated to the config and catalog files used to generate them.
    """

    digest = hashlib.sha256()
    with open(path, "rb") as input_file:
        for chunk in iter(lambda: input_file.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_run_metadata(params, run_folder):
    """
    Write provenance for code, environment, and input files.
    """

    metadata = {
        "restore_version": RESTORE.__version__,
        "python_version": platform.python_version(),
        "config_file": params["config_file"],
        "config_sha256": _file_sha256(params["config_file"]),
        "catalog_file": params["name_catalog"],
        "catalog_sha256": _file_sha256(params["name_catalog"]),
    }
    with open(os.path.join(run_folder, "metadata.json"), "w") as output:
        json.dump(metadata, output, indent=2)


def read_replenished_catalog(catalog_file):
    catalog = pd.read_csv(
        catalog_file,
        sep=r"\s+",
        header=None,
        names=CATALOG_COLUMNS,
        engine="python",
    )
    for column in CATALOG_COLUMNS:
        catalog[column] = pd.to_numeric(catalog[column])

    catalog["datetime"] = pd.to_datetime(
        dict(
            year=catalog["year"].astype(int),
            month=catalog["month"].astype(int),
            day=catalog["day"].astype(int),
            hour=catalog["hour"].astype(int),
            minute=catalog["minute"].astype(int),
        )
    ) + pd.to_timedelta(catalog["second"], unit="s")
    return catalog


def _format_min_median_max(values):
    if len(values) == 0:
        return ""

    min_value = np.nanmin(values)
    median_value = np.nanmedian(values)
    max_value = np.nanmax(values)
    return f"min={min_value:.4g}; median={median_value:.4g}; max={max_value:.4g}"


def add_reconstructed_stats(row, data, columns):
    """
    Add reconstructed eqs statistics.
    """

    for column in columns:
        values = data[column].to_numpy(dtype=float)
        row["reconstructed_" + column] = _format_min_median_max(values)
        row["_reconstructed_" + column + "_median"] = np.nanmedian(values) if len(values) else np.nan


def read_run_summary(run_folder):
    """
    Read the JSON summary written by RESTORE for one realization.
    """

    summary_file = os.path.join(run_folder, "RESTORE_run_summary.json")
    if not os.path.exists(summary_file):
        return {}
    with open(summary_file, "r", encoding="utf-8") as input_file:
        return json.load(input_file)


def summarize_catalog(catalog, seed, run_summary=None):
    """
    Summarize one replenished catalog realization.

    The function returns a dictionary for a given realization 
    (a row of the ``realization_metrics.csv`` file that will be created later).
    It aggregates information taken from ``run_summary`` with reconstructed eqs statistics
    for magnitude, longitude, latitude, and hypocentral depth.
    """

    run_summary = run_summary or {}
    reconstructed = catalog[catalog["flag"] == 1].copy()
    row = {
        "seed": seed,
        "status": "ok",
        "error_type": "",
        "error_message": "",
        "n_events_total": len(catalog),
        "n_events_reconstructed": len(reconstructed),
        "fraction_reconstructed": len(reconstructed) / len(catalog) if len(catalog) else np.nan,
        "n_stai_gaps": run_summary.get("n_stai_gaps", np.nan),
        "total_stai_duration_days": run_summary.get("total_stai_duration_days", np.nan),
        "reference_mc": run_summary.get("reference_mc", np.nan),
        "b_value": run_summary.get("b_value", np.nan),
        "b_value_source": run_summary.get("b_value_source", ""),
    }

    metric_columns = ["magnitude", "longitude", "latitude", "depth_km"]
    add_reconstructed_stats(row, reconstructed, metric_columns)

    return row


def metric_distance(frame, center, metrics):
    """
    Calculate how far each row (= realization) is from the ensemble central behavior.

    For every metric, the difference from the central value is normalized by that metric's
    interquartile range (IQR). This prevents variables with larger numerical
    units from dominating the comparison.
    The final distance is the mean absolute normalized difference across
    the representative metrics. Smaller values indicate behavior closer to the
    ensemble center.

    NOTE: If a metric has zero or undefined IQR, a scale of 1 is used to avoid
    division by zero.
    """

    data = frame[metrics].apply(pd.to_numeric, errors="coerce")
    center = center[metrics].apply(pd.to_numeric, errors="coerce")
    scales = data.quantile(0.75) - data.quantile(0.25)
    scales = scales.replace(0, 1.0).fillna(1.0)

    normalized_distances = ((data - center) / scales).abs()
    weights = pd.Series({metric: REPRESENTATIVE_METRIC_WEIGHTS[metric] for metric in metrics})

    return normalized_distances.mul(weights).sum(axis=1) / weights.sum()


def choose_representative(realization_metrics):
    """
    Select the realization closest to the central behavior of the ensemble.

    Rationale:
    every successful (size, step, seed) combination = one RESTORE realization
    --> calculate ensemble center (= median of each representative metric across all successful realizations)
    --> calculate distance of every realization from that center across the representative metrics after 
    normalizing each metric by its interquartile range (IQR)
    --> the realization with the smallest mean normalized distance is selected as representative
    """

    successful = realization_metrics[realization_metrics["status"] == "ok"].copy()
    if len(successful) == 0:
        return None, realization_metrics

    metrics = list(REPRESENTATIVE_METRIC_WEIGHTS)
    missing_metrics = [metric for metric in metrics if metric not in successful.columns]

    if missing_metrics:
        raise ValueError(
            "Missing representative metrics: " + ", ".join(missing_metrics))

    # Central behavior of the ensemble
    central_behavior = successful[metrics].median()

    # Distance of every successful realization from the ensemble center
    distances = metric_distance(successful, central_behavior, metrics)

    # Store the distance in "realization_metrics.csv"
    successful["distance_from_ensemble_center"] = distances

    realization_metrics = realization_metrics.copy()
    realization_metrics["distance_from_ensemble_center"] = np.nan

    realization_metrics.loc[successful.index, "distance_from_ensemble_center"] = successful["distance_from_ensemble_center"]

    # Select the most central realization
    representative_index = distances.idxmin()
    representative = successful.loc[representative_index].to_dict()
    representative["selection_metrics"] = ",".join(metrics)

    return representative, realization_metrics


def format_settings(params, seed=None):
    """
    Return RESTORE settings as a text report.
    """

    lines = [
        "----- RESTORE SETTINGS -----",
        f"Configuration file = {params['config_file']}",
        f"Catalog = {params['name_catalog']}",
        f"Catalog delimiter = {params['delimiter']!r}",
        f"Window size = {params['size']} events",
        f"Window step = {params['step']} events",
        f"Seed = {seed}",
        f"End of quiet period = {params['t_end_quiet']}",
        f"Reference Mc = {params['mc'] if params['mc'] is not None else 'estimated'}",
        f"b-value = {params['b-value'] if params['b-value'] is not None else '1.0 (default)'}",
        f"Confidence level for STAI gaps = {params['st_dev_multiplier']} * Mc Sigma",
        f"Sigma (smoothing kernel) = {params['Sigma']}",
        f"sbin (spatial grid resolution) = {params['sbin']}°",
        f"Fault-length multiplier = {params['fault_length_multiplier']}",
        f"alpha Lilliefors test = {params['alpha']}",
        f"Hypocentral-depth distribution = {params['depth_distribution']}",
        f"Depth-distribution initial params = {params['p0']}",
        f"Output folder = {params['output_folder']}",
    ]

    return "\n".join(lines) + "\n"


def format_run_overview(params, n_parameter_sets):
    """
    Return an overview of the requested run(s).
    """

    lines = [
        "---------------------",
        "--- RESTORE INPUT ---",
        "---------------------",
        f"Configuration file: {params['config_file']}",
        f"Catalog: {params['name_catalog']}",
        f"Output folder: {params['output_folder']}",
        f"Window sizes: {params['sizes']}",
        f"Window steps: {params['steps']}",
        f"Seeds: {params['seeds']}",
        f"Parameter combinations: {n_parameter_sets}",
        f"Total realizations: {n_parameter_sets * len(params['seeds'])}",
    ]
    return "\n".join(lines) + "\n"


def write_settings_report(params, seed, run_folder):
    """
    Write the settings for one realization.
    """

    report_file = os.path.join(run_folder, "settings.txt")
    with open(report_file, "w", encoding="utf-8") as output:
        output.write(format_settings(params, seed=seed))
    return report_file


def run_restore_once(params, seed, run_folder):
    """
    Run one RESTORE realization for a given size, step, and random seed.
    """

    make_folder(run_folder)
    rng = np.random.default_rng(seed)

    shutil.copy2(params["config_file"], os.path.join(run_folder, "input_file.txt"))
    write_settings_report(params, seed, run_folder)
    write_run_metadata(params, run_folder)

    catalog_sel = RESTORE.read_catalog_zmap(params["name_catalog"], params["delimiter"])
    serial_times = RESTORE.serial_time(catalog_sel)
    t_start = RESTORE.serial_time(params["t_end_quiet"])

    old_cwd = os.getcwd()
    os.chdir(run_folder)
    try:
        RESTORE.replenished_catalog(
            catalog_sel,
            t_start,
            params["mc"],
            params["b-value"],
            params["size"],
            params["step"],
            params["st_dev_multiplier"],
            params["alpha"],
            serial_times,
            params["Sigma"],
            params["sbin"],
            params["fault_length_multiplier"],
            params["depth_distribution"],
            params["p0"],
            rng=rng,
        )
    finally:
        os.chdir(old_cwd)

    return os.path.join(run_folder, "Replenished_catalog.txt")


def parameter_label(params):
    return "size_" + str(params["size"]) + "_step_" + str(params["step"])


def run_seed(params, seed, realizations_folder):
    """
    Run and summarize one RESTORE realization for a given size, step, and seed.
    """
        
    param_label = parameter_label(params)
    run_folder = os.path.join(realizations_folder, param_label, "seed_" + str(seed))
    make_folder(run_folder)
    print("-----------------------------")
    print(f"Running realization:\n"
        f"size={params['size']}, step={params['step']}, seed={seed}")
    print("-----------------------------")
    try:
        catalog_file = run_restore_once(params, seed, run_folder)
        catalog = read_replenished_catalog(catalog_file)
        if (catalog["flag"] == 1).sum() == 0:
            raise ValueError(
                "STAI gaps were detected, but no events were reconstructed."
            )
        row = summarize_catalog(catalog, seed, read_run_summary(run_folder))
        row["catalog_file"] = os.path.abspath(catalog_file)
    except RESTORE.NoSTAIGapsFound as err:
        print(f"Seed {seed}: no STAI gaps found.")
        print(err)
        row = {
            "seed": seed,
            "status": "no_gaps_found",
            "error_type": type(err).__name__,
            "error_message": str(err),
            "catalog_file": "",
        }
    except Exception as err:
        print(f"Seed {seed}: RESTORE failed.")
        print(f"{type(err).__name__}: {err}")
        row = {
            "seed": seed,
            "status": "failed",
            "error_type": type(err).__name__,
            "error_message": str(err),
            "catalog_file": "",
            "traceback": traceback.format_exc(),
        }
    row["size"] = params["size"]
    row["step"] = params["step"]
    return row


def build_parameter_sets(params):
    parameter_sets = []
    for size, step in product(params["sizes"], params["steps"]):
        item = dict(params)
        item["size"] = int(size)
        item["step"] = int(step)
        parameter_sets.append(item)
    return parameter_sets


def print_run_overview(params, n_parameter_sets):
    print(format_run_overview(params, n_parameter_sets), end="")


def main():

    parser = argparse.ArgumentParser(
        description="Run RESTORE stochastic realizations."
    )
    parser.add_argument(
        "--config",
        default=config.DEFAULT_CONFIG_FILE,
        help="Path to the RESTORE input file.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="RESTORE " + RESTORE.__version__,
    )
    args = parser.parse_args()

    # Load user-defined input params
    params = config.load_config(args.config)

    make_folder(params["output_folder"])
    realizations_folder = os.path.join(params["output_folder"], "realizations")
    make_folder(realizations_folder)

    # Build all requested size/step combinations 
    # --> each run uses one size, one step and one seed
    parameter_sets = build_parameter_sets(params)

    print_run_overview(params, len(parameter_sets))

    rows = [
        run_seed(run_params, seed, realizations_folder)
        for run_params in parameter_sets
        for seed in params["seeds"]]

    realization_metrics = pd.DataFrame(rows)
    n_successful = (realization_metrics["status"] == "ok").sum()
    n_total = len(realization_metrics)

    representative, realization_metrics = choose_representative(realization_metrics)
    realization_metrics_file = os.path.join(params["output_folder"], "realization_metrics.csv")

    # Internal numeric columns prefixed with "_" are used for representative
    # selection only and are removed before writing realization_metrics.csv.
    realization_metrics = realization_metrics[
        [column for column in realization_metrics.columns if not column.startswith("_")]]

    first_columns = ["size", "step", "seed", "status"]
    other_columns = [
        column
        for column in realization_metrics.columns
        if column not in first_columns]
    realization_metrics = realization_metrics[first_columns + other_columns]
    realization_metrics.to_csv(realization_metrics_file, index=False)

    representative_file = os.path.join(params["output_folder"], "representative_realization.txt")

    with open(representative_file, "w") as output:
        if representative is None:
            output.write("No successful RESTORE realization available.\n")
        else:
            output.write("Representative RESTORE realization\n")
            output.write("size: " + str(int(representative["size"])) + "\n")
            output.write("step: " + str(int(representative["step"])) + "\n")
            output.write("seed: " + str(int(representative["seed"])) + "\n")
            output.write("distance_from_ensemble_center: " + str(representative["distance_from_ensemble_center"]) + "\n")
            output.write("catalog_file: " + str(representative["catalog_file"]) + "\n")
            output.write("selection_metrics: " + str(representative["selection_metrics"]) + "\n")

    print("------------------------")
    print("RESTORE run(s): complete")
    print("------------------------")
    print(f"Successful realizations: {n_successful}/{n_total}")

    if n_successful == 0:
        print(
            "Warning: no RESTORE realization completed successfully."
        )
    if representative is not None:
        print(
            "Representative realization: "
            f"size={int(representative['size'])}, "
            f"step={int(representative['step'])}, "
            f"seed={int(representative['seed'])}"
        )

    print(f"Realization metrics: {realization_metrics_file}")
    print(f"Representative file: {representative_file}")

if __name__ == "__main__":
    main()
