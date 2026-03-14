# %% imports
import logging
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
from scipy.stats import t
from config import rdir, replication_seed
from scow.paths import coefs_path, loglik_path


def _raise_if_nan(frame, columns, run_name, artifact_name):
    nan_mask = frame[columns].isna()
    if not nan_mask.any().any():
        return

    failing_rows = frame.loc[nan_mask.any(axis=1), :]
    preview = failing_rows.head(5).to_dict(orient="records")
    raise ValueError(
        f"NaN detected in {artifact_name} for run {run_name}. "
        f"Affected columns: {', '.join(columns)}. "
        f"Sample rows: {preview}"
    )


def _raise_if_nan_series(series, labels, run_name, artifact_name):
    failing = [label for label in labels if pd.isna(series[label])]
    if not failing:
        return
    raise ValueError(
        f"NaN detected in {artifact_name} for run {run_name}. "
        f"Affected fields: {', '.join(failing)}"
    )


def _cached_run_matches_terms(data_file, run_name, expected_terms, output_label=None):
    coef_file = coefs_path(data_file, _artifact_run_name(run_name), output_label=output_label)
    loglik_file = loglik_path(data_file, _artifact_run_name(run_name), output_label=output_label)

    if not coef_file.exists() or not loglik_file.exists():
        return False

    ests = pd.read_csv(coef_file)
    if {'estimate', 'std.error', 'statistic', 'p.value'} - set(ests.columns):
        return False
    if ests[['estimate', 'std.error', 'statistic', 'p.value']].isna().any().any():
        return False

    cached_terms = set(ests['term'].dropna().unique())
    if cached_terms != set(expected_terms):
        return False

    loglik = pd.read_csv(loglik_file)
    required_stats = ['logLik', 'null.logLik', 'AIC', 'BIC', 'r.squared', 'adj.r.squared', 'nobs']
    if set(required_stats) - set(loglik.columns):
        return False
    if loglik[required_stats].isna().any().any():
        return False

    return True


def _artifact_run_name(run_name):
    return str(run_name).replace(":", "_")


def subprocess_rscript(
    work_dir,
    datafile,
    num_runs,
    run_name,
    last_integer,
    variables,
    seed=replication_seed,
    output_label=None,
):
    """
    Runs rscript in a subprocess
    :param work_dir:
    :param datafile:
    :param num_runs:
    :param run_name:
    :param last_integer:
    :param variables:
    :return:
    """
    artifact_run_name = _artifact_run_name(run_name)
    temp_run_name = artifact_run_name if output_label is None else f"tmp__{output_label}__{artifact_run_name}"

    command = [
        str(rdir),
        "--vanilla",
        "discrete_choice.R",
        str(work_dir),
        datafile,
        num_runs,
        temp_run_name,
        last_integer,
        str(seed),
    ] + variables

    # Execute the command
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # logging.info("Starting R script ...")
    print("Starting R script ...")

    while True:
        output = process.stdout.readline()
        if output == "" and process.poll() is not None:
            break
        if output:
            #logging.info(output.strip())
            print(output.strip())

    stderr = process.stderr.read()
    return_code= process.returncode
    if return_code != 0:
        # logging.info(f"Error: {stderr.strip()}")
        raise RuntimeError(f"R script failed with error: {stderr}")
    else:
        # logging.info("R script finished successfully")
        print("R script finished successfully")

    if output_label is not None:
        source_coef = Path(work_dir) / "data" / "results" / f"spatialdc_coefs_{datafile}_{temp_run_name}.csv"
        source_loglik = Path(work_dir) / "data" / "results" / f"spatialdc_loglik_{datafile}_{temp_run_name}.csv"
        target_coef = coefs_path(datafile, artifact_run_name, output_label=output_label)
        target_loglik = loglik_path(datafile, artifact_run_name, output_label=output_label)
        source_coef.replace(target_coef)
        source_loglik.replace(target_loglik)


def postprocess_spatialdc(work_dir, datafile, run_name, output_label=None):
    artifact_run_name = _artifact_run_name(run_name)
    ests = pd.read_csv(coefs_path(datafile, artifact_run_name, output_label=output_label))
    _raise_if_nan(ests, ['estimate', 'std.error', 'statistic', 'p.value'], run_name, 'raw coefficient output')
    n = len(ests.groupby('mod').count())

    if ests['term'].str.contains('scalePar').any():
        ests.loc[ests['term'].str.contains('scalePar'), 'term'] = 'lcoe'

    coefs = ests[['term', 'estimate']].groupby(['term']).mean()
    coefs['std'] = ests[['term', 'estimate']].groupby(['term']).std()
    coefs['tstat'] = coefs['estimate'] / coefs['std']
    coefs['pval'] = t.sf(np.abs(coefs['tstat']), n - 1)
    if 'min_lcoe' in coefs.index:
        coefs['valuation'] = coefs['estimate'] / coefs.loc['min_lcoe', 'estimate']
    else:
        raise ValueError(f"min_lcoe coefficient missing for run {run_name}")
    coefs.loc[coefs['pval'] <= 0.01, 'significance'] = '***'
    coefs.loc[(coefs['pval'] > 0.01) & (coefs['pval'] <= 0.05), 'significance'] = '**'
    coefs.loc[(coefs['pval'] > 0.05) & (coefs['pval'] <= 0.1), 'significance'] = '*'
    coefs.loc[coefs['pval'] > 0.1, 'significance'] = ''
    _raise_if_nan(coefs.reset_index(), ['estimate', 'std', 'tstat', 'pval', 'valuation'], run_name, 'processed coefficient summary')

    loglik = pd.read_csv(loglik_path(datafile, artifact_run_name, output_label=output_label))
    _raise_if_nan(loglik, ['logLik', 'null.logLik', 'AIC', 'BIC', 'r.squared', 'adj.r.squared', 'nobs'], run_name, 'raw model statistics')
    loglik = loglik.mean()
    _raise_if_nan_series(loglik, ['logLik', 'null.logLik', 'AIC', 'BIC', 'r.squared', 'adj.r.squared', 'nobs'], run_name, 'mean model statistics')

    return coefs, loglik


def leave_one_out(
    variables,
    integer_variables,
    name,
    work_dir,
    data_file,
    num_runs,
    compute=True,
    protected_variables=None,
    output_label=None,
):
    """
    Iteratively runs spatial discrete choice leaving one variable or interaction term out
    :param compute:
    :param variables:
    :param integer_variables:
    :param name:
    :param work_dir:
    :param data_file:
    :param num_runs:
    :return:
    """
    # generate all leave-one-out variable combinations
    protected_variables = set() if protected_variables is None else set(protected_variables)

    iter_dict = {}
    removable_variables = [var for var in variables if var not in protected_variables]
    for i, var in enumerate(removable_variables):
        reduced_vars = [v for v in variables if v != var]
        iter_dict[i] = [var, reduced_vars]

    # run each leave-one-out combination
    if compute:
        for i in range(len(iter_dict)):
            var = iter_dict[i][0]
            reduced_vars = iter_dict[i][1]
            run_name = f"{name}_{var}"
            last_int = str(len([ivar for ivar in integer_variables if ivar != var]))
            artifact_name = _artifact_run_name(run_name)
            cached = _cached_run_matches_terms(
                data_file,
                artifact_name,
                reduced_vars,
                output_label=output_label,
            )
            if not cached:
                subprocess_rscript(
                    work_dir,
                    data_file,
                    num_runs,
                    artifact_name,
                    last_int,
                    reduced_vars,
                    output_label=output_label,
                )

    # process results of runs
    run_stats = []
    run_coeffs = []
    for i in range(len(iter_dict)):
        var = iter_dict[i][0]
        run_name = f"{name}_{var}"
        coeff, loglik = postprocess_spatialdc(
            work_dir,
            data_file,
            _artifact_run_name(run_name),
            output_label=output_label,
        )
        run_stat = pd.concat([pd.Series([var]), loglik], ignore_index=True)
        run_stats.append(run_stat)
        run_coeffs.append(coeff)

    selection_criteria = pd.concat(run_stats, axis=1).T.reset_index(drop=True)
    selection_criteria.columns = (["variable"] + list(loglik.index))

    return selection_criteria, run_coeffs


def calculate_num_turbines(num_turbines, lcoe_map, lcoe_thresh, space_px, max_turbines, min_turbines):
    if num_turbines == 'auto':
        nt = np.ceil((lcoe_map < lcoe_thresh).sum().sum() / 1000) * 1000
        ntmax = np.floor(lcoe_map.size / space_px**3 * 0.5 / 1000) * 1000
        nturb = max(min(nt, ntmax, max_turbines), min_turbines)
    else:
        nturb = num_turbines
    return nturb


def _extract_values_at_locations(array, l_idx, b_idx):
    """Extract pairwise raster values without materializing an NxN intermediate."""
    return array.to_numpy()[l_idx, b_idx]


def locations_to_gdf(cost_array, locations, cost_name="LCOE", energy_array=None, power_array=None):
    if cost_array.chunks is not None:
        cost_array = cost_array.compute()

    l_idx = locations["l"].to_numpy(dtype=np.int64)
    b_idx = locations["b"].to_numpy(dtype=np.int64)
    cost_gdf = gpd.GeoDataFrame(
        data={cost_name: _extract_values_at_locations(cost_array, l_idx, b_idx)},
        geometry=gpd.points_from_xy(cost_array["x"].to_numpy()[b_idx], cost_array["y"].to_numpy()[l_idx]),
        crs=cost_array.rio.crs,
    )
    if energy_array is not None:
        if energy_array.chunks is not None:
            energy_array = energy_array.compute()
        cost_gdf["Energy"] = _extract_values_at_locations(energy_array, l_idx, b_idx)

    if power_array is not None:
        if power_array.chunks is not None:
            power_array = power_array.compute()
        cost_gdf["Power"] = _extract_values_at_locations(power_array, l_idx, b_idx)

    return cost_gdf.sort_values(by=cost_name).reset_index(drop=True)
