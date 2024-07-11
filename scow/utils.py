# %% imports
import logging
import os
import subprocess
import numpy as np
import pandas as pd
import geopandas as gpd
import gams.transfer as gt
from pathlib import Path
from scipy.stats import t


def subprocess_rscript(work_dir, datafile, num_runs, run_name, last_integer, variables):
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
    command = [
        "C:/myprogs/R/R-4.2.2/bin/Rscript",
        "discrete_choice.R",
        str(work_dir),
        datafile,
        num_runs,
        run_name,
        last_integer,
    ] + variables

    # Execute the command
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # Wait for the process to finish
    stdout, stderr = process.communicate()
    # Check for errors
    if process.returncode != 0:
        print("Error:", stderr.decode())
    else:
        print("Output:", stdout.decode())


def postprocess_spatialdc(work_dir, datafile, run_name):
    ests = pd.read_csv(Path(work_dir) / 'data' / 'results' / f'spatialdc_coefs_{datafile}_{run_name}.csv')
    n = len(ests.groupby('mod').count())

    if ests['term'].str.contains('scalePar').any():
        ests.loc[ests['term'].str.contains('scalePar'), 'term'] = 'lcoe'

    coefs = ests[['term', 'estimate']].groupby(['term']).mean()
    coefs['std'] = ests[['term', 'estimate']].groupby(['term']).std()
    coefs['tstat'] = coefs['estimate'] / coefs['std']
    coefs['pval'] = t.sf(np.abs(coefs['tstat']), n - 1)
    if 'min_lcoe' in coefs.columns:
        coefs['valuation'] = coefs['estimate'] / coefs.loc['min_lcoe', 'estimate']
    else:
        coefs['valuation'] = np.nan
    coefs.loc[coefs['pval'] <= 0.01, 'significance'] = '***'
    coefs.loc[(coefs['pval'] > 0.01) & (coefs['pval'] <= 0.05), 'significance'] = '**'
    coefs.loc[(coefs['pval'] > 0.05) & (coefs['pval'] <= 0.1), 'significance'] = '*'
    coefs.loc[coefs['pval'] > 0.1, 'significance'] = ''

    loglik = pd.read_csv(Path(work_dir) / 'data' / 'results' / f'spatialdc_loglik_{datafile}_{run_name}.csv')
    loglik = loglik.mean()

    return coefs, loglik


def leave_one_out(variables, integer_variables, name, work_dir, data_file, num_runs, compute=True):
    """
    iteratively runs spatial discrete choice leaving one variable out
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
    iter_dict = {}
    for i, var in enumerate(variables):
        reduced_vars = [v for v in variables if var not in v]
        iter_dict[i] = [var, reduced_vars]

    # run each leave-one-out combination
    if compute:
        for i in range(0, len(iter_dict)):
            var = iter_dict[i][0]
            reduced_vars = iter_dict[i][1]
            run_name = f"{name}_{var}"
            last_int = str(len([ivar for ivar in integer_variables if ivar != var]))
            if not (Path(work_dir) / "data" / "results" / f"spatialdc_coefs_{data_file}_{run_name}.csv").exists():
                subprocess_rscript(work_dir, data_file, num_runs, run_name, last_int, reduced_vars)

    # process results of runs
    run_stats = []
    run_coeffs = []
    for i in range(0, len(iter_dict)):
        var = iter_dict[i][0]
        run_name = f"{name}_{var}"
        coeff, loglik = postprocess_spatialdc(work_dir, data_file, run_name.replace(":", "_"))
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


def clear_symbols(container):
    for sym in ['l', 'b', 'i', 'j', 'lcoe', 'num_turbines']:
        if container.hasSymbols(sym):
            container.removeSymbols(sym)


def add_symbols(container, lcoe_map, nturb, space_px):
    laenge = container.addSet('l', records=list(lcoe_map.index), description='laenge')
    breite = container.addSet('b', records=list(lcoe_map.columns), description='breite')
    abstd_l = container.addSet('i', records=list(range(0, space_px)), description='abstand in laenge')
    abstd_b = container.addSet('j', records=list(range(0, space_px)), description='abstand in  breite')
    container.addParameter('lcoe', domain=[laenge, breite], records=lcoe_map.stack().reset_index())
    container.addParameter('num_turbines', domain=[], records=nturb)


def run_optimization(gams_dict, gdx_out, n):
    gms_exe_dir = gams_dict['gams_exe']
    gms_model = gams_dict['gams_model']
    subprocess.run(f'{gms_exe_dir}\\gams {gms_model} gdx={gdx_out} lo=3 o=nul --fnamelp=oploc{n}.lp')


def sliced_location_optimization(gams_dict, gams_container, lcoe_array, num_slices, num_turbines, space_px,
                                 min_turbines=1000, max_turbines=10000, lcoe_thresh=85, gdx_out_string='base',
                                 read_only=False, axis=0):

    locations = pd.DataFrame()
    lcoe_df = pd.DataFrame(data=lcoe_array.data)
    num_pixels = np.ceil(np.round(lcoe_array.count().data, 0) / num_slices)
    os.chdir(gams_dict['gdx_output'])

    dim_max = lcoe_array.sizes['x'] if axis == 0 else lcoe_array.sizes['y']

    j = 0
    for _ in range(num_slices):
        size = 0
        i = j
        while size <= num_pixels and i < dim_max:
            if axis == 0:
                lcoe_slice = lcoe_df.iloc[:, i]
            elif axis == 1:
                lcoe_slice = lcoe_df.iloc[i, :]
            else:
                raise ValueError('axis must be 0 or 1')
            lcoe_slice = lcoe_slice.dropna()
            i += 1
            size += len(lcoe_slice)

        if axis == 0:
            lcoe_map = lcoe_df.iloc[:, j:i]
        else:
            lcoe_map = lcoe_df.iloc[j:i, :]

        nturb = calculate_num_turbines(num_turbines, lcoe_map, lcoe_thresh, space_px, max_turbines, min_turbines)
        logging.info(f'Number of turbines: {nturb}')

        gdx_out = gams_dict['gdx_output'] / f'locations_{gdx_out_string}_{_}.gdx'
        if not read_only:
            clear_symbols(gams_container)
            add_symbols(gams_container, lcoe_map, nturb, space_px)
            gams_container.write(str(gams_dict['gdx_input']))
            run_optimization(gams_dict, gdx_out, _)

        results = gt.Container()
        results.read(str(gdx_out), 'build')
        locs = results.data['build'].records[['l', 'b', 'level']]
        locs = locs.loc[locs['level'] > 0]
        locations = pd.concat([locations, locs])
        j = i

    locations[['l', 'b']] = locations[['l', 'b']].astype(int)
    return locations


def locations_to_gdf(cost_array, locations, cost_name="LCOE", energy_array=None, power_array=None):
    if cost_array.chunks is not None:
        cost_array = cost_array.compute()

    cost_at_locations = cost_array[locations['l'], locations['b']]
    cost_gdf = gpd.GeoDataFrame(
        data=cost_at_locations.data.diagonal(),
        columns=[cost_name],
        geometry=gpd.points_from_xy(cost_at_locations.x, cost_at_locations.y),
        crs=cost_array.rio.crs
    )
    if energy_array is not None:
        if energy_array.chunks is not None:
            energy_array = energy_array.compute()
        cost_gdf['Energy'] = energy_array[locations['l'], locations['b']].data.diagonal()

    if power_array is not None:
        if power_array.chunks is not None:
            power_array = power_array.compute()
        cost_gdf['Power'] = power_array[locations['l'], locations['b']].data.diagonal()

    return cost_gdf.sort_values(by=cost_name).reset_index(drop=True)

