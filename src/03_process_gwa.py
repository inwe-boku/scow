# %% imports
import os
import logging
import numpy as np
import pandas as pd
import xarray as xr

from config import ROOTDIR, turbines, TURBINE_YEARS, country
from config import BIAS_CORRECTION, OM_FIXED, OM_VARIABLE, DISCOUNT_RATE, LIFETIME
from src.vala import simulate_capacity_factors, compute_lcoe
# from src.vala import turbine_overnight_cost, grid_connect_cost, levelized_cost
from logging_config import setup_logging

# %% calculate power per pixel in LCOE optimum
if __name__ == "__main__":
    setup_logging()

    dir_results = ROOTDIR / 'data/preprocessed'
    if not os.path.exists(dir_results):
        os.mkdir(dir_results)

    logging.info('Processing Global Wind Atlas')
    logging.info(f'Parameters: Bias Correction Factor={BIAS_CORRECTION}, fixed O&M={OM_FIXED}, var O&M={OM_VARIABLE}, '
                 f'discount rate={DISCOUNT_RATE}, lifetime={LIFETIME}')

    for turbine_year in TURBINE_YEARS:
        turbine_sel = {t: val for t, val in turbines.items() if val[3] <= turbine_year}
        power_curves = pd.read_csv(ROOTDIR / 'data/preprocessed/powercurves.csv', sep=";", decimal=',')
        power_curves.set_index('speed', drop=True, inplace=True)

        cap_factors = simulate_capacity_factors(power_curves, turbine_sel, ROOTDIR, country, BIAS_CORRECTION)

        lcoe = compute_lcoe(cap_factors, power_curves, turbine_sel, BIAS_CORRECTION, OM_FIXED, OM_VARIABLE, DISCOUNT_RATE, LIFETIME)
        lcoe.to_netcdf(dir_results / f'lcoe_{country}_{turbine_year}.nc',
                       format='NETCDF4', engine='netcdf4')
        logging.info(f'LCOE for model year {turbine_year} computed and exported')

        least_cost_turbines = lcoe.idxmin(dim='turbine_models')
        least_cost_turbines.to_netcdf(dir_results / f'least_cost_turbines_{country}_{turbine_year}.nc',
                                      format='NETCDF4', engine='netcdf4')
        least_cost_idx = lcoe.fillna(999).argmin(dim='turbine_models')
        least_cost_idx.to_netcdf(dir_results / f'least_cost_index_{country}_{turbine_year}.nc',
                                 format='NETCDF4', engine='netcdf4')

        power = least_cost_turbines.copy()
        for name, char in turbines.items():
            power = xr.where(power == name, char[0], power)
        power.data = power.data.astype('float')
        power.to_netcdf(path=ROOTDIR / f'data/preprocessed/installed_power_{country}_{turbine_year}.nc',
                        format='NETCDF4', engine='netcdf4')

        least_cost_capacity_factor = cap_factors.isel(turbine_models=least_cost_idx).drop_vars('turbine_models')
        least_cost_capacity_factor = least_cost_capacity_factor.assign_coords({'turbine_models': 'min_lcoe'})

        cap_facs = xr.concat([cap_factors, least_cost_capacity_factor], dim='turbine_models', fill_value=np.nan)
        cap_facs.name = 'capacity_factor'
        cap_facs = cap_facs.rio.write_crs('epsg:3416')
        cap_facs.to_netcdf(path=ROOTDIR / f'data/preprocessed/capacity_factors_{country}_{turbine_year}.nc',
                           format='NETCDF4', engine='netcdf4')

        # energy per pixel = capacity_factor * power * hours_in_year / 1000
        energy_giga = least_cost_capacity_factor * power * 8760 / 1000000
        energy_giga.to_netcdf(path=ROOTDIR / f'data/preprocessed/energy_generation_{country}_{turbine_year}.nc',
                              format='NETCDF4', engine='netcdf4')

    logging.info('Processing of GWA complete')
