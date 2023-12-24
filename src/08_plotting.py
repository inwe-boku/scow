# %% imports
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely import wkt
import contextily as cx
import xarray as xr
from config import ROOTDIR
from logging_config import setup_logging

setup_logging()

# %% settings
BL = 'Niederösterreich'
country = 'AUT'
model_spec = 'base'  # 'parsimonious'
scenario = 'cost'
year = 2014
GENERATION_TARGET = 8.00

# %% read optimization results and maps
state = gpd.read_file(ROOTDIR / 'data/vgd/vgd_oesterreich.shp')
state = state[['BL', 'geometry']].dissolve(by='BL')
state.reset_index(inplace=True)
state = state.loc[state.BL == BL]

energy = xr.open_dataarray(ROOTDIR / f'data/preprocessed/energy_generation_{country}_{year}.nc')
energy = energy.drop_vars('turbine_models')
energy = energy.rio.reproject(state.crs)
energy_noe = energy.rio.clip(state.geometry, state.crs)

power = xr.open_dataarray(ROOTDIR / f'data/preprocessed/installed_power_{country}_{year}.nc')
power = power.rio.write_crs('epsg:3416')
power = power.rio.reproject(state.crs)
power_noe = power.rio.clip(state.geometry, state.crs)

social_cost = xr.open_dataset(ROOTDIR / f'data/results/social_cost_{scenario}_{year}_{model_spec}.nc')

optloc_soco = pd.read_csv(ROOTDIR / f'data/results/optimal_turbines_noe_6_13400_{scenario}_soco.csv')
optloc_lcoe = pd.read_csv(ROOTDIR / f'data/results/optimal_turbines_noe_6_13400_{scenario}_lcoe.csv')
# locations_lcoe = pd.read_csv(ROOTDIR / 'data/results/opt_locations_soco_min_lcoe.csv')

# %% convert to geopandas geodataframe
optloc_soco.geometry = optloc_soco.geometry.apply(wkt.loads)
optloc_soco = gpd.GeoDataFrame(geometry=optloc_soco['geometry'],
                               data=optloc_soco[['LCOE', 'Power', 'Energy', 'CumEnergy']])
optloc_soco = optloc_soco.rename(columns={'LCOE': 'Social Cost'})

optloc_lcoe.geometry = optloc_lcoe.geometry.apply(wkt.loads)
optloc_lcoe = gpd.GeoDataFrame(geometry=optloc_lcoe['geometry'],
                               data=optloc_lcoe[['LCOE', 'Power', 'Energy', 'CumEnergy']])

# %% update power, energy, and cum energy
optloc_soco['Energy'] = energy_noe.sel(x=optloc_soco.geometry.x.to_xarray(), y=optloc_soco.geometry.y.to_xarray(),
                                       method='nearest').to_pandas()
optloc_soco['Power'] = power_noe.sel(x=optloc_soco.geometry.x.to_xarray(), y=optloc_soco.geometry.y.to_xarray(),
                                     method='nearest').to_pandas()
optloc_soco['CumEnergy'] = optloc_soco['Energy'].cumsum() / 1000

optloc_lcoe['Energy'] = energy_noe.sel(x=optloc_lcoe.geometry.x.to_xarray(), y=optloc_lcoe.geometry.y.to_xarray(),
                                       method='nearest').to_pandas()
optloc_lcoe['Power'] = power_noe.sel(x=optloc_lcoe.geometry.x.to_xarray(), y=optloc_lcoe.geometry.y.to_xarray(),
                                     method='nearest').to_pandas()
optloc_lcoe['CumEnergy'] = optloc_lcoe['Energy'].cumsum() / 1000

# %% compute LCOE at locations minimizing social cost
optloc_soco['LCOE'] = social_cost['lcoe_min_lcoe'].sel(x=optloc_soco.geometry.x.to_xarray(),
                                                       y=optloc_soco.geometry.y.to_xarray(),
                                                       method="nearest").to_pandas()

optloc_soco['Local Cost'] = optloc_soco['Social Cost'] - optloc_soco['LCOE']

# %% compute total social cost at (a) socially optimal and (b) individually/LCOE optimal locations
optloc_lcoe['Local Social Cost'] = social_cost['Local Social Cost'].sel(x=optloc_lcoe.geometry.x.to_xarray(),
                                                                        y=optloc_lcoe.geometry.y.to_xarray(),
                                                                        method="nearest").to_pandas()
optloc_lcoe['Annual Local Social Cost'] = optloc_lcoe['Local Social Cost'] * optloc_lcoe['Energy'] * 1000
optloc_lcoe['Annual Private Cost'] = optloc_lcoe['LCOE'] * optloc_lcoe['Energy'] * 1000
optloc_lcoe.loc[0:1199, ['Annual Local Social Cost', 'Annual Private Cost', 'Energy']].sum()

optloc_soco['LCOE'] = social_cost['lcoe_min_lcoe'].sel(x=optloc_soco.geometry.x.to_xarray(),
                                                       y=optloc_soco.geometry.y.to_xarray(),
                                                       method="nearest").to_pandas()
optloc_soco['Annual Local Social Cost'] = optloc_soco['Social Cost'] * optloc_soco['Energy'] * 1000
optloc_soco['Annual Private Cost'] = optloc_soco['LCOE'] * optloc_soco['Energy'] * 1000

opt_alloc_soco = optloc_soco.loc[optloc_soco['CumEnergy'] <= GENERATION_TARGET, :]
opt_alloc_lcoe = optloc_lcoe.loc[optloc_lcoe['CumEnergy'] <= GENERATION_TARGET, :]

logging.info(f"Social cost optimal allocation: \n {opt_alloc_soco[['Annual Local Social Cost', 'Annual Private Cost']].sum() / 1000000}")
logging.info(f"Private cost optimal allocation: \n {opt_alloc_lcoe[['Annual Local Social Cost', 'Annual Private Cost']].sum() / 1000000}")
# optloc_soco.loc[0:1199, ['Annual Local Social Cost', 'Annual Private Cost', 'Energy']].sum()
# optloc_soco.loc[0:1585, ['Annual Local Social Cost', 'Annual Private Cost', 'Energy']].sum()

# %% plot maps of spatial wind turbine allocations
fig, ax = plt.subplots(figsize=(8, 5))
state.plot(ax=ax, alpha=0.3, edgecolor='k')
opt_alloc_lcoe.plot(ax=ax, color='red', markersize=2)
cx.add_basemap(ax, crs=state.crs, source=cx.providers.BasemapAT.terrain, zoom=10)  # source=cx.providers.BasemapAT.terrain, zoom=10)
ax.set_axis_off()
plt.tight_layout()
plt.savefig(ROOTDIR / 'figures/map_turbines_lcoe.png', dpi=200)

fig, ax = plt.subplots(figsize=(8, 5))
state.plot(ax=ax, alpha=0.3, edgecolor='k')
opt_alloc_soco.plot(ax=ax, color='red', markersize=2)
cx.add_basemap(ax, crs=state.crs, source=cx.providers.BasemapAT.terrain, zoom=10)
ax.set_axis_off()
plt.tight_layout()
plt.savefig(ROOTDIR / 'figures/map_turbines_soco.png', dpi=200)

logging.info('Maps of optimal wind turbine sites saved')

# %% plot cost curve
fig, ax = plt.subplots(figsize=(8, 5))
plt.stackplot(optloc_soco['CumEnergy'], optloc_soco['LCOE'], optloc_soco['Local Cost'], alpha=0.825)
ax.set_xlabel('Energy [TWh]')
ax.set_ylabel('Social Cost [€/MWh]')
ax.set_xlim([0, 60])
plt.grid()
ax.set_axisbelow(True)
plt.legend(['Quasi-LCOE', 'Local Social Cost'], loc='upper left')
plt.tight_layout()
plt.savefig(ROOTDIR / 'figures/cost_curve.png', dpi=200)
plt.close()

logging.info('Social cost curve saved')

# %% postprocess least-cost turbines
lc_turbines = xr.open_dataarray(ROOTDIR / f'data/preprocessed/least_cost_turbines_{country}_{year}.nc')
state_3416 = state.to_crs(lc_turbines.rio.crs)
lc_turbines_noe = lc_turbines.rio.clip(state_3416.geometry, state_3416.crs).to_pandas()
unique_lc_noe = lc_turbines_noe.stack().unique()

num_lc = pd.DataFrame(index=unique_lc_noe, columns=['count'])
for tt in unique_lc_noe:
    num_lc.loc[tt, 'count'] = (lc_turbines_noe == tt).sum().sum()
    print(f'{tt}: {(lc_turbines_noe == tt).sum().sum()}')
