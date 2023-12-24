# %% imports
import logging

import pandas as pd
import geopandas as gpd
import xarray as xr
import gams.transfer as gt

from config import ROOTDIR, gams_conf, country, nslices, nturbines, min_turb, max_turb, spacing
from src.vala import sliced_location_optimization, locations_to_gdf
from logging_config import setup_logging

setup_logging()

# %% input & data
sc_name = ['cost']
years = ['2014']  # , '2016']
mnames = ['base']  # , 'full']
campaigns = [f'{sc}_{yr}_{mn}' for sc in sc_name for yr in years for mn in mnames]

read_locations = True
nslices = 6

for camp in campaigns:
    scenario, year, model = camp.split('_')

    austria = gpd.read_file(ROOTDIR / 'data/vgd/vgd_oesterreich.shp')
    austria = austria[['BL', 'geometry']].dissolve(by='BL')
    austria.reset_index(inplace=True)
    noe = austria.loc[austria.BL == 'Niederösterreich', :]

    energy = xr.open_dataarray(ROOTDIR / f'data/preprocessed/energy_generation_{country}_{year}.nc')
    energy = energy.drop_vars('turbine_models')
    energy = energy.rio.reproject(austria.crs)
    energy_noe = energy.rio.clip(noe.geometry, noe.crs)

    power = xr.open_dataarray(ROOTDIR / f'data/preprocessed/installed_power_{country}_{year}.nc')
    power = power.rio.write_crs('epsg:3416')
    power = power.rio.reproject(austria.crs)
    power_noe = power.rio.clip(noe.geometry, noe.crs)

# %% optimal turbine allocation
    social_cost = xr.open_dataset(ROOTDIR / f'data/results/social_cost_{camp}.nc')

# %% commence computations
    objectives = ['soco']  # , 'lcoe']
    for obj in objectives:
        if obj == 'soco':
            tcost_array = social_cost['lcoe_min_lcoe'] + social_cost['Local Social Cost']
        elif obj == 'lcoe':
            tcost_array = social_cost['lcoe_min_lcoe']
        else:
            print('objective misspecified')
        tcost_array = tcost_array.rio.reproject(austria.crs)

        energy_noe = energy_noe.interp_like(tcost_array)
        power_noe = power_noe.interp_like(tcost_array)

        if read_locations is False:
            gams_transfer_container = gt.Container()
            locations = sliced_location_optimization(gams_conf, gams_transfer_container, tcost_array, min_turbines=1200,
                                                 max_turbines=6000, lcoe_thresh=120, num_slices=6, space_px=spacing,
                                                 num_turbines=nturbines, gdx_out_string=obj, read_only=False)
            locations.to_csv(ROOTDIR / f'data/results/opt_locations_{camp}_{obj}.csv')
        else:
            locations = pd.read_csv(ROOTDIR / f'data/results/opt_locations_{camp}_{obj}.csv', index_col=[0])

        nturbines = len(locations)

    # %% uncover total social cost at optimal locations
        soco = [locations_to_gdf(tcost_array, locations[i*nturbines:(i+1)*nturbines], energy_array=energy_noe,
                                 power_array=power_noe) for i in range(0, nslices)]
        soco = pd.concat(soco)
        soco = soco.dropna(how='any', axis=0)
        soco = soco.sort_values(by='LCOE')

        soco['CumEnergy'] = soco['Energy'].cumsum() / 1000
        fname = ROOTDIR / f'data/results/optimal_turbines_noe_{nslices}_{nturbines}_{scenario}_{obj}.csv'
        soco.to_csv(fname)
        logging.info(f'Optimal wind turbine sites determined and written to {fname}')

"""
# %% uncover LCOE at optimal locations
lcoe = social_cost['lcoe']
lcoe = lcoe.rio.reproject(austria.crs)

lcoe_at_soco = gpd.GeoDataFrame(data=soco['CumEnergy'], geometry=soco.geometry, crs=austria.crs)
lcoe_at_soco['LCOE'] = lcoe.sel(y=lcoe_at_soco.geometry.y.to_xarray(), x=lcoe_at_soco.geometry.x.to_xarray()).data

# %% uncover local social cost at optimal locations
loco = social_cost['Total Social Cost']
loco = loco.rio.reproject(austria.crs)

loco_at_soco = gpd.GeoDataFrame(data=soco['CumEnergy'], geometry=soco.geometry, crs=austria.crs)
loco_at_soco['Local Social Cost'] = loco.sel(y=loco_at_soco.geometry.y.to_xarray(), x=loco_at_soco.geometry.x.to_xarray()).data

# %% plot social cost supply curve and LCOE supply curve at optimal locations
fig, ax = plt.subplots(1, 1)
ax.plot(soco['CumEnergy'], soco['LCOE'])
# ax.plot(lcoe_at_soco['CumEnergy'], lcoe_at_soco['LCOE'])
# plt.title("Wind Power Supply Curve, Lower Austria")
plt.ylabel('Estimated Private + Local Social Cost [€/MWh]')
plt.xlabel('Expected Generation [TWh/a]')
# plt.xlim([0, 100])
plt.grid()
plt.savefig(ROOTDIR / 'doc/figures/supply_curve_noe')

# %% plot distributions of LCOE versus social cost at optimal locations
cost_distribution = pd.DataFrame(data=soco['LCOE'].values, columns=['Cost'])
cost_distribution['target'] = 0
cost_distribution['Cost Component'] = 'Total Social Cost'

lcoe_at_soco['target'] = 1
lcoe_at_soco['Cost Component'] = 'LCOE'
cost_distribution = pd.concat([cost_distribution, lcoe_at_soco[['LCOE', 'target', 'Cost Component']].rename({'LCOE': 'Cost'}, axis=1)])

loco_at_soco['target'] = 2
loco_at_soco['Cost Component'] = 'Local Social Cost'
cost_distribution = pd.concat([cost_distribution, loco_at_soco[['Local Social Cost', 'target', 'Cost Component']].rename({'Local Social Cost': 'Cost'}, axis=1)])

sns.displot(data=cost_distribution, x='Cost', hue='Cost Component', kind='kde', fill=True,
            palette=sns.color_palette('bright')[:3], height=5, aspect=1.5)
plt.savefig(ROOTDIR / 'figures/cost_distribution.png', dpi=200)

# %% plot distributions of social cost at actual locations versus social cost at optimal locations
"""
"""
# %% get private cost at locations
social_cost = xr.open_dataset(ROOTDIR / 'data/results/total_social_cost.nc')

socomp = {}
for component in social_cost:
    cmptmp = social_cost[component]
    cmptmp = cmptmp.rio.reproject(tcost_array.rio.crs)
    cmptmp = cmptmp.interp_like(tcost_array)
    socomp[component] = locations_to_gdf(cmptmp, locations, energy_array=energy_noe, power_array=power_noe)
    socomp[component].rename({'LCOE': component}, axis=1, inplace=True)
    if component == list(social_cost.keys())[0]:
        df = socomp[component]
    else:
        df = df.merge(socomp[component][['geometry', component]], how='left', on='geometry')

df = df.dropna(axis=0, how='any')
df['Total Social Cost'] = df[['private cost', 'disamenities (loglin)', 'ecological cost', 'interdependency']].sum(axis=1)
df = df.sort_values('Total Social Cost')
df['CumEnergy'] = df['Energy'].cumsum() / 1000

# %%
import matplotlib.pyplot as plt

plt.stackplot(df['CumEnergy'], df['private cost'], df['disamenities (loglin)'], df['ecological cost'], df['interdependency'])

# %% plot social cost and LCOE

# %
goco = soco.to_crs('epsg:4326')
# turbine in standing waters (Neusiedler See),
goco['x'] = goco.geometry.x
goco['y'] = goco.geometry.y
goco[['y', 'x']].to_csv(ROOTDIR / f'data/results/turbine_locations_noe_{nslices}_{nturbines}.csv')
"""
