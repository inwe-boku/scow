# %% Imports
import logging
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
import contextily as cx
import xarray as xr
from config import repo, atlas_path, data_ver
from scow.paths import (
    cost_curve_figure_path,
    latex_cost_comparison_path,
    map_turbines_figure_path,
    opt_cost_path,
    soco_path,
)
from scow.postprocess_sites_helpers import annual_cost, build_cost_comparison, sites_to_target, to_geopandas
from scow.stages import load_atlas_for_analysis

# %% Settings
BL = 'Niederösterreich'
GENERATION_TARGET = 8.00

# %% Read optimization results and maps
state = gpd.read_file(repo / 'data' / 'site' / 'vgd' / 'vgd_oesterreich.shp')
state = state[['BL', 'geometry']].dissolve(by='BL').reset_index()
state = state[state.BL == BL]

atlas = load_atlas_for_analysis({"atlas_path": atlas_path, "scenario": data_ver})
social_cost = xr.open_dataset(soco_path(data_ver), engine="netcdf4")
opt_cost_soco = pd.read_csv(opt_cost_path(data_ver, "soco"))
opt_cost_loco = pd.read_csv(opt_cost_path(data_ver, "loco"))
opt_cost_lcoe = pd.read_csv(opt_cost_path(data_ver, "lcoe"))

opt_cost_soco = to_geopandas(opt_cost_soco)
opt_cost_loco = to_geopandas(opt_cost_loco)
opt_cost_lcoe = to_geopandas(opt_cost_lcoe)


opt_cost_soco['lcoe'] = social_cost['min_lcoe'].sel(x=opt_cost_soco.geometry.x.to_xarray(),
                                                    y=opt_cost_soco.geometry.y.to_xarray(),
                                                    method="nearest").to_pandas()

opt_cost_soco['loco'] = opt_cost_soco['soco'] - opt_cost_soco['lcoe']
opt_cost_soco = annual_cost(opt_cost_soco)

# %% compute cost at locations minimizing lcoe
opt_cost_lcoe['loco'] = social_cost['Local Social Cost'].sel(x=opt_cost_lcoe.geometry.x.to_xarray(),
                                                             y=opt_cost_lcoe.geometry.y.to_xarray(),
                                                             method='nearest').to_pandas()

opt_cost_lcoe['soco'] = opt_cost_lcoe['lcoe'] + opt_cost_lcoe['loco']
opt_cost_lcoe = annual_cost(opt_cost_lcoe)

# %% compute cost at locations minimizing local social cost
opt_cost_loco['lcoe'] = social_cost['min_lcoe'].sel(x=opt_cost_loco.geometry.x.to_xarray(),
                                                    y=opt_cost_loco.geometry.y.to_xarray(),
                                                    method='nearest').to_pandas()

opt_cost_loco['soco'] = opt_cost_loco['lcoe'] + opt_cost_loco['loco']
opt_cost_loco = annual_cost(opt_cost_loco)

# %% Compute total social cost
cost_comparison = build_cost_comparison(opt_cost_lcoe, opt_cost_soco, GENERATION_TARGET)

latex_coco = cost_comparison.to_latex(float_format="%.1f")
coco_lines = latex_coco.split('\n')[4:-3]
latex_coco_without_header_footer = '\n'.join(coco_lines)

latex_cost_comparison_path().write_text(latex_coco_without_header_footer)


# %% Plot maps of spatial wind turbine allocations
for objective, df in zip(["lcoe", "loco", "soco"], [opt_cost_lcoe, opt_cost_loco, opt_cost_soco]):
    fig, ax = plt.subplots(figsize=(8, 5))
    state.plot(ax=ax, alpha=0.3, edgecolor='k')
    df.loc[sites_to_target(df, GENERATION_TARGET)].plot(ax=ax, color='red', markersize=2)
    try:
        cx.add_basemap(ax, crs=state.crs, source=cx.providers.BasemapAT.terrain, zoom=10)
    except Exception as exc:
        logging.warning("Basemap fetch failed for %s: %s", objective, exc)
    ax.set_axis_off()
    plt.tight_layout()
    plt.savefig(map_turbines_figure_path(objective), dpi=300)
    plt.close()

logging.info('Maps of optimal wind turbine sites saved')

# %% Plot social cost curve
fig, ax = plt.subplots(figsize=(8, 5))
plt.stackplot(opt_cost_soco['CumSumEnergy'], opt_cost_soco['lcoe'], opt_cost_soco['loco'], alpha=0.825)
ax.set_xlabel('Energy [TWh]')
ax.set_ylabel('Social Cost [€/MWh]')
ax.set_xlim([0, 60])
plt.grid()
ax.set_axisbelow(True)
plt.legend(['Quasi-LCOE', 'Local Social Cost'], loc='upper left')
plt.tight_layout()
plt.savefig(cost_curve_figure_path(), dpi=300)
plt.close()

logging.info('Social cost curve saved')
