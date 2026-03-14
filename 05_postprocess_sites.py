# %% Imports
import logging
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
import pandas as pd
import geopandas as gpd
from shapely import wkt
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
from scow.stages import load_atlas_for_analysis

# %% Settings
BL = 'Niederösterreich'
country = 'AUT'
model_spec = 'base'  # ['parsimonious', 'base']
scenario = 'cost'
year = 2014  # [2014, 'modern']
GENERATION_TARGET = 8.00
sens_name = "" if data_ver == "2014" else f"_{data_ver}"

# %% Read optimization results and maps
state = gpd.read_file(repo / 'data' / 'site' / 'vgd' / 'vgd_oesterreich.shp')
state = state[['BL', 'geometry']].dissolve(by='BL').reset_index()
state = state[state.BL == BL]

atlas = load_atlas_for_analysis({"atlas_path": atlas_path, "scenario": data_ver})
social_cost = xr.open_dataset(soco_path(data_ver), engine="netcdf4")
opt_cost_soco = pd.read_csv(opt_cost_path(data_ver, "soco"))
opt_cost_loco = pd.read_csv(opt_cost_path(data_ver, "loco"))
opt_cost_lcoe = pd.read_csv(opt_cost_path(data_ver, "lcoe"))


# %% Convert to GeoDataFrame
def to_geopandas(df):
    df.geometry = df.geometry.apply(wkt.loads)
    cols = [c for c in list(df.columns) if c not in ['geometry', 'Unnamed: 0']]
    df = gpd.GeoDataFrame(geometry=df['geometry'], data=df[cols])
    return df

opt_cost_soco = to_geopandas(opt_cost_soco)
opt_cost_loco = to_geopandas(opt_cost_loco)
opt_cost_lcoe = to_geopandas(opt_cost_lcoe)


# %%  Compute cost at locations minimizing social cost
def annual_cost(df):
    """
    Calculate the annual cost in million euros
    :param df:
    :return:
    """
    df['Annual Social Cost'] = df['soco'] * df['Energy'] / 1000
    df['Annual Local Social Cost'] = df['loco'] * df['Energy'] / 1000
    df['Annual Private Cost'] = df['lcoe'] * df['Energy'] / 1000
    return df


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
cost_comparison = pd.DataFrame(index=["Private Cost", "Social Cost"],
                               columns=["Wind Turbines", "Local Social Cost", "Private Cost", "Total Local Social Cost"])
# number of wind turbines
turbine_sites_lcoe = opt_cost_lcoe['CumSumEnergy'] <= GENERATION_TARGET
cost_comparison.loc["Private Cost", "Wind Turbines"] = turbine_sites_lcoe.sum()
cost_comparison.loc["Private Cost", "Private Cost"] = opt_cost_lcoe.loc[turbine_sites_lcoe, "Annual Private Cost"].sum()
cost_comparison.loc["Private Cost", "Local Social Cost"] = opt_cost_lcoe.loc[turbine_sites_lcoe, "Annual Local Social Cost"].sum()
cost_comparison.loc["Private Cost", "Total Local Social Cost"] = opt_cost_lcoe.loc[turbine_sites_lcoe, "Annual Social Cost"].sum()

turbine_sites_soco = opt_cost_soco['CumSumEnergy'] <= GENERATION_TARGET
cost_comparison.loc["Social Cost", "Wind Turbines"] = turbine_sites_soco.sum()
cost_comparison.loc["Social Cost", "Private Cost"] = opt_cost_soco.loc[turbine_sites_soco, "Annual Private Cost"].sum()
cost_comparison.loc["Social Cost", "Local Social Cost"] = opt_cost_soco.loc[turbine_sites_soco, "Annual Local Social Cost"].sum()
cost_comparison.loc["Social Cost", "Total Local Social Cost"] = opt_cost_soco.loc[turbine_sites_soco, "Annual Social Cost"].sum()

latex_coco = cost_comparison.to_latex(float_format="%.1f")
coco_lines = latex_coco.split('\n')[4:-3]
latex_coco_without_header_footer = '\n'.join(coco_lines)

latex_cost_comparison_path().write_text(latex_coco_without_header_footer)


# %% Plot maps of spatial wind turbine allocations
for df, filename in zip([opt_cost_lcoe, opt_cost_loco, opt_cost_soco],
                        ['map_turbines_lcoe.png', 'map_turbines_loco.png', 'map_turbines_soco.png']):
    fig, ax = plt.subplots(figsize=(8, 5))
    state.plot(ax=ax, alpha=0.3, edgecolor='k')
    df[df['CumSumEnergy'] <= GENERATION_TARGET].plot(ax=ax, color='red', markersize=2)
    try:
        cx.add_basemap(ax, crs=state.crs, source=cx.providers.BasemapAT.terrain, zoom=10)
    except Exception as exc:
        logging.warning("Basemap fetch failed for %s: %s", filename, exc)
    ax.set_axis_off()
    plt.tight_layout()
    objective = filename.removeprefix("map_turbines_").removesuffix(".png")
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
