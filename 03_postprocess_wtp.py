# %% imports
import logging
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns
import cleo

from scow.utils import postprocess_spatialdc
from config import repo, data_ver

# %% settings
sens_name = "" if data_ver == "2014" else "_modern"

with open(f"parsimonious_name{sens_name}.txt", "r") as f:
    parsimonious_name = f.read()

# %% load estimated model
parsimonious_coeffs, parsimonious_stats = postprocess_spatialdc(repo, data_ver, parsimonious_name)
parsimonious_coeffs["valuation"] = parsimonious_coeffs["estimate"] / parsimonious_coeffs.loc["min_lcoe", "estimate"]

# %% load site data
atlas = cleo.Atlas(repo, "AUT",  crs="epsg:31287").load(region='Niederösterreich', scenario='2014', timestamp='20241210T140021')

if not atlas.wind_turbines:
    atlas.wind_turbines = [
        "Enercon.E40.500", "Enercon.E82.3000", "Enercon.E101.3050", "Enercon.E115.3000",
        "Vestas.V100.1800", "Vestas.V100.2000", "Vestas.V112.3075"]

if 'lcoe' not in atlas.wind.data.data_vars:
    atlas.wind.compute_lcoe(turbine_cost_share=0.7)
if 'min_lcoe' not in atlas.wind.data.data_vars:
    atlas.wind.minimum_lcoe()
if 'min_lcoe' not in atlas.landscape.data.data_vars:
    atlas.landscape.add(atlas.wind.data.min_lcoe, name='min_lcoe')

# generate interaction data
for var in list(parsimonious_coeffs.index):
    if ":" in var:
        vars = var.split(":")
        atlas.landscape.data[var] = atlas.landscape.data[vars[0]] * atlas.landscape.data[vars[1]]

atlas.landscape.data = atlas.landscape.data.merge(atlas.wind.data["min_lcoe"].compute())
spatial_data = atlas.landscape.data[list(parsimonious_coeffs.index)]

# %% compute social cost
soco = xr.Dataset()
for var in list(parsimonious_coeffs.index):
    soco[var] = spatial_data[var] * parsimonious_coeffs.loc[var, "valuation"]
    if var != "min_lcoe":
        if "Local Social Cost" not in soco.data_vars:
            soco["Local Social Cost"] = soco[var]
        else:
            soco["Local Social Cost"] = soco[var] + soco["Local Social Cost"]

if soco["Local Social Cost"].min() < 0:
    soco["Local Social Cost"] -= soco["Local Social Cost"].min()

soco["Local Social Cost"].attrs['unit'] = "€/MWh"
soco.to_netcdf(atlas.path / 'data' / 'results' / f'soco{sens_name}.nc')

# %% Plot map of local social cost and LCOE
lcoe = soco['min_lcoe']
lcoe.name = "Quasi-LCOE"
lcoe.attrs['unit'] = '€/MWh'

figures_dir = atlas.path / 'figures'
figures_dir.mkdir(parents=True, exist_ok=True)

fig, axs = plt.subplots(figsize=(7, 5))
soco["Local Social Cost"].plot(ax=axs, robust=True)
axs.set_title("")
axs.axis('off')
axs.set_title("")
axs.axis('off')
plt.tight_layout()
plt.savefig(figures_dir / f'map_soco_{data_ver}.png', dpi=300)
plt.close()

fig, axs = plt.subplots(figsize=(7, 5))
lcoe.plot(ax=axs, robust=True)
axs.set_title("")
axs.axis('off')
axs.set_title("")
axs.axis('off')
plt.tight_layout()
plt.savefig(figures_dir / f'map_lcoe_{data_ver}.png', dpi=300)
plt.close()

logging.info('Maps of Local Social Cost and LCOE plotted')

# %% Plot distributions of LCOE, local social cost, and total cost
cut_off = 0.999  # 0.999
cost_by_turbine = pd.DataFrame()

cost_distribution = pd.DataFrame(data=soco['Local Social Cost'].to_pandas().stack().reset_index(drop=True), columns=['Cost [€/MWh]'])
cost_distribution['target'] = 0
cost_distribution['Cost Component'] = f'Local Social Cost'

lcoe_at_all = pd.DataFrame(data=soco[f'min_lcoe'].to_pandas().stack().reset_index(drop=True), columns=['Cost [€/MWh]'])
lcoe_at_all['target'] = 1
lcoe_at_all['Cost Component'] = f'Quasi-LCOE'
lcoe_at_all = lcoe_at_all.loc[lcoe_at_all['Cost [€/MWh]'] <= lcoe_at_all['Cost [€/MWh]'].quantile(cut_off), :]
cost_distribution = pd.concat([cost_distribution, lcoe_at_all])

loco = soco['Local Social Cost'] + soco[f'min_lcoe']
loco = pd.DataFrame(data=loco.to_pandas().stack().reset_index(drop=True), columns=['Cost [€/MWh]'])
loco_at_all = loco.loc[loco['Cost [€/MWh]'] <= loco['Cost [€/MWh]'].quantile(cut_off), :]
loco_at_all['target'] = 2
loco_at_all['Cost Component'] = f'Total Local Social Cost'
cost_distribution = pd.concat([cost_distribution, loco_at_all])
cost_by_turbine = pd.concat([cost_by_turbine, cost_distribution])

# %% compute summary stats of each distribution
destat_qlcoe = cost_by_turbine.loc[cost_by_turbine["Cost Component"] == "Quasi-LCOE", "Cost [€/MWh]"].describe()
destat_loco = cost_by_turbine.loc[cost_by_turbine["Cost Component"] == "Local Social Cost", "Cost [€/MWh]"].describe()
destat_soco = cost_by_turbine.loc[cost_by_turbine["Cost Component"] == "Total Local Social Cost", "Cost [€/MWh]"].describe()


# %% generate distributions plot
fig, ax = plt.subplots(figsize=(8, 5))
sns.kdeplot(data=cost_by_turbine, x='Cost [€/MWh]', hue='Cost Component', fill=True)
plt.grid(alpha=0.4)
ax.set_axisbelow(True)
plt.tight_layout()
plt.savefig(atlas.path / f'figures/dist_soco_{data_ver}.png', dpi=300)
plt.close()

logging.info('Cost distributions plotted')

# %% analyse social cost components
import numpy as np
qus = [0.9, 0.95, 0.99]
soqu = pd.DataFrame(index=list(soco.data_vars), columns=qus)
for var in soco.data_vars:
    for q in qus:
        soqu.loc[var, q] = np.abs(soco[var]).quantile(q).values


#  plot maps of all social cost components
num_vars = len(soco.data_vars)
# Calculate the grid size for the subplots
num_cols = 4
num_rows = (num_vars + num_cols - 1) // num_cols
# Create a figure and axes for the subplots
fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(26, num_rows * 5))
# Flatten the axes array for easy iteration
axes = axes.flatten()
# Plot each data variable in a subplot
for i, var in enumerate(soco.data_vars):
    soco[var].plot(ax=axes[i])
    axes[i].set_title(var)
    axes[i].tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)  # Disable ticks

# Remove any unused subplots
for j in range(i+1, num_rows * num_cols):
    fig.delaxes(axes[j])
# Adjust layout
plt.tight_layout()
plt.savefig(atlas.path / f'figures/soco_facet{sens_name}.png', dpi=300)
plt.close()
logging.info("Facet of social cost variables plotted")
