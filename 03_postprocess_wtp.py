# %% imports
import logging
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns

from scow.paths import (
    dist_soco_figure_path,
    map_lcoe_figure_path,
    map_soco_figure_path,
    parsimonious_name_path,
    soco_facet_figure_path,
    soco_path,
)
from scow.stages import load_atlas_for_analysis
from scow.utils import postprocess_spatialdc
from config import repo, atlas_path, data_ver


def _flat_finite(da):
    values = np.asarray(da).ravel()
    return values[np.isfinite(values)]


with open(parsimonious_name_path(data_ver), "r") as f:
    parsimonious_name = f.read()

# %% load estimated model
parsimonious_coeffs, parsimonious_stats = postprocess_spatialdc(repo, data_ver, parsimonious_name)
parsimonious_coeffs["valuation"] = parsimonious_coeffs["estimate"] / parsimonious_coeffs.loc["min_lcoe", "estimate"]

# %% load site data
atlas = load_atlas_for_analysis({"atlas_path": atlas_path, "scenario": data_ver})

# generate interaction data
landscape_data = atlas.landscape.data
for var in list(parsimonious_coeffs.index):
    if ":" in var:
        parts = var.split(":")
        landscape_data[var] = landscape_data[parts[0]] * landscape_data[parts[1]]

if "min_lcoe" not in landscape_data.data_vars:
    landscape_data = landscape_data.merge(atlas.wind.data["min_lcoe"].compute())
spatial_data = landscape_data[list(parsimonious_coeffs.index)]

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
soco.to_netcdf(soco_path(data_ver))

# %% Plot map of local social cost and LCOE
lcoe = soco['min_lcoe']
lcoe.name = "Quasi-LCOE"
lcoe.attrs['unit'] = '€/MWh'

fig, axs = plt.subplots(figsize=(7, 5))
soco["Local Social Cost"].plot(ax=axs, robust=True)
axs.set_title("")
axs.axis('off')
plt.tight_layout()
plt.savefig(map_soco_figure_path(data_ver), dpi=300)
plt.close()

fig, axs = plt.subplots(figsize=(7, 5))
lcoe.plot(ax=axs, robust=True)
axs.set_title("")
axs.axis('off')
plt.tight_layout()
plt.savefig(map_lcoe_figure_path(data_ver), dpi=300)
plt.close()

logging.info('Maps of Local Social Cost and LCOE plotted')

# %% Plot distributions of LCOE, local social cost, and total cost
cut_off = 0.999  # 0.999
cost_by_turbine = pd.DataFrame()

loco_values = _flat_finite(soco["Local Social Cost"])
cost_distribution = pd.DataFrame(data=loco_values, columns=['Cost [€/MWh]'])
cost_distribution['target'] = 0
cost_distribution['Cost Component'] = f'Local Social Cost'

lcoe_values = _flat_finite(soco["min_lcoe"])
lcoe_at_all = pd.DataFrame(data=lcoe_values, columns=['Cost [€/MWh]'])
lcoe_at_all['target'] = 1
lcoe_at_all['Cost Component'] = f'Quasi-LCOE'
lcoe_at_all = lcoe_at_all.loc[lcoe_at_all['Cost [€/MWh]'] <= lcoe_at_all['Cost [€/MWh]'].quantile(cut_off), :]
cost_distribution = pd.concat([cost_distribution, lcoe_at_all], ignore_index=True)

total_local_social_cost = _flat_finite(soco["Local Social Cost"] + soco["min_lcoe"])
loco = pd.DataFrame(data=total_local_social_cost, columns=['Cost [€/MWh]'])
loco_at_all = loco.loc[loco['Cost [€/MWh]'] <= loco['Cost [€/MWh]'].quantile(cut_off), :]
loco_at_all['target'] = 2
loco_at_all['Cost Component'] = f'Total Local Social Cost'
cost_distribution = pd.concat([cost_distribution, loco_at_all], ignore_index=True)
cost_by_turbine = pd.concat([cost_by_turbine, cost_distribution], ignore_index=True)

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
plt.savefig(dist_soco_figure_path(data_ver), dpi=300)
plt.close()

logging.info('Cost distributions plotted')

# %% analyse social cost components
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
plt.savefig(soco_facet_figure_path(data_ver), dpi=300)
plt.close()
logging.info("Facet of social cost variables plotted")
