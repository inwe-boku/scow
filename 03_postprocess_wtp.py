# %% imports
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from scow.postprocess_wtp_helpers import (
    build_interaction_layers,
    compute_social_cost_dataset,
    cost_distribution_frame,
)
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


with open(parsimonious_name_path(data_ver), "r") as f:
    parsimonious_name = f.read()

# %% load estimated model
parsimonious_coeffs, parsimonious_stats = postprocess_spatialdc(repo, data_ver, parsimonious_name)
parsimonious_coeffs["valuation"] = parsimonious_coeffs["estimate"] / parsimonious_coeffs.loc["min_lcoe", "estimate"]

# %% load site data
atlas = load_atlas_for_analysis({"atlas_path": atlas_path, "scenario": data_ver})

# generate interaction data
landscape_data = build_interaction_layers(atlas.landscape.data, list(parsimonious_coeffs.index))

if "min_lcoe" not in landscape_data.data_vars:
    landscape_data = landscape_data.merge(atlas.wind.data["min_lcoe"].compute())
spatial_data = landscape_data[list(parsimonious_coeffs.index)]

# %% compute social cost
soco = compute_social_cost_dataset(spatial_data, parsimonious_coeffs["valuation"])
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
cost_by_turbine = cost_distribution_frame(soco, cut_off)

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
