# %% imports
import logging
import numpy as np
import pandas as pd
import xarray as xr
from scow.optimization import solve_objective_chunked
from scow.paths import opt_cost_path, opt_locations_path, soco_path
from scow.stages import load_atlas_for_analysis
from scow.utils import locations_to_gdf
from config import atlas_path, data_ver, spacing, nturbines

slice_count = 6
read_locations = False


def load_energy_for_stage4(atlas, power):
    """Return annual energy in GWh, deriving it when CLEO `optimal_energy` is unusable."""
    energy = atlas.wind.data.optimal_energy.compute()
    if not energy.isnull().all():
        return energy

    lcoe = atlas.wind.data.lcoe.compute()
    capacity_factors = atlas.wind.data.capacity_factors.compute()
    valid = np.isfinite(lcoe).any(dim="turbine")
    argmin = lcoe.fillna(np.inf).argmin(dim="turbine")

    selected_cf = xr.apply_ufunc(
        lambda cf_vals, idx_vals: np.take_along_axis(
            cf_vals,
            idx_vals[None, ...].astype(np.int64),
            axis=0,
        )[0],
        capacity_factors,
        argmin,
        input_core_dims=[["turbine", "y", "x"], ["y", "x"]],
        output_core_dims=[["y", "x"]],
        vectorize=False,
        dask="allowed",
        output_dtypes=[capacity_factors.dtype],
    )
    # `optimal_power` is in kW and downstream SCOW accounting expects annual energy in GWh.
    return (power * selected_cf * 8760 / 1_000_000).where(valid).rename("optimal_energy")


def objective_outputs_complete(locations_path, cost_path):
    """Return True only when both stage-4 outputs exist and contain rows."""
    if not locations_path.exists() or not cost_path.exists():
        return False

    try:
        locations = pd.read_csv(locations_path, index_col=0)
        cost = pd.read_csv(cost_path, index_col=0)
    except Exception:
        return False

    return (not locations.empty) and (not cost.empty)

# %%
atlas = load_atlas_for_analysis(
    {"atlas_path": atlas_path, "scenario": data_ver},
    require_optimal_power=True,
    require_optimal_energy=True,
)

power = atlas.wind.data.optimal_power.compute()
energy = load_energy_for_stage4(atlas, power)

social_cost = xr.open_dataset(soco_path(data_ver))

objectives = ["loco", "lcoe", "soco"]

for objective in objectives:
    if objective == "lcoe":
        cost_array = social_cost["min_lcoe"].compute()
        cost_threshold = 75
    elif objective == "loco":
        cost_array = social_cost["Local Social Cost"].compute()
        cost_threshold = 100
    elif objective == "soco":
        cost_array = (social_cost["Local Social Cost"] + social_cost["min_lcoe"]).compute()
        cost_threshold = 125
    else:
        raise ValueError(f"Unknown stage-4 objective: {objective}")

    locations_path = opt_locations_path(data_ver, objective)
    cost_path = opt_cost_path(data_ver, objective)
    if objective_outputs_complete(locations_path, cost_path):
        logging.info("Skipping completed objective: %s", objective)
        continue

    if not read_locations:
        locations = solve_objective_chunked(
            cost_array,
            scenario=data_ver,
            objective=objective,
            num_slices=slice_count,
            num_turbines=nturbines,
            space_px=spacing,
            lcoe_thresh=cost_threshold,
            min_turbines=1200,
            max_turbines=6000,
            axis=0,
        )

        locations.to_csv(locations_path)
    else:
        locations = pd.read_csv(locations_path, index_col=[0])

    # %% uncover total social cost at optimal locations
    cost_at_best_locations = locations_to_gdf(
        cost_array,
        locations,
        cost_name=objective,
        energy_array=energy,
        power_array=power,
    )
    cost_at_best_locations = cost_at_best_locations.dropna(how="any", axis=0)
    cost_at_best_locations = cost_at_best_locations.sort_values(by=objective)

    cost_at_best_locations["CumSumEnergy"] = cost_at_best_locations["Energy"].cumsum() / 1000
    cost_at_best_locations.to_csv(cost_path)
    logging.info(f"Optimal cost for {objective} written to: "
                 f"{str(cost_path)}")
