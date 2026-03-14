# %% imports
import json
import logging
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import xarray as xr
from scow.optimization import final_output_metadata, solve_objective_chunked
from scow.paths import opt_cost_path, opt_locations_path, opt_metadata_path, soco_path
from scow.stages import load_atlas_for_analysis
from scow.utils import locations_to_gdf
from config import atlas_path, data_ver, spacing, nturbines

slice_count = 6
MIN_TURBINES = 1200
MAX_TURBINES = 6000
# Retain the current stage-4 screening thresholds in EUR/MWh, but name them
# explicitly so the objective-specific assumptions are visible in one place.
OBJECTIVE_COST_THRESHOLDS_EUR_PER_MWH = {
    "loco": 100,
    "lcoe": 75,
    "soco": 125,
}


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


def objective_cost_array(social_cost, objective):
    """Return the stage-4 cost surface for one named objective."""
    if objective == "lcoe":
        return social_cost["min_lcoe"].compute()
    if objective == "loco":
        return social_cost["Local Social Cost"].compute()
    if objective == "soco":
        return (social_cost["Local Social Cost"] + social_cost["min_lcoe"]).compute()
    raise ValueError(f"Unknown stage-4 objective: {objective}")

# %%
atlas = load_atlas_for_analysis(
    {"atlas_path": atlas_path, "scenario": data_ver},
    require_optimal_power=True,
    require_optimal_energy=True,
)

power = atlas.wind.data.optimal_power.compute()
energy = load_energy_for_stage4(atlas, power)

social_cost = xr.open_dataset(soco_path(data_ver))

for objective, cost_threshold in OBJECTIVE_COST_THRESHOLDS_EUR_PER_MWH.items():
    cost_array = objective_cost_array(social_cost, objective)
    locations_path = opt_locations_path(data_ver, objective)
    cost_path = opt_cost_path(data_ver, objective)
    metadata_path = opt_metadata_path(data_ver, objective)

    locations = solve_objective_chunked(
        cost_array,
        scenario=data_ver,
        objective=objective,
        num_slices=slice_count,
        num_turbines=nturbines,
        space_px=spacing,
        lcoe_thresh=cost_threshold,
        min_turbines=MIN_TURBINES,
        max_turbines=MAX_TURBINES,
        axis=0,
    )

    locations.to_csv(locations_path)

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
    metadata = final_output_metadata(
        cost_array,
        scenario=data_ver,
        objective=objective,
        num_slices=slice_count,
        num_turbines=nturbines,
        space_px=spacing,
        lcoe_thresh=cost_threshold,
        min_turbines=MIN_TURBINES,
        max_turbines=MAX_TURBINES,
        axis=0,
    )
    metadata["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
    metadata["locations_file"] = locations_path.name
    metadata["cost_file"] = cost_path.name
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    logging.info(f"Optimal cost for {objective} written to: "
                 f"{str(cost_path)}")
