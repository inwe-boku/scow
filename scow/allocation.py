from __future__ import annotations

import json
from datetime import datetime, timezone
import logging
from pathlib import Path

import numpy as np
import xarray as xr

from scow.logging_utils import configure_script_logging
from scow.optimization import final_output_metadata, solve_objective_chunked
from scow.paths import opt_cost_path, opt_locations_path, opt_metadata_path
from scow.social_cost import load_social_cost_atlas
from scow.stages import load_atlas_for_analysis
from scow.utils import locations_to_gdf


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
    return (power * selected_cf * 8760 / 1_000_000).where(valid).rename("optimal_energy")


def objective_cost_array(social_cost, objective: str):
    """Return the current stage-4 cost surface for one named objective."""

    if objective == "lcoe":
        return social_cost["min_lcoe"].compute()
    if objective == "loco":
        return social_cost["Local Social Cost"].compute()
    if objective == "soco":
        return (social_cost["Local Social Cost"] + social_cost["min_lcoe"]).compute()
    raise ValueError(f"Unknown stage-4 objective: {objective}")


def load_allocation_inputs(*, atlas_file: Path, scenario: str):
    """Load the atlas-derived inputs consumed by stage 4."""

    atlas = load_atlas_for_analysis(
        {"atlas_path": atlas_file, "scenario": scenario},
        require_optimal_power=True,
        require_optimal_energy=True,
    )
    power = atlas.wind.data.optimal_power.compute()
    energy = load_energy_for_stage4(atlas, power)
    social_cost = load_social_cost_atlas(scenario)
    return power, energy, social_cost


def optimize_allocation(
    *,
    objective: str,
    scenario: str,
    spacing: int,
    num_turbines,
    power,
    energy,
    social_cost,
    cost_threshold: float,
    slice_count: int,
    min_turbines: int,
    max_turbines: int,
    logger: logging.Logger | logging.LoggerAdapter | None = None,
) -> None:
    """Solve one stage-4 siting objective and write the canonical artifacts."""

    if logger is None:
        logger = configure_script_logging("stage4")

    logger.info("Solving objective %s with threshold %s EUR/MWh", objective, cost_threshold)
    cost_array = objective_cost_array(social_cost, objective)
    locations_path = opt_locations_path(scenario, objective)
    cost_path = opt_cost_path(scenario, objective)
    metadata_path = opt_metadata_path(scenario, objective)

    locations = solve_objective_chunked(
        cost_array,
        scenario=scenario,
        objective=objective,
        num_slices=slice_count,
        num_turbines=num_turbines,
        space_px=spacing,
        lcoe_thresh=cost_threshold,
        min_turbines=min_turbines,
        max_turbines=max_turbines,
        axis=0,
    )
    locations.to_csv(locations_path)

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
        scenario=scenario,
        objective=objective,
        num_slices=slice_count,
        num_turbines=num_turbines,
        space_px=spacing,
        lcoe_thresh=cost_threshold,
        min_turbines=min_turbines,
        max_turbines=max_turbines,
        axis=0,
    )
    metadata["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
    metadata["locations_file"] = locations_path.name
    metadata["cost_file"] = cost_path.name
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    logger.info("Wrote stage-4 costs for %s to %s", objective, cost_path)
