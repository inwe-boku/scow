from __future__ import annotations

from pathlib import Path

import xarray as xr

from scow.artifact_validators import validate_soco_dataset
from scow.logging_utils import configure_script_logging
from scow.paths import parsimonious_name_path, soco_path
from scow.postprocess_wtp_helpers import (
    build_interaction_layers,
    compute_social_cost_dataset,
)
from scow.stages import load_atlas_for_analysis
from scow.utils import postprocess_spatialdc


def build_social_cost_atlas(*, repo_root: Path, atlas_file: Path, scenario: str) -> Path:
    """Build the current stage-3 social-cost artifact for one scenario.

    Parameters
    ----------
    repo_root:
        Repository root used for stage-2 artifact lookup.
    atlas_file:
        Path to the CLEO atlas consumed by the analysis stages.
    scenario:
        Scenario label, currently `2014` or `modern`.

    Returns
    -------
    Path
        Path to the written `soco.nc` artifact.
    """

    logger = configure_script_logging("stage3")
    parsimonious_name = parsimonious_name_path(scenario).read_text().strip()

    parsimonious_coeffs, _ = postprocess_spatialdc(repo_root, scenario, parsimonious_name)
    parsimonious_coeffs["valuation"] = (
        parsimonious_coeffs["estimate"] / parsimonious_coeffs.loc["min_lcoe", "estimate"]
    )

    atlas = load_atlas_for_analysis({"atlas_path": atlas_file, "scenario": scenario})
    landscape_data = build_interaction_layers(
        atlas.landscape.data,
        list(parsimonious_coeffs.index),
    )

    if "min_lcoe" not in landscape_data.data_vars:
        landscape_data = landscape_data.merge(atlas.wind.data["min_lcoe"].compute())
    spatial_data = landscape_data[list(parsimonious_coeffs.index)]

    social_cost = compute_social_cost_dataset(spatial_data, parsimonious_coeffs["valuation"])
    output_path = soco_path(scenario)
    social_cost.to_netcdf(output_path)
    logger.info("Wrote social-cost dataset to %s", output_path)
    return output_path


def load_social_cost_atlas(scenario: str, *, engine: str | None = None) -> xr.Dataset:
    """Load and validate the current `soco.nc` artifact for one scenario."""

    open_kwargs = {} if engine is None else {"engine": engine}
    return validate_soco_dataset(xr.open_dataset(soco_path(scenario), **open_kwargs))
