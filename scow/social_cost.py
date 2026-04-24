from __future__ import annotations

from pathlib import Path

import xarray as xr

from scow.artifact_validators import validate_soco_dataset
from scow.analysis_atlas import load_analysis_atlas
from scow.logging_utils import configure_script_logging
from scow.paths import parsimonious_name_path, soco_path
from scow.postprocess_wtp_helpers import (
    build_interaction_layers,
    compute_social_cost_dataset,
)
from scow.study import ScowStudy
from scow.utils import postprocess_spatialdc
from social_cost_inference.atlas import SocialCostAtlas


def build_social_cost_atlas(
    *,
    repo_root: Path,
    atlas_file: Path,
    study: ScowStudy,
) -> SocialCostAtlas:
    """Build the current stage-3 social-cost artifact for one explicit study.

    Parameters
    ----------
    repo_root:
        Repository root used for stage-2 artifact lookup.
    atlas_file:
        Path to the CLEO atlas consumed by the analysis stages.
    study:
        Explicit project-local study selection. Its runtime scenario still controls
        turbine and atlas loading, while its social-cost key controls the stage-2
        lookup namespace and the written stage-3 artifact name.

    Returns
    -------
    SocialCostAtlas
        Path-based handle to the written `soco.nc` artifact.
    """

    logger = configure_script_logging("stage3")
    artifact_key = study.artifact_key
    parsimonious_name = parsimonious_name_path(artifact_key).read_text().strip()

    parsimonious_coeffs, _ = postprocess_spatialdc(repo_root, artifact_key, parsimonious_name)
    parsimonious_coeffs["valuation"] = (
        parsimonious_coeffs["estimate"] / parsimonious_coeffs.loc["min_lcoe", "estimate"]
    )

    atlas = load_analysis_atlas(study, atlas_file)
    landscape_data = build_interaction_layers(
        atlas.landscape.data,
        list(parsimonious_coeffs.index),
    )

    if "min_lcoe" not in landscape_data.data_vars:
        landscape_data = landscape_data.merge(atlas.wind.data["min_lcoe"].compute())
    spatial_data = landscape_data[list(parsimonious_coeffs.index)]

    social_cost = compute_social_cost_dataset(spatial_data, parsimonious_coeffs["valuation"])
    output_path = soco_path(artifact_key)
    social_cost.to_netcdf(output_path)
    logger.info("Wrote social-cost dataset to %s", output_path)
    return SocialCostAtlas(
        scenario=study.scenario,
        artifact_key=artifact_key,
        netcdf_path=output_path,
    )


def get_social_cost_atlas(
    scenario: str | ScowStudy,
    *,
    runtime_scenario: str | None = None,
) -> SocialCostAtlas:
    """Return the current path-based social-cost artifact handle for one study key."""

    if isinstance(scenario, ScowStudy):
        runtime_scenario = scenario.scenario
        artifact_key = scenario.artifact_key
    else:
        runtime_scenario = scenario if runtime_scenario is None else runtime_scenario
        artifact_key = scenario
    return SocialCostAtlas(
        scenario=runtime_scenario,
        artifact_key=artifact_key,
        netcdf_path=soco_path(artifact_key),
    )


def open_social_cost_atlas(
    social_cost_atlas: SocialCostAtlas | Path,
    *,
    engine: str | None = None,
) -> xr.Dataset:
    """Open and validate the current `soco.nc` artifact from an explicit atlas handle."""

    netcdf_path = (
        social_cost_atlas.netcdf_path
        if isinstance(social_cost_atlas, SocialCostAtlas)
        else Path(social_cost_atlas)
    )
    open_kwargs = {} if engine is None else {"engine": engine}
    return validate_soco_dataset(xr.open_dataset(netcdf_path, **open_kwargs))
