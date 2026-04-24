"""Validation helpers for stage artifacts and cross-stage handoff contracts."""

from __future__ import annotations

from typing import Mapping

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr

from scow.paths import opt_cost_path, opt_locations_path


STAGE4_OBJECTIVES = {"lcoe", "loco", "soco"}
STAGE4_METADATA_REQUIRED_KEYS = {
    "artifact_key",
    "objective",
    "spacing",
    "num_slices",
    "num_turbines_mode",
    "min_turbines",
    "max_turbines",
    "cost_threshold_eur_per_mwh",
    "axis",
    "final_fingerprint",
    "generated_at_utc",
    "locations_file",
    "cost_file",
}


def validate_soco_dataset(dataset: xr.Dataset) -> xr.Dataset:
    """Validate the current stage-3 social-cost dataset contract."""
    required_vars = {"min_lcoe", "Local Social Cost"}
    missing = sorted(required_vars - set(dataset.data_vars))
    if missing:
        raise ValueError(f"soco dataset is missing required variables: {', '.join(missing)}")

    for name, da in dataset.data_vars.items():
        unit = da.attrs.get("units") or da.attrs.get("unit")
        if unit != "€/MWh":
            raise ValueError(f"{name} must declare unit '€/MWh'")
        if not str(da.attrs.get("long_name", "")).strip():
            raise ValueError(f"{name} must declare a non-empty long_name")

    return dataset


def validate_stage4_metadata(
    metadata: Mapping[str, object],
    *,
    artifact_key: str | None = None,
    objective: str | None = None,
) -> Mapping[str, object]:
    """Validate the final stage-4 metadata sidecar."""
    missing = sorted(STAGE4_METADATA_REQUIRED_KEYS - set(metadata))
    if missing:
        raise ValueError(f"stage-4 metadata is missing required keys: {', '.join(missing)}")

    meta_objective = str(metadata["objective"])
    if meta_objective not in STAGE4_OBJECTIVES:
        raise ValueError(f"Unknown stage-4 objective in metadata: {meta_objective}")

    if artifact_key is not None and metadata["artifact_key"] != artifact_key:
        raise ValueError(
            f"Metadata artifact_key {metadata['artifact_key']} does not match expected {artifact_key}"
        )
    if objective is not None and meta_objective != objective:
        raise ValueError(f"Metadata objective {meta_objective} does not match expected {objective}")

    expected_locations = opt_locations_path(str(metadata["artifact_key"]), meta_objective).name
    expected_cost = opt_cost_path(str(metadata["artifact_key"]), meta_objective).name
    if metadata["locations_file"] != expected_locations:
        raise ValueError(
            f"Metadata locations_file {metadata['locations_file']} does not match expected {expected_locations}"
        )
    if metadata["cost_file"] != expected_cost:
        raise ValueError(f"Metadata cost_file {metadata['cost_file']} does not match expected {expected_cost}")

    if not str(metadata["final_fingerprint"]).strip():
        raise ValueError("stage-4 metadata final_fingerprint must be non-empty")

    return metadata


def validate_stage4_locations(df: pd.DataFrame) -> pd.DataFrame:
    """Validate the final stage-4 locations table."""
    missing = sorted({"l", "b"} - set(df.columns))
    if missing:
        raise ValueError(f"stage-4 locations are missing required columns: {', '.join(missing)}")

    for col in ("l", "b"):
        numeric = pd.to_numeric(df[col], errors="raise")
        if numeric.isna().any() or not np.equal(numeric.to_numpy(), numeric.astype(int).to_numpy()).all():
            raise ValueError(f"stage-4 locations column {col} must contain integer-valued coordinates")

    return df


def validate_stage4_costs(df: pd.DataFrame | gpd.GeoDataFrame, objective: str) -> pd.DataFrame | gpd.GeoDataFrame:
    """Validate the final stage-4 cost table for one objective."""
    if objective not in STAGE4_OBJECTIVES:
        raise ValueError(f"Unknown stage-4 objective: {objective}")

    required_columns = {"geometry", "Energy", "CumSumEnergy", objective}
    missing = sorted(required_columns - set(df.columns))
    if missing:
        raise ValueError(f"stage-4 costs for {objective} are missing columns: {', '.join(missing)}")

    if isinstance(df, gpd.GeoDataFrame):
        if df.geometry.isna().any():
            raise ValueError(f"stage-4 costs for {objective} contain missing geometries")
    elif df["geometry"].isna().any():
        raise ValueError(f"stage-4 costs for {objective} contain missing geometry values")

    if df[["Energy", "CumSumEnergy", objective]].isna().any().any():
        raise ValueError(f"stage-4 costs for {objective} contain missing Energy, CumSumEnergy, or cost values")

    cumsum = df["CumSumEnergy"].to_numpy()
    if len(cumsum) > 1 and not (cumsum[1:] >= cumsum[:-1]).all():
        raise ValueError(f"stage-4 costs for {objective} must have nondecreasing CumSumEnergy")

    return df


def validate_stage4_bundle(
    locations: pd.DataFrame,
    costs: pd.DataFrame | gpd.GeoDataFrame,
    metadata: Mapping[str, object],
    *,
    artifact_key: str,
    objective: str,
) -> None:
    """Validate one in-memory final stage-4 output bundle."""
    validate_stage4_metadata(metadata, artifact_key=artifact_key, objective=objective)
    validated_locations = validate_stage4_locations(locations)
    validated_costs = validate_stage4_costs(costs, objective)

    if len(validated_locations) != len(validated_costs):
        raise ValueError(
            "stage-4 locations and costs must contain the same number of selected turbines "
            f"(got {len(validated_locations)} locations and {len(validated_costs)} costs)"
        )
