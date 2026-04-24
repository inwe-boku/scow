"""Load final stage-4 outputs through one canonical path."""

from __future__ import annotations

import json

import pandas as pd

from scow.artifact_validators import validate_stage4_bundle
from scow.paths import opt_cost_path, opt_locations_path, opt_metadata_path
from scow.postprocess_sites_helpers import to_geopandas


def _require_file(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing stage-4 artifact: {path}")
    return path


def load_stage4_locations(artifact_key: str, objective: str) -> pd.DataFrame:
    """Load the final stage-4 locations CSV as the plain `l`/`b` table."""
    path = _require_file(opt_locations_path(artifact_key, objective))
    locations = pd.read_csv(path)
    if "Unnamed: 0" in locations.columns:
        locations = locations.drop(columns=["Unnamed: 0"])
    if {"l", "b"}.issubset(locations.columns):
        locations[["l", "b"]] = locations[["l", "b"]].astype(int)
    return locations


def load_stage4_costs(artifact_key: str, objective: str):
    """Load the final stage-4 cost CSV as a GeoDataFrame."""
    path = _require_file(opt_cost_path(artifact_key, objective))
    return to_geopandas(pd.read_csv(path))


def load_stage4_metadata(artifact_key: str, objective: str) -> dict:
    """Load the final stage-4 metadata JSON."""
    path = _require_file(opt_metadata_path(artifact_key, objective))
    return json.loads(path.read_text())


def load_stage4_bundle(artifact_key: str, objective: str, validate: bool = True) -> dict:
    """Load final stage-4 locations, costs, and metadata together."""
    locations = load_stage4_locations(artifact_key, objective)
    costs = load_stage4_costs(artifact_key, objective)
    metadata = load_stage4_metadata(artifact_key, objective)

    if validate:
        validate_stage4_bundle(locations, costs, metadata, artifact_key=artifact_key, objective=objective)

    return {"locations": locations, "costs": costs, "metadata": metadata}
