import json

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import Point

from scow import paths
from scow.stage4_outputs import load_stage4_bundle, load_stage4_costs


def _write_stage4_bundle(tmp_path, monkeypatch, *, objective="lcoe", include_objective_column=True):
    monkeypatch.setattr(paths, "repo", tmp_path)
    locations = pd.DataFrame({"l": [1, 2], "b": [3, 4]})
    costs = pd.DataFrame(
        {
            objective if include_objective_column else "wrong_name": [50.0, 55.0],
            "Energy": [2.0, 3.0],
            "CumSumEnergy": [0.002, 0.005],
            "geometry": [Point(15.0, 47.0).wkt, Point(15.2, 47.1).wkt],
        }
    )
    metadata = {
        "scenario": "2014",
        "objective": objective,
        "spacing": 3,
        "num_slices": 6,
        "num_turbines_mode": "auto",
        "min_turbines": 1200,
        "max_turbines": 6000,
        "cost_threshold_eur_per_mwh": 75.0,
        "axis": 0,
        "final_fingerprint": "abc123",
        "generated_at_utc": "2026-03-14T00:00:00+00:00",
        "locations_file": f"opt_locations_2014_{objective}.csv",
        "cost_file": f"opt_cost_2014_{objective}.csv",
    }
    locations.to_csv(paths.opt_locations_path("2014", objective))
    costs.to_csv(paths.opt_cost_path("2014", objective), index=False)
    paths.opt_metadata_path("2014", objective).write_text(json.dumps(metadata))


def test_load_stage4_bundle_returns_valid_bundle(monkeypatch, tmp_path):
    _write_stage4_bundle(tmp_path, monkeypatch, objective="lcoe")

    bundle = load_stage4_bundle("2014", "lcoe")

    assert bundle["locations"]["l"].tolist() == [1, 2]
    assert isinstance(bundle["costs"], gpd.GeoDataFrame)
    assert bundle["metadata"]["objective"] == "lcoe"


def test_load_stage4_costs_returns_geodataframe(monkeypatch, tmp_path):
    _write_stage4_bundle(tmp_path, monkeypatch, objective="soco")

    costs = load_stage4_costs("2014", "soco")

    assert isinstance(costs, gpd.GeoDataFrame)
    assert costs.geometry.x.round(1).tolist() == [15.0, 15.2]


def test_load_stage4_costs_does_not_require_metadata(monkeypatch, tmp_path):
    _write_stage4_bundle(tmp_path, monkeypatch, objective="lcoe")
    paths.opt_metadata_path("2014", "lcoe").unlink()

    costs = load_stage4_costs("2014", "lcoe")

    assert list(costs["lcoe"]) == [50.0, 55.0]


def test_load_stage4_bundle_rejects_missing_objective_column(monkeypatch, tmp_path):
    _write_stage4_bundle(tmp_path, monkeypatch, objective="loco", include_objective_column=False)

    with pytest.raises(ValueError, match="missing columns"):
        load_stage4_bundle("2014", "loco")


def test_load_stage4_bundle_rejects_row_count_mismatch(monkeypatch, tmp_path):
    _write_stage4_bundle(tmp_path, monkeypatch, objective="soco")
    pd.DataFrame({"l": [1], "b": [3]}).to_csv(paths.opt_locations_path("2014", "soco"))

    with pytest.raises(ValueError, match="same number of selected turbines"):
        load_stage4_bundle("2014", "soco")
