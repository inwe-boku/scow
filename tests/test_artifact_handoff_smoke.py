import json
import subprocess

import pandas as pd
import xarray as xr
from shapely.geometry import Point

from conftest import R_SCRIPT, ensure_stage2_runtime
from scow import paths
from scow.artifact_validators import validate_soco_dataset
from scow.postprocess_sites_helpers import build_cost_comparison
from scow.stage4_outputs import load_stage4_bundle

def test_stage2_artifact_handoff_smoke(tmp_path):
    ensure_stage2_runtime()
    (tmp_path / "data" / "processed").mkdir(parents=True)
    design = pd.DataFrame(
        {
            "min_lcoe": [1.0, 1.1, 1.2, 2.0, 2.1, 2.2],
            "zoning": [1, 1, 1, 0, 0, 0],
        }
    )
    design.to_csv(tmp_path / "data" / "processed" / "dc_data_test.csv", index=False)

    subprocess.run(
        [
            "Rscript",
            "--vanilla",
            str(R_SCRIPT),
            str(tmp_path),
            "test",
            "2",
            "smoke",
            "0",
            "42",
            "min_lcoe",
        ],
        cwd=R_SCRIPT.parent,
        check=True,
        capture_output=True,
        text=True,
    )

    assert (tmp_path / "data" / "results" / "spatialdc_coefs_test_smoke.csv").exists()
    assert (tmp_path / "data" / "results" / "spatialdc_loglik_test_smoke.csv").exists()


def test_stage4_to_stage5_artifact_handoff_smoke(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "repo", tmp_path)

    soco = xr.Dataset(
        {
            "min_lcoe": xr.DataArray([[50.0, 55.0]]),
            "Local Social Cost": xr.DataArray([[10.0, 12.0]], attrs={"unit": "€/MWh"}),
        }
    )
    soco.to_netcdf(paths.soco_path("2014"))
    validate_soco_dataset(xr.open_dataset(paths.soco_path("2014")))

    for objective, values in {
        "lcoe": [50.0, 55.0],
        "soco": [60.0, 67.0],
    }.items():
        pd.DataFrame({"l": [0, 0], "b": [0, 1]}).to_csv(paths.opt_locations_path("2014", objective))
        pd.DataFrame(
            {
                objective: values,
                "Energy": [4000.0, 5000.0],
                "CumSumEnergy": [4.0, 9.0],
                "geometry": [Point(15.0, 47.0).wkt, Point(15.2, 47.1).wkt],
            }
        ).to_csv(paths.opt_cost_path("2014", objective), index=False)
        paths.opt_metadata_path("2014", objective).write_text(
            json.dumps(
                {
                    "scenario": "2014",
                    "objective": objective,
                    "spacing": 3,
                    "num_slices": 6,
                    "num_turbines_mode": "auto",
                    "min_turbines": 1200,
                    "max_turbines": 6000,
                    "cost_threshold_eur_per_mwh": 75.0 if objective == "lcoe" else 125.0,
                    "axis": 0,
                    "final_fingerprint": f"{objective}-fingerprint",
                    "generated_at_utc": "2026-03-14T00:00:00+00:00",
                    "locations_file": f"opt_locations_2014_{objective}.csv",
                    "cost_file": f"opt_cost_2014_{objective}.csv",
                }
            )
        )

    lcoe_bundle = load_stage4_bundle("2014", "lcoe")
    soco_bundle = load_stage4_bundle("2014", "soco")
    lcoe_costs = lcoe_bundle["costs"].copy()
    soco_costs = soco_bundle["costs"].copy()
    lcoe_costs["loco"] = [10.0, 12.0]
    lcoe_costs["soco"] = lcoe_costs["lcoe"] + lcoe_costs["loco"]
    soco_costs["lcoe"] = [50.0, 55.0]
    soco_costs["loco"] = soco_costs["soco"] - soco_costs["lcoe"]
    for df in (lcoe_costs, soco_costs):
        df["Annual Private Cost"] = df["lcoe"] * df["Energy"] / 1000
        df["Annual Local Social Cost"] = df["loco"] * df["Energy"] / 1000
        df["Annual Social Cost"] = df["soco"] * df["Energy"] / 1000

    comparison = build_cost_comparison(lcoe_costs, soco_costs, 8.0)

    assert comparison.loc["Private Cost", "Wind Turbines"] == 2
    assert comparison.loc["Social Cost", "Total Local Social Cost"] > 0
