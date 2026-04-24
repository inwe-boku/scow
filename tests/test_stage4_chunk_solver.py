import numpy as np
import pandas as pd
import pytest
import xarray as xr

from scow.optimization import (
    candidate_frame,
    conflict_edges,
    final_output_metadata,
    reconstruct_chunk,
    solve_chunk_gurobi,
)


def test_candidate_frame_keeps_finite_zero_costs_and_excludes_nan():
    chunk_map = pd.DataFrame(
        [[0.0, np.nan, 2.5], [1.0, 0.5, np.nan]],
        index=[0, 1],
        columns=[0, 1, 2],
    )

    candidates = candidate_frame(chunk_map)

    assert candidates.to_dict(orient="records") == [
        {"l": 0, "b": 0, "cost": 0.0},
        {"l": 0, "b": 2, "cost": 2.5},
        {"l": 1, "b": 0, "cost": 1.0},
        {"l": 1, "b": 1, "cost": 0.5},
    ]


def test_conflict_edges_match_current_spacing_stencil_for_space_px_2():
    candidates = pd.DataFrame(
        [
            {"l": 0, "b": 0, "cost": 1.0},
            {"l": 0, "b": 1, "cost": 2.0},
            {"l": 1, "b": 0, "cost": 3.0},
            {"l": 1, "b": 1, "cost": 4.0},
        ]
    )

    conflicts = conflict_edges(candidates, space_px=2)

    assert conflicts == [(0, 1), (0, 2), (1, 3), (2, 3)]


def test_reconstruct_chunk_preserves_current_chunk_boundaries():
    cost_array = xr.DataArray(
        np.array(
            [
                [1.0, 2.0, np.nan, 5.0],
                [1.5, 2.5, 3.5, np.nan],
            ]
        ),
        dims=("y", "x"),
        coords={"y": [0, 1], "x": [0, 1, 2, 3]},
    )

    first = reconstruct_chunk(
        cost_array,
        chunk_index=0,
        num_slices=2,
        num_turbines="auto",
        space_px=2,
        lcoe_thresh=10.0,
        min_turbines=1,
        max_turbines=10,
    )
    second = reconstruct_chunk(
        cost_array,
        chunk_index=1,
        num_slices=2,
        num_turbines="auto",
        space_px=2,
        lcoe_thresh=10.0,
        min_turbines=1,
        max_turbines=10,
    )

    assert first["x_start"] == 0
    assert first["x_stop"] == 2
    assert second["x_start"] == 2
    assert second["x_stop"] == 4


def test_solve_chunk_gurobi_prefers_low_cost_non_conflicting_sites():
    gp = pytest.importorskip("gurobipy")
    candidates = pd.DataFrame(
        [
            {"l": 0, "b": 0, "cost": 1.0},
            {"l": 0, "b": 1, "cost": 9.0},
            {"l": 1, "b": 0, "cost": 2.0},
            {"l": 1, "b": 1, "cost": 3.0},
        ]
    )
    conflicts = conflict_edges(candidates, space_px=2)

    try:
        selected, info = solve_chunk_gurobi(
            candidates,
            conflicts,
            num_turbines=2,
            threads=1,
        )
    except gp.GurobiError as exc:
        pytest.skip(f"Gurobi unavailable for solver smoke test: {exc}")

    assert selected.to_dict(orient="records") == [{"l": 0, "b": 0}, {"l": 1, "b": 1}]
    assert info["solution_count"] >= 1


def test_final_output_metadata_fingerprint_is_stable_and_sensitive_to_settings():
    cost_array = xr.DataArray(
        np.array(
            [
                [1.0, 2.0, np.nan, 5.0],
                [1.5, 2.5, 3.5, np.nan],
            ]
        ),
        dims=("y", "x"),
        coords={"y": [0, 1], "x": [0, 1, 2, 3]},
    )

    first = final_output_metadata(
        cost_array,
        artifact_key="2014",
        objective="lcoe",
        num_slices=2,
        num_turbines="auto",
        space_px=2,
        lcoe_thresh=75.0,
        min_turbines=1,
        max_turbines=10,
    )
    second = final_output_metadata(
        cost_array,
        artifact_key="2014",
        objective="lcoe",
        num_slices=2,
        num_turbines="auto",
        space_px=2,
        lcoe_thresh=75.0,
        min_turbines=1,
        max_turbines=10,
    )
    changed = final_output_metadata(
        cost_array,
        artifact_key="2014",
        objective="lcoe",
        num_slices=2,
        num_turbines="auto",
        space_px=3,
        lcoe_thresh=75.0,
        min_turbines=1,
        max_turbines=10,
    )

    assert first["final_fingerprint"] == second["final_fingerprint"]
    assert first["final_fingerprint"] != changed["final_fingerprint"]
    assert first["min_turbines"] == 1
    assert first["max_turbines"] == 10
    assert first["cost_threshold_eur_per_mwh"] == 75.0
