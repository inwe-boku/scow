"""Python-native stage-4 chunked optimization helpers."""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path

import gurobipy as gp
import numpy as np
import pandas as pd
from gurobipy import GRB

from scow.paths import chunk_locations_cache_path, chunk_locations_meta_path
from scow.utils import calculate_num_turbines


def reconstruct_chunk(
    cost_array,
    *,
    chunk_index: int,
    num_slices: int,
    num_turbines,
    space_px: int,
    lcoe_thresh: float,
    min_turbines: int,
    max_turbines: int,
    axis: int = 0,
):
    """Reconstruct one stage-4 chunk using the current production chunk logic."""
    cost_df = pd.DataFrame(data=cost_array.data)
    num_pixels = np.ceil(np.round(cost_array.count().data, 0) / num_slices)
    dim_max = cost_array.sizes["x"] if axis == 0 else cost_array.sizes["y"]

    j = 0
    for current_chunk in range(num_slices):
        size = 0
        i = j
        while size <= num_pixels and i < dim_max:
            if axis == 0:
                cost_slice = cost_df.iloc[:, i]
            elif axis == 1:
                cost_slice = cost_df.iloc[i, :]
            else:
                raise ValueError("axis must be 0 or 1")
            cost_slice = cost_slice.dropna()
            i += 1
            size += len(cost_slice)

        if axis == 0:
            chunk_map = cost_df.iloc[:, j:i]
        else:
            chunk_map = cost_df.iloc[j:i, :]

        resolved_num_turbines = calculate_num_turbines(
            num_turbines,
            chunk_map,
            lcoe_thresh,
            space_px,
            max_turbines,
            min_turbines,
        )

        if current_chunk == chunk_index:
            return {
                "chunk_map": chunk_map,
                "chunk_index": current_chunk,
                "x_start": int(j),
                "x_stop": int(i),
                "axis": int(axis),
                "resolved_num_turbines": int(resolved_num_turbines),
            }
        j = i

    raise ValueError(f"Chunk {chunk_index} not found")


def candidate_frame(chunk_map):
    """Return a candidate dataframe with full-array integer `l`, `b`, and `cost` columns."""
    records = []
    for l in chunk_map.index:
        for b in chunk_map.columns:
            cost = float(chunk_map.at[l, b])
            if not np.isfinite(cost):
                continue
            records.append({"l": int(l), "b": int(b), "cost": cost})
    candidates = pd.DataFrame.from_records(records)
    if candidates.empty:
        return pd.DataFrame(columns=["l", "b", "cost"])
    return candidates


def conflict_edges(candidates, *, space_px: int):
    """Return a deduplicated undirected conflict set for the current spacing stencil."""
    coord_to_node = {
        (int(row.l), int(row.b)): int(row.Index)
        for row in candidates.reset_index().itertuples()
    }

    offsets = set()
    for di in range(space_px):
        for dj in range(space_px):
            if (di == 0 and dj == 0) or (di == space_px - 1 and dj == space_px - 1):
                continue
            for sign_l, sign_b in ((1, 1), (1, -1), (-1, -1), (-1, 1)):
                offsets.add((sign_l * di, sign_b * dj))

    conflicts = set()
    for row in candidates.itertuples():
        node = int(row.Index)
        coord = (int(row.l), int(row.b))
        for dl, db in offsets:
            neighbor = (coord[0] + dl, coord[1] + db)
            neighbor_node = coord_to_node.get(neighbor)
            if neighbor_node is None or neighbor_node == node:
                continue
            left, right = sorted((node, neighbor_node))
            conflicts.add((left, right))

    return sorted(conflicts)


def solve_chunk_gurobi(
    candidates,
    conflicts,
    *,
    num_turbines: int,
    threads: int = 4,
    mip_gap: float = 0.01,
    time_limit_s: int | None = None,
):
    """Solve one chunk directly with gurobipy using the current baseline semantics."""
    model = gp.Model("scow_stage4_chunk")
    model.Params.Threads = threads
    model.Params.MIPGap = mip_gap
    model.Params.OutputFlag = 1
    if time_limit_s is not None:
        model.Params.TimeLimit = time_limit_s

    build = model.addVars(
        range(len(candidates)),
        vtype=GRB.BINARY,
        obj=candidates["cost"].to_numpy(),
        name="build",
    )
    model.addConstr(build.sum() >= num_turbines, name="turbine_count")
    model.addConstrs(
        (build[left] + build[right] <= 1 for left, right in conflicts),
        name="spacing",
    )

    model.optimize()

    if model.SolCount == 0:
        raise RuntimeError(f"No feasible chunk solution found (status={model.Status})")

    selected_nodes = [idx for idx in range(len(candidates)) if build[idx].X > 0.5]
    selected = candidates.iloc[selected_nodes][["l", "b"]].copy()
    selected[["l", "b"]] = selected[["l", "b"]].astype(int)
    return selected.reset_index(drop=True), {
        "status": int(model.Status),
        "solution_count": int(model.SolCount),
        "objective_value": float(model.ObjVal),
        "objective_bound": float(model.ObjBound),
        "mip_gap": float(model.MIPGap),
    }


def _chunk_cost_fingerprint(chunk_map) -> str:
    chunk_values = np.asarray(chunk_map.to_numpy(dtype=np.float64, copy=False))
    return hashlib.sha256(chunk_values.tobytes()).hexdigest()


def _hash_chunk_metadata_list(metadata_list) -> str:
    payload = json.dumps(metadata_list, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _chunk_cache_metadata(
    *,
    scenario: str,
    objective: str,
    chunk_index: int,
    chunk_info,
    space_px: int,
    num_slices: int,
    num_turbines,
    lcoe_thresh: float,
    mip_gap: float,
    time_limit_s: int | None,
):
    return {
        "scenario": scenario,
        "objective": objective,
        "chunk_index": int(chunk_index),
        "axis": int(chunk_info["axis"]),
        "x_start": int(chunk_info["x_start"]),
        "x_stop": int(chunk_info["x_stop"]),
        "spacing": int(space_px),
        "num_slices": int(num_slices),
        "num_turbines_mode": str(num_turbines),
        "resolved_num_turbines": int(chunk_info["resolved_num_turbines"]),
        "lcoe_thresh": float(lcoe_thresh),
        "mip_gap": float(mip_gap),
        "time_limit_s": None if time_limit_s is None else int(time_limit_s),
        "cost_fingerprint": _chunk_cost_fingerprint(chunk_info["chunk_map"]),
    }


def _load_cached_chunk(cache_path: Path, meta_path: Path, expected_metadata):
    if not cache_path.exists() or not meta_path.exists():
        return None

    try:
        cached_metadata = json.loads(meta_path.read_text())
        if cached_metadata != expected_metadata:
            return None
        cached = pd.read_parquet(cache_path)
    except Exception as exc:
        logging.warning("Could not read cached chunk %s: %s", cache_path, exc)
        return None

    required_columns = {"l", "b"}
    if set(cached.columns) != required_columns:
        return None

    cached[["l", "b"]] = cached[["l", "b"]].astype(int)
    return cached.sort_values(["l", "b"]).reset_index(drop=True)


def _write_chunk_cache(cache_path: Path, meta_path: Path, locations, metadata):
    locations[["l", "b"]].astype(int).to_parquet(cache_path, index=False)
    meta_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")


def final_output_metadata(
    cost_array,
    *,
    scenario: str,
    objective: str,
    num_slices: int,
    num_turbines,
    space_px: int,
    lcoe_thresh: float,
    min_turbines: int,
    max_turbines: int,
    axis: int = 0,
    mip_gap: float = 0.01,
    time_limit_s: int | None = None,
):
    """Build final stage-4 output metadata from the current chunk-cache contract."""
    chunk_metadata = []
    for chunk_index in range(num_slices):
        chunk_info = reconstruct_chunk(
            cost_array,
            chunk_index=chunk_index,
            num_slices=num_slices,
            num_turbines=num_turbines,
            space_px=space_px,
            lcoe_thresh=lcoe_thresh,
            min_turbines=min_turbines,
            max_turbines=max_turbines,
            axis=axis,
        )
        chunk_metadata.append(
            _chunk_cache_metadata(
                scenario=scenario,
                objective=objective,
                chunk_index=chunk_index,
                chunk_info=chunk_info,
                space_px=space_px,
                num_slices=num_slices,
                num_turbines=num_turbines,
                lcoe_thresh=lcoe_thresh,
                mip_gap=mip_gap,
                time_limit_s=time_limit_s,
            )
        )

    return {
        "scenario": scenario,
        "objective": objective,
        "spacing": int(space_px),
        "num_slices": int(num_slices),
        "num_turbines_mode": str(num_turbines),
        "min_turbines": int(min_turbines),
        "max_turbines": int(max_turbines),
        "cost_threshold_eur_per_mwh": float(lcoe_thresh),
        "axis": int(axis),
        "final_fingerprint": _hash_chunk_metadata_list(chunk_metadata),
    }


def solve_objective_chunked(
    cost_array,
    *,
    scenario: str,
    objective: str,
    num_slices: int,
    num_turbines,
    space_px: int,
    lcoe_thresh: float,
    min_turbines: int,
    max_turbines: int,
    axis: int = 0,
    threads: int = 4,
    mip_gap: float = 0.01,
    time_limit_s: int | None = None,
):
    """Solve one objective across all chunks using Python-native chunk caching."""
    chunk_results = []

    for chunk_index in range(num_slices):
        chunk_info = reconstruct_chunk(
            cost_array,
            chunk_index=chunk_index,
            num_slices=num_slices,
            num_turbines=num_turbines,
            space_px=space_px,
            lcoe_thresh=lcoe_thresh,
            min_turbines=min_turbines,
            max_turbines=max_turbines,
            axis=axis,
        )
        metadata = _chunk_cache_metadata(
            scenario=scenario,
            objective=objective,
            chunk_index=chunk_index,
            chunk_info=chunk_info,
            space_px=space_px,
            num_slices=num_slices,
            num_turbines=num_turbines,
            lcoe_thresh=lcoe_thresh,
            mip_gap=mip_gap,
            time_limit_s=time_limit_s,
        )
        cache_path = chunk_locations_cache_path(scenario, objective, chunk_index)
        meta_path = chunk_locations_meta_path(scenario, objective, chunk_index)
        cached = _load_cached_chunk(cache_path, meta_path, metadata)

        if cached is not None:
            logging.info("Reusing Python-native chunk cache: %s", cache_path.name)
            chunk_results.append(cached)
            continue

        candidates = candidate_frame(chunk_info["chunk_map"])
        if candidates.empty:
            if chunk_info["resolved_num_turbines"] > 0:
                raise RuntimeError(
                    f"Chunk {chunk_index} for {objective} has no feasible candidates, "
                    f"but requires {chunk_info['resolved_num_turbines']} turbines."
                )
            locations = pd.DataFrame(columns=["l", "b"])
            solve_info = {"objective_value": 0.0, "mip_gap": 0.0}
        else:
            conflicts = conflict_edges(candidates, space_px=space_px)
            locations, solve_info = solve_chunk_gurobi(
                candidates,
                conflicts,
                num_turbines=chunk_info["resolved_num_turbines"],
                threads=threads,
                mip_gap=mip_gap,
                time_limit_s=time_limit_s,
            )
        logging.info(
            "Solved chunk %s/%s for %s: selected=%s objective=%.6f gap=%.6f",
            chunk_index + 1,
            num_slices,
            objective,
            len(locations),
            solve_info["objective_value"],
            solve_info["mip_gap"],
        )
        _write_chunk_cache(cache_path, meta_path, locations, metadata)
        chunk_results.append(locations)

    locations = pd.concat(chunk_results, ignore_index=True)
    locations[["l", "b"]] = locations[["l", "b"]].astype(int)
    return locations.sort_values(["l", "b"]).reset_index(drop=True)
