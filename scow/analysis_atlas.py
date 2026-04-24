from __future__ import annotations

from pathlib import Path

import cleo

from scow.study import ScowStudy


def load_analysis_atlas(
    study: ScowStudy,
    atlas_root: Path,
    *,
    require_optimal_power: bool = False,
    require_optimal_energy: bool = False,
) -> cleo.Atlas:
    """Load existing atlas stores for downstream analysis stages.

    This is a load path for stages 3 to 5. It is not a second build path.
    """
    wind_store = atlas_root / "wind.zarr"
    landscape_store = atlas_root / "landscape.zarr"
    if not wind_store.exists():
        raise FileNotFoundError(f"wind.zarr not found at {wind_store}. Run stage 1 first.")
    if not landscape_store.exists():
        raise FileNotFoundError(f"landscape.zarr not found at {landscape_store}. Run stage 1 first.")

    atlas = cleo.Atlas(atlas_root, "AUT", "EPSG:31287")
    atlas.select(area=study.area_name, inplace=True)
    atlas.wind.select(turbines=study.turbines)

    if "capacity_factors" not in atlas.wind.data.data_vars:
        atlas.wind.compute("capacity_factors", loss_factor=study.loss_factor).materialize(
            overwrite=True,
            allow_method_change=True,
        )
    if "lcoe" not in atlas.wind.data.data_vars:
        atlas.wind.compute(
            "lcoe",
            cf={"loss_factor": study.loss_factor},
            economics=study.economics,
        ).materialize(overwrite=True)
    if "min_lcoe" not in atlas.wind.data.data_vars:
        atlas.wind.data["min_lcoe"] = atlas.wind.data["lcoe"].min(dim="turbine").rename("min_lcoe")
    if require_optimal_power and "optimal_power" not in atlas.wind.data.data_vars:
        atlas.wind.compute(
            "optimal_power",
            cf={"loss_factor": study.loss_factor},
            economics=study.economics,
        ).materialize(overwrite=True)
    if require_optimal_energy and "optimal_energy" not in atlas.wind.data.data_vars:
        atlas.wind.compute(
            "optimal_energy",
            cf={"loss_factor": study.loss_factor},
            economics=study.economics,
        ).materialize(overwrite=True)

    return atlas
