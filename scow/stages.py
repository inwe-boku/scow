from pathlib import Path
from typing import Mapping

import cleo

AREA_NAME = "Niederösterreich"

TURBINES_2014 = [
    "Enercon.E40.500",
    "Enercon.E82.3000",
    "Enercon.E101.3050",
    "Enercon.E115.3000",
    "Vestas.V100.1800",
    "Vestas.V100.2000",
    "Vestas.V112.3075",
]
TURBINES_MODERN = [
    "Enercon.E138.3500",
    "Enercon.E160.5560",
    "Vestas.V150.4200",
]

LOSS_FACTOR = 0.71197
ECONOMICS = {
    "discount_rate": 0.05,
    "lifetime_a": 25,
    "om_fixed_eur_per_kw_a": 20.0,
    "om_variable_eur_per_kwh": 0.008,
    "bos_cost_share": 0.3,
    "grid_connect_cost_eur_per_kw": 0.0,
}


def _turbines_for_scenario(scenario: str) -> list[str]:
    if scenario == "2014":
        return TURBINES_2014
    if scenario == "modern":
        return TURBINES_MODERN
    raise ValueError(f"Unknown scenario: {scenario}")


def load_atlas_for_analysis(
    config: Mapping[str, object],
    *,
    require_optimal_power: bool = False,
    require_optimal_energy: bool = False,
) -> cleo.Atlas:
    """Create an analysis atlas from existing CLEO stores.

    This helper centralizes the current baseline-safe CLEO loading path used by
    stages 3 to 5. It intentionally avoids fragile calls such as
    `min_lcoe_turbine`, which are known to break on the current stack.
    """
    atlas_root = Path(config["atlas_path"])
    scenario = str(config["scenario"])

    wind_store = atlas_root / "wind.zarr"
    landscape_store = atlas_root / "landscape.zarr"
    if not wind_store.exists():
        raise FileNotFoundError(
            f"wind.zarr not found at {wind_store}. Run stage 1 first."
        )
    if not landscape_store.exists():
        raise FileNotFoundError(
            f"landscape.zarr not found at {landscape_store}. Run stage 1 first."
        )

    atlas = cleo.Atlas(atlas_root, "AUT", "EPSG:31287")
    atlas.select(area=AREA_NAME, inplace=True)
    atlas.wind.select(turbines=_turbines_for_scenario(scenario))

    if "capacity_factors" not in atlas.wind.data.data_vars:
        atlas.wind.compute("capacity_factors", loss_factor=LOSS_FACTOR).materialize(
            overwrite=True,
            allow_method_change=True,
        )
    if "lcoe" not in atlas.wind.data.data_vars:
        atlas.wind.compute(
            "lcoe",
            cf={"loss_factor": LOSS_FACTOR},
            economics=ECONOMICS,
        ).materialize(overwrite=True)
    if "min_lcoe" not in atlas.wind.data.data_vars:
        atlas.wind.data["min_lcoe"] = atlas.wind.data["lcoe"].min(dim="turbine").rename("min_lcoe")
    if require_optimal_power and "optimal_power" not in atlas.wind.data.data_vars:
        atlas.wind.compute(
            "optimal_power",
            cf={"loss_factor": LOSS_FACTOR},
            economics=ECONOMICS,
        ).materialize(overwrite=True)
    if require_optimal_energy and "optimal_energy" not in atlas.wind.data.data_vars:
        atlas.wind.compute(
            "optimal_energy",
            cf={"loss_factor": LOSS_FACTOR},
            economics=ECONOMICS,
        ).materialize(overwrite=True)

    return atlas
