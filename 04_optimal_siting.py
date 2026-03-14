# %% imports
from config import atlas_path, data_ver, nturbines, spacing
from scow.allocation import (
    load_allocation_inputs,
    optimize_allocation,
)
from scow.logging_utils import configure_script_logging

SLICE_COUNT = 6
MIN_TURBINES = 1200
MAX_TURBINES = 6000
OBJECTIVE_COST_THRESHOLDS_EUR_PER_MWH = {
    "loco": 100,
    "lcoe": 75,
    "soco": 125,
}


def main() -> None:
    """Run the canonical stage-4 batch over the current objective family."""

    logger = configure_script_logging("stage4")
    logger.info("Loading atlas and social-cost inputs for stage 4")
    power, energy, social_cost = load_allocation_inputs(
        atlas_file=atlas_path,
        scenario=data_ver,
    )

    for objective, cost_threshold in OBJECTIVE_COST_THRESHOLDS_EUR_PER_MWH.items():
        optimize_allocation(
            objective=objective,
            scenario=data_ver,
            spacing=spacing,
            num_turbines=nturbines,
            power=power,
            energy=energy,
            social_cost=social_cost,
            cost_threshold=cost_threshold,
            slice_count=SLICE_COUNT,
            min_turbines=MIN_TURBINES,
            max_turbines=MAX_TURBINES,
            logger=logger,
        )


if __name__ == "__main__":
    main()
