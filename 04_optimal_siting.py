# %% imports
from config import atlas_path, study
from scow.allocation import load_allocation_inputs, optimize_allocation
from scow.logging_utils import configure_script_logging


def main() -> None:
    """Run the canonical stage-4 batch over the current objective family."""

    logger = configure_script_logging("stage4")
    logger.info("Loading atlas and social-cost inputs for stage 4")
    power, energy, social_cost = load_allocation_inputs(
        atlas_file=atlas_path,
        scenario=study.scenario,
    )

    for objective in study.objective_family:
        optimize_allocation(
            objective=objective,
            scenario=study.scenario,
            spacing=study.spacing,
            num_turbines=study.num_turbines,
            power=power,
            energy=energy,
            social_cost=social_cost,
            cost_threshold=study.objective_threshold(objective),
            slice_count=study.slice_count,
            min_turbines=study.min_turbines,
            max_turbines=study.max_turbines,
            logger=logger,
        )


if __name__ == "__main__":
    main()
