# %% imports
from config import atlas_path, study
from scow.allocation import load_allocation_inputs, optimize_allocation
from scow.logging_utils import configure_script_logging
from scow.social_cost import open_social_cost_atlas


def main() -> None:
    """Run the canonical stage-4 batch over the current objective family."""

    logger = configure_script_logging("stage4")
    logger.info("Loading atlas and social-cost inputs for stage 4")
    power, energy, social_cost_atlas = load_allocation_inputs(
        atlas_file=atlas_path,
        scenario=study.scenario,
    )
    social_cost = open_social_cost_atlas(social_cost_atlas)
    try:
        for objective in study.objective_family:
            optimize_allocation(
                request=study.optimization_request(objective),
                power=power,
                energy=energy,
                social_cost=social_cost,
                logger=logger,
            )
    finally:
        social_cost.close()


if __name__ == "__main__":
    main()
