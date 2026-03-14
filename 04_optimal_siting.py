# %% imports
from config import atlas_path, study
from scow.allocation import load_allocation_inputs, optimize_allocation
from scow.logging_utils import configure_script_logging
from scow.social_cost import open_social_cost_atlas
from social_cost_inference.optimization import OptimizationRequest


def build_optimization_request(objective: str) -> OptimizationRequest:
    """Map current study defaults onto the reusable optimization request vocabulary."""

    return OptimizationRequest(
        scenario=study.scenario,
        objective=objective,
        spacing=study.spacing,
        num_turbines=study.num_turbines,
        slice_count=study.slice_count,
        min_turbines=study.min_turbines,
        max_turbines=study.max_turbines,
        cost_threshold=study.objective_threshold(objective),
    )


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
            request = build_optimization_request(objective)
            optimize_allocation(
                request=request,
                power=power,
                energy=energy,
                social_cost=social_cost,
                logger=logger,
            )
    finally:
        social_cost.close()


if __name__ == "__main__":
    main()
