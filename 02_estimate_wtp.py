# %% imports
import argparse

from config import study
from scow.estimation import estimate_preferences, presets
from scow.logging_utils import configure_script_logging


def build_default_spec():
    """Build the canonical stage-2 spec from the current config module."""
    return study.estimation_spec()


def parse_args(argv=None):
    """Parse the thin stage-2 CLI surface."""
    parser = argparse.ArgumentParser(description="Run stage-2 estimation.")
    parser.add_argument(
        "--preset",
        choices=("baseline", "single_check", "fast_debug", "modern_baseline"),
        help="Named starting configuration. Defaults to the current config-driven baseline.",
    )
    parser.add_argument("--scenario", choices=("2014", "modern"))
    parser.add_argument(
        "--method",
        help="Econometric method override. Canonical execution currently supports mnl_logitr only.",
    )
    parser.add_argument("--num-runs", type=int)
    parser.add_argument("--selection-criterion", choices=("AIC", "BIC"))
    parser.add_argument("--seed", type=int)

    search_group = parser.add_mutually_exclusive_group()
    search_group.add_argument("--single", action="store_true")
    search_group.add_argument("--full-significant", action="store_true")
    search_group.add_argument("--gts-leave-one-out", action="store_true")

    return parser.parse_args(argv)


def build_spec_from_args(args):
    """Construct an estimation spec from CLI arguments."""
    if args.preset is None:
        spec = build_default_spec()
    else:
        spec = getattr(presets, args.preset)()

    if args.scenario is not None:
        spec = spec.for_scenario(args.scenario)
    if args.method is not None:
        spec = spec.using_method(args.method)
    if args.num_runs is not None:
        spec = spec.with_num_runs(args.num_runs)
    if args.selection_criterion is not None:
        spec = spec.with_selection_criterion(args.selection_criterion)
    if args.seed is not None:
        spec = spec.with_seed(args.seed)

    if args.single:
        spec = spec.single()
    elif args.full_significant:
        spec = spec.full_significant()
    elif args.gts_leave_one_out:
        spec = spec.gts_leave_one_out()

    return spec


def main(argv=None):
    """Run the canonical stage-2 estimation workflow."""
    logger = configure_script_logging("stage2")
    args = parse_args(argv)
    spec = build_spec_from_args(args)
    logger.info("Running stage 2 for scenario=%s method=%s num_runs=%s", spec.scenario, spec.method, spec.num_runs)
    estimate_preferences(spec, compute=True, report=True)


if __name__ == "__main__":
    main()
