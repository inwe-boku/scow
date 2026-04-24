from __future__ import annotations

import argparse

from scow.estimation import estimate_preferences_for_study, presets
from scow.logging_utils import configure_script_logging
from social_cost_inference.estimation import EstimationRequest
from social_cost_inference.studies import SpatialDecisionStudy


def build_default_spec(study: SpatialDecisionStudy):
    """Build the canonical stage-2 spec from the selected study."""
    return study.estimation_spec()


def build_request_from_spec(spec) -> EstimationRequest:
    """Map the current stage-2 spec onto the reusable request vocabulary."""

    return EstimationRequest(
        scenario=spec.scenario,
        artifact_key=spec.artifact_key,
        method=spec.method,
        search_mode=spec.search_mode,
        num_runs=spec.num_runs,
        selection_criterion=spec.selection_criterion,
        seed=spec.seed,
        output_label=spec.output_label,
    )


def parse_args(argv=None):
    """Parse the thin stage-2 CLI surface."""
    parser = argparse.ArgumentParser(description="Run stage-2 estimation.")
    parser.add_argument(
        "--preset",
        choices=("baseline", "single_check", "fast_debug", "modern_baseline"),
        help="Named starting configuration. Defaults to the current study-driven baseline.",
    )
    parser.add_argument(
        "--scenario",
        choices=("2014", "modern"),
        help="Legacy runtime scenario label. Use the explicit study registry for project-level study IDs.",
    )
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


def build_spec_from_args(study: SpatialDecisionStudy, args):
    """Construct an estimation spec from CLI arguments."""
    study_id = getattr(study, "resolved_study_id", None)
    if args.preset is None:
        spec = build_default_spec(study)
    else:
        spec = getattr(presets, args.preset)()

    if study_id is not None:
        if args.scenario is not None and args.scenario != study.scenario:
            raise ValueError(
                f"Explicit study {study_id!r} uses runtime scenario {study.scenario!r}; "
                f"got conflicting override {args.scenario!r}"
            )
        spec = spec.for_scenario(study.scenario).with_artifact_key(study.artifact_key)
    elif args.scenario is not None:
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


def run_estimation_cli(study: SpatialDecisionStudy, argv=None):
    """Run the canonical stage-2 workflow from the current script CLI."""
    logger = configure_script_logging("stage2")
    args = parse_args(argv)
    spec = build_spec_from_args(study, args)
    request = build_request_from_spec(spec)
    logger.info(
        "Running stage 2 for scenario=%s artifact_key=%s method=%s num_runs=%s",
        spec.scenario,
        spec.artifact_key,
        spec.method,
        spec.num_runs,
    )
    return estimate_preferences_for_study(study, request, compute=True, report=True)
