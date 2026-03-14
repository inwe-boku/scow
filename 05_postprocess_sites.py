from config import repo, study
from scow.reporting import generate_reports
from social_cost_inference.reporting import ReportingRequest


def build_reporting_request() -> ReportingRequest:
    """Map the current study defaults onto the reusable reporting request vocabulary."""

    return ReportingRequest(
        scenario=study.scenario,
        area_name=study.area_name,
        generation_target=study.generation_target,
    )


def main() -> None:
    """Run the canonical stage-5 reporting step."""

    generate_reports(
        repo_root=repo,
        request=build_reporting_request(),
    )


if __name__ == "__main__":
    main()
