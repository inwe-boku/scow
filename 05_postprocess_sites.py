from config import repo, study
from scow.reporting import generate_reports


def main() -> None:
    """Run the canonical stage-5 reporting step."""

    generate_reports(
        repo_root=repo,
        scenario=study.scenario,
        area_name=study.area_name,
        generation_target=study.generation_target,
    )


if __name__ == "__main__":
    main()
