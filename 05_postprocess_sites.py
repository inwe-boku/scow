from config import repo, study
from scow.reporting import generate_reports


def main() -> None:
    """Run the canonical stage-5 reporting step."""

    generate_reports(
        repo_root=repo,
        request=study.reporting_request(),
    )


if __name__ == "__main__":
    main()
