from config import data_ver, repo
from scow.reporting import generate_reports

AREA_NAME = "Niederösterreich"
GENERATION_TARGET = 8.0


def main() -> None:
    """Run the canonical stage-5 reporting step."""

    generate_reports(
        repo_root=repo,
        scenario=data_ver,
        area_name=AREA_NAME,
        generation_target=GENERATION_TARGET,
    )


if __name__ == "__main__":
    main()
