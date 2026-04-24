from config import study
from scow.estimation_cli import run_estimation_cli


def main(argv=None):
    """Run the canonical stage-2 estimation workflow."""
    run_estimation_cli(study, argv=argv)


if __name__ == "__main__":
    main()
