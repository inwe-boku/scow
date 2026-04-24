# %% imports
from config import repo, study
from scow.spatial import prepare_spatial_data


def main() -> None:
    """Run the project-local stage-1 spatial data build and export workflow."""
    prepare_spatial_data(study, repo)


if __name__ == "__main__":
    main()
