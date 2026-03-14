from config import atlas_path, repo, study
from scow.social_cost import build_social_cost_atlas


def main() -> None:
    """Run the canonical stage-3 social-cost construction step."""

    build_social_cost_atlas(
        repo_root=repo,
        atlas_file=atlas_path,
        scenario=study.scenario,
    )


if __name__ == "__main__":
    main()
