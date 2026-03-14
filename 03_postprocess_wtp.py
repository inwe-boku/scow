from config import atlas_path, data_ver, repo
from scow.social_cost import build_social_cost_atlas


def main() -> None:
    """Run the canonical stage-3 social-cost construction step."""

    build_social_cost_atlas(
        repo_root=repo,
        atlas_file=atlas_path,
        scenario=data_ver,
    )


if __name__ == "__main__":
    main()
