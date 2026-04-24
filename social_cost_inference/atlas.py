from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class SocialCostAtlas:
    """Path-based handle for the current social-cost artifact seam."""

    scenario: str
    artifact_key: str
    netcdf_path: Path


SocialCostArtifacts = SocialCostAtlas
