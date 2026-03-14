from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class SocialCostArtifacts:
    """Path-based handle for the current social-cost artifact seam."""

    scenario: str
    netcdf_path: Path


SocialCostAtlas = SocialCostArtifacts
