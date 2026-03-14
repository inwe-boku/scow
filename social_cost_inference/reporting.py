from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ReportingRequest:
    """Minimal reusable description of one reporting run."""

    scenario: str
    area_name: str
    generation_target: float


@dataclass(frozen=True)
class ReportArtifacts:
    """Minimal path-based handle for reporting outputs."""

    figure_paths: tuple[Path, ...]
    table_paths: tuple[Path, ...]
