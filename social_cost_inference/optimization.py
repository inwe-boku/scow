from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class OptimizationRequest:
    """Minimal reusable description of one optimization solve."""

    scenario: str
    objective: str
    spacing: int
    num_turbines: int | str
    slice_count: int
    min_turbines: int
    max_turbines: int
    cost_threshold: float


@dataclass(frozen=True)
class OptimizationArtifacts:
    """Minimal path-based handle for one optimization output family."""

    scenario: str
    objective: str
    locations_path: Path
    cost_path: Path
    metadata_path: Path
