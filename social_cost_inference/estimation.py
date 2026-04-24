from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

from social_cost_inference.studies import SpatialDecisionStudy


@dataclass(frozen=True)
class EstimationRequest:
    """Minimal reusable description of one estimation run."""

    scenario: str
    artifact_key: str
    method: str
    search_mode: str
    num_runs: int
    selection_criterion: str
    seed: int
    output_label: str | None = None


@dataclass(frozen=True)
class EstimationArtifacts:
    """Minimal reusable metadata about one estimation result set."""

    scenario: str
    artifact_key: str
    selected_run_name: str
    output_label: str | None = None


class DecisionEstimator(Protocol):
    """Minimal estimator seam for study-specific implementations."""

    name: str

    def estimate_preferences(
        self,
        study: SpatialDecisionStudy,
        request: EstimationRequest,
        *,
        compute: bool = True,
        report: bool = True,
    ) -> EstimationArtifacts:
        """Estimate preferences for one study/request pair."""
