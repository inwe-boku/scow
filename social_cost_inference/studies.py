from __future__ import annotations

from typing import Any, Protocol, runtime_checkable

from social_cost_inference.optimization import OptimizationRequest
from social_cost_inference.reporting import ReportingRequest


@runtime_checkable
class SpatialDecisionStudy(Protocol):
    """Minimal study seam for the future reusable core."""

    scenario: str

    def estimation_spec(self) -> Any:
        """Return the project-local estimation specification for the current study."""

    def optimization_request(self, objective: str) -> OptimizationRequest:
        """Return the project-local optimization request for one named objective."""

    def reporting_request(self) -> ReportingRequest:
        """Return the project-local reporting request for the current study."""
