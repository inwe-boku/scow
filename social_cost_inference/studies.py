from __future__ import annotations

from typing import Any, Protocol, runtime_checkable


@runtime_checkable
class SpatialDecisionStudy(Protocol):
    """Minimal study seam for the future reusable core."""

    scenario: str

    def estimation_spec(self) -> Any:
        """Return the project-local estimation specification for the current study."""
