"""Minimal shared logging setup for numbered pipeline scripts."""

from __future__ import annotations

import logging


class _StageFilter(logging.Filter):
    def __init__(self, stage_name: str):
        super().__init__()
        self.stage_name = stage_name

    def filter(self, record: logging.LogRecord) -> bool:
        # The adapter already sets `stage` for direct script logs; this backfills
        # the field for records emitted outside the adapter so the formatter stays stable.
        record.stage = self.stage_name
        return True


def configure_script_logging(stage_name: str, level: int = logging.INFO) -> logging.LoggerAdapter:
    """Configure consistent stderr logging for one numbered pipeline script."""
    root = logging.getLogger()
    root.handlers.clear()
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(asctime)s [%(stage)s] %(levelname)s %(message)s"))
    handler.addFilter(_StageFilter(stage_name))
    root.addHandler(handler)
    root.setLevel(level)
    return logging.LoggerAdapter(logging.getLogger(stage_name), {"stage": stage_name})
