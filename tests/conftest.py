from pathlib import Path
import shutil
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


R_SCRIPT = ROOT / "discrete_choice.R"


def ensure_stage2_runtime():
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript is not available")

    command = [
        "Rscript",
        "--vanilla",
        "-e",
        (
            'if (!requireNamespace("renv", quietly = TRUE)) quit(status = 1); '
            'renv::load(project = ".", quiet = TRUE); '
            'for (pkg in c("logitr", "tidyverse", "broom")) { '
            'if (!requireNamespace(pkg, quietly = TRUE)) quit(status = 1) }'
        ),
    ]
    result = subprocess.run(command, cwd=R_SCRIPT.parent, capture_output=True, text=True)
    if result.returncode != 0:
        pytest.skip("Stage-2 R runtime is not provisioned")
