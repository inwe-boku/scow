from pathlib import Path
import shutil
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


R_SCRIPT = ROOT / "discrete_choice.R"
STAGE2_RUNTIME_PROBE_TIMEOUT_SECONDS = 15


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
            'for (pkg in c("logitr", "tidyverse", "broom", "jsonlite")) { '
            'if (!requireNamespace(pkg, quietly = TRUE)) quit(status = 1) }; '
            'rscript_path <- Sys.which("Rscript"); '
            'if (rscript_path == "") quit(status = 1); '
            'cl <- tryCatch('
            'parallel::makePSOCKcluster('
            '1, master = "127.0.0.1", setup_strategy = "sequential", outfile = "", rscript = rscript_path'
            '), '
            'error = function(e) NULL'
            '); '
            'if (is.null(cl)) quit(status = 1); '
            'parallel::stopCluster(cl)'
        ),
    ]
    try:
        result = subprocess.run(
            command,
            cwd=R_SCRIPT.parent,
            capture_output=True,
            text=True,
            timeout=STAGE2_RUNTIME_PROBE_TIMEOUT_SECONDS,
        )
    except subprocess.TimeoutExpired:
        pytest.skip("Stage-2 R runtime probe timed out")
    if result.returncode != 0:
        pytest.skip("Stage-2 R runtime is not provisioned")
