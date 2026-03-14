import subprocess
import shutil
from pathlib import Path

import pandas as pd
import pytest


BASELINE_TERMS = ["min_lcoe"]
R_SCRIPT = Path(__file__).resolve().parents[1] / "discrete_choice.R"


def _ensure_stage2_runtime():
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


def _write_design_matrix(path: Path):
    zone_in = pd.DataFrame(
        {
            "min_lcoe": [1.0, 1.2, 1.4, 1.7, 1.9, 2.1, 2.3, 2.6],
            "zoning": [1] * 8,
        }
    )
    zone_out = pd.DataFrame(
        {
            "min_lcoe": [1.5, 1.8, 2.0, 2.2, 2.5, 2.7, 3.0, 3.2, 3.5, 3.7, 4.0, 4.2, 4.5, 4.8, 5.0, 5.3],
            "zoning": [0] * 16,
        }
    )
    pd.concat([zone_in, zone_out], ignore_index=True).to_csv(path, index=False)


def _run_discrete_choice(work_dir: Path, run_name: str):
    command = [
        "Rscript",
        "--vanilla",
        str(R_SCRIPT),
        str(work_dir),
        "test",
        "2",
        run_name,
        "0",
        "42",
        *BASELINE_TERMS,
    ]
    subprocess.run(command, check=True, cwd=R_SCRIPT.parent, capture_output=True, text=True)


def test_stage2_same_seed_produces_identical_outputs(tmp_path):
    _ensure_stage2_runtime()
    (tmp_path / "data" / "processed").mkdir(parents=True)
    _write_design_matrix(tmp_path / "data" / "processed" / "dc_data_test.csv")

    _run_discrete_choice(tmp_path, "repro_a")
    _run_discrete_choice(tmp_path, "repro_b")

    coef_a = (tmp_path / "data" / "results" / "spatialdc_coefs_test_repro_a.csv").read_bytes()
    coef_b = (tmp_path / "data" / "results" / "spatialdc_coefs_test_repro_b.csv").read_bytes()
    loglik_a = (tmp_path / "data" / "results" / "spatialdc_loglik_test_repro_a.csv").read_bytes()
    loglik_b = (tmp_path / "data" / "results" / "spatialdc_loglik_test_repro_b.csv").read_bytes()

    assert coef_a == coef_b
    assert loglik_a == loglik_b


def test_stage2_rejects_invalid_seed_early(tmp_path):
    _ensure_stage2_runtime()
    (tmp_path / "data" / "processed").mkdir(parents=True)
    _write_design_matrix(tmp_path / "data" / "processed" / "dc_data_test.csv")

    command = [
        "Rscript",
        "--vanilla",
        str(R_SCRIPT),
        str(tmp_path),
        "test",
        "2",
        "bad_seed",
        "0",
        "not_an_int",
        *BASELINE_TERMS,
    ]
    result = subprocess.run(command, cwd=R_SCRIPT.parent, capture_output=True, text=True)

    assert result.returncode != 0
    assert "seed must be a valid integer" in result.stderr
