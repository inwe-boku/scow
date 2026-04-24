import json
from pathlib import Path

import pandas as pd

from scow import utils


def test_leave_one_out_forwards_seed(monkeypatch):
    calls = []

    def fake_subprocess_rscript(*args, **kwargs):
        calls.append(kwargs["seed"])

    def fake_postprocess_spatialdc(*args, **kwargs):
        coeff = pd.DataFrame(
            {
                "estimate": [1.0],
                "std": [0.1],
                "tstat": [10.0],
                "pval": [0.01],
                "valuation": [1.0],
            },
            index=["min_lcoe"],
        )
        loglik = pd.Series(
            {
                "logLik": -1.0,
                "null.logLik": -2.0,
                "AIC": 10.0,
                "BIC": 11.0,
                "r.squared": 0.1,
                "adj.r.squared": 0.1,
                "nobs": 2.0,
            }
        )
        return coeff, loglik

    monkeypatch.setattr(utils, "_cached_run_matches_terms", lambda *args, **kwargs: False)
    monkeypatch.setattr(utils, "subprocess_rscript", fake_subprocess_rscript)
    monkeypatch.setattr(utils, "postprocess_spatialdc", fake_postprocess_spatialdc)

    utils.leave_one_out(
        variables=("min_lcoe", "distance_roads"),
        integer_variables=(),
        name="exclude_1",
        work_dir=Path("."),
        data_file="2014",
        num_runs="2",
        seed=123,
        compute=True,
        protected_variables=(),
    )

    assert calls == [123, 123]


def test_cached_run_matches_terms_requires_matching_metadata(monkeypatch, tmp_path):
    def fake_coef_path(scenario, run_name, output_label=None):
        return tmp_path / f"spatialdc_coefs_{scenario}_{run_name}.csv"

    def fake_loglik_path(scenario, run_name, output_label=None):
        return tmp_path / f"spatialdc_loglik_{scenario}_{run_name}.csv"

    def fake_meta_path(scenario, run_name, output_label=None):
        return tmp_path / f"spatialdc_meta_{scenario}_{run_name}.json"

    monkeypatch.setattr(utils, "coefs_path", fake_coef_path)
    monkeypatch.setattr(utils, "loglik_path", fake_loglik_path)
    monkeypatch.setattr(utils, "stage2_metadata_path", fake_meta_path)

    pd.DataFrame(
        {
            "term": ["min_lcoe"],
            "estimate": [1.0],
            "std.error": [0.1],
            "statistic": [10.0],
            "p.value": [0.01],
        }
    ).to_csv(fake_coef_path("2014", "full"), index=False)
    pd.DataFrame(
        {
            "logLik": [-1.0],
            "null.logLik": [-2.0],
            "AIC": [10.0],
            "BIC": [11.0],
            "r.squared": [0.1],
            "adj.r.squared": [0.1],
            "nobs": [2.0],
        }
    ).to_csv(fake_loglik_path("2014", "full"), index=False)
    fake_meta_path("2014", "full").write_text(
        json.dumps(
            {
                "artifact_key": "2014",
                "run_name": "full",
                "num_runs": 5,
                "seed": 42,
                "integer_count": 0,
                "variables": ["min_lcoe"],
            }
        )
    )

    assert utils._cached_run_matches_terms(
        "2014",
        "full",
        ["min_lcoe"],
        num_runs=5,
        seed=42,
        integer_count=0,
    )

    fake_meta_path("2014", "full").write_text(
        json.dumps(
            {
                "artifact_key": "2014",
                "run_name": "full",
                "num_runs": 5,
                "seed": 99,
                "integer_count": 0,
                "variables": ["min_lcoe"],
            }
        )
    )

    assert not utils._cached_run_matches_terms(
        "2014",
        "full",
        ["min_lcoe"],
        num_runs=5,
        seed=42,
        integer_count=0,
    )


def test_cached_run_matches_terms_treats_invalid_metadata_as_cache_miss(monkeypatch, tmp_path):
    def fake_coef_path(scenario, run_name, output_label=None):
        return tmp_path / f"spatialdc_coefs_{scenario}_{run_name}.csv"

    def fake_loglik_path(scenario, run_name, output_label=None):
        return tmp_path / f"spatialdc_loglik_{scenario}_{run_name}.csv"

    def fake_meta_path(scenario, run_name, output_label=None):
        return tmp_path / f"spatialdc_meta_{scenario}_{run_name}.json"

    monkeypatch.setattr(utils, "coefs_path", fake_coef_path)
    monkeypatch.setattr(utils, "loglik_path", fake_loglik_path)
    monkeypatch.setattr(utils, "stage2_metadata_path", fake_meta_path)

    pd.DataFrame(
        {
            "term": ["min_lcoe"],
            "estimate": [1.0],
            "std.error": [0.1],
            "statistic": [10.0],
            "p.value": [0.01],
        }
    ).to_csv(fake_coef_path("2014", "full"), index=False)
    pd.DataFrame(
        {
            "logLik": [-1.0],
            "null.logLik": [-2.0],
            "AIC": [10.0],
            "BIC": [11.0],
            "r.squared": [0.1],
            "adj.r.squared": [0.1],
            "nobs": [2.0],
        }
    ).to_csv(fake_loglik_path("2014", "full"), index=False)
    fake_meta_path("2014", "full").write_text("{not valid json")

    assert not utils._cached_run_matches_terms(
        "2014",
        "full",
        ["min_lcoe"],
        num_runs=5,
        seed=42,
        integer_count=0,
    )
