from argparse import Namespace

import pytest

from scow.estimation_cli import build_spec_from_args
from scow.pipeline import build_pipeline
from scow.social_cost import get_social_cost_atlas
from scow.study_registry import get_study, study_ids
from scow.study import ScowStudy


def test_study_registry_exposes_explicit_case_studies():
    assert study_ids() == (
        "lower_austria_pv_2022",
        "lower_austria_wind_2014",
        "lower_austria_wind_2024",
    )


def test_wind_2014_is_the_only_executable_pipeline():
    pipeline = build_pipeline("lower_austria_wind_2014")

    assert pipeline.study.executable is True
    assert pipeline.executable is True
    assert all(stage.status == "implemented" for stage in pipeline.stages)


def test_wind_2024_pipeline_has_stage1_payload_but_is_not_executable():
    pipeline = build_pipeline("lower_austria_wind_2024")

    assert pipeline.study.status == "stage5_namespace_implemented"
    assert pipeline.executable is False
    assert all(stage.status == "implemented" for stage in pipeline.stages)
    assert [stage.name for stage in pipeline.blocked_stages] == [
        "spatial_data",
        "estimate_wtp",
        "postprocess_wtp",
        "optimal_siting",
        "postprocess_sites",
    ]


def test_pv_2022_pipeline_is_blocked_on_missing_payload():
    pipeline = build_pipeline("lower_austria_pv_2022")
    study = get_study("lower_austria_pv_2022")

    assert study.executable is False
    assert pipeline.executable is False
    assert pipeline.stages[0].status == "blocked"
    assert "APG grid data delivery used in the PV zoning process" in study.data_state.major_missing_inputs


def test_scow_study_from_study_id_builds_2024_stage1_payload(tmp_path):
    study = ScowStudy.from_study_id("lower_austria_wind_2024")
    payload = study.stage1_payload(tmp_path)

    assert study.scenario == "modern"
    assert payload.zoning_source.name == "official_wind_zones_2024_epsg4326.zip"
    assert payload.municipal_zoning_source.name == "noe_widmungshuellen_rru_wi_huelle_epsg31259.zip"
    assert payload.buildings_in_greenland_source.name == "noe_erhaltenswerte_gebaeude_rru_wi_geb_epsg31259.zip"
    assert payload.existing_turbines_built_before_year == 2025
    assert payload.design_matrix_specs[0].name == "lower_austria_wind_2024"


def test_legacy_scenario_mode_keeps_supported_2014_and_modern_artifact_names():
    baseline = ScowStudy(scenario="2014")
    modern = ScowStudy(scenario="modern")

    assert baseline.artifact_key == "2014"
    assert modern.artifact_key == "modern"


def test_estimation_specs_keep_2014_namespace_and_use_explicit_2024_namespace():
    study_2014 = ScowStudy.from_study_id("lower_austria_wind_2014")
    study_2024 = ScowStudy.from_study_id("lower_austria_wind_2024")

    spec_2014 = study_2014.estimation_spec()
    spec_2024 = study_2024.estimation_spec()
    legacy_modern = ScowStudy(scenario="modern").estimation_spec()

    assert spec_2014.scenario == "2014"
    assert spec_2014.artifact_key == "2014"
    assert spec_2024.scenario == "modern"
    assert spec_2024.artifact_key == "lower_austria_wind_2024"
    assert legacy_modern.scenario == "modern"
    assert legacy_modern.artifact_key == "modern"


def test_study_artifact_keys_keep_legacy_and_explicit_names_distinct():
    study_2014 = ScowStudy.from_study_id("lower_austria_wind_2014")
    study_2024 = ScowStudy.from_study_id("lower_austria_wind_2024")
    modern = ScowStudy(scenario="modern")

    assert study_2014.artifact_key == "2014"
    assert study_2024.artifact_key == "lower_austria_wind_2024"
    assert modern.artifact_key == "modern"


def test_social_cost_artifact_handle_keeps_runtime_scenario_and_artifact_key_split():
    atlas = get_social_cost_atlas(ScowStudy.from_study_id("lower_austria_wind_2024"))

    assert atlas.scenario == "modern"
    assert atlas.artifact_key == "lower_austria_wind_2024"
    assert atlas.netcdf_path.name == "soco_lower_austria_wind_2024.nc"


def test_stage4_and_stage5_requests_keep_runtime_scenario_and_split_artifact_key():
    study_2014 = ScowStudy.from_study_id("lower_austria_wind_2014")
    study_2024 = ScowStudy.from_study_id("lower_austria_wind_2024")

    optimization_request_2014 = study_2014.optimization_request("soco")
    optimization_request_2024 = study_2024.optimization_request("soco")
    reporting_request_2014 = study_2014.reporting_request()
    reporting_request_2024 = study_2024.reporting_request()

    assert optimization_request_2014.scenario == "2014"
    assert optimization_request_2014.artifact_key == "2014"
    assert optimization_request_2024.scenario == "modern"
    assert optimization_request_2024.artifact_key == "lower_austria_wind_2024"
    assert reporting_request_2014.scenario == "2014"
    assert reporting_request_2014.artifact_key == "2014"
    assert reporting_request_2024.scenario == "modern"
    assert reporting_request_2024.artifact_key == "lower_austria_wind_2024"


def test_build_spec_from_args_rejects_conflicting_explicit_study_scenario_override():
    study = ScowStudy.from_study_id("lower_austria_wind_2024")

    with pytest.raises(ValueError, match="conflicting override"):
        build_spec_from_args(
            study,
            Namespace(
                preset=None,
                scenario="2014",
                method=None,
                num_runs=None,
                selection_criterion=None,
                seed=None,
                single=False,
                full_significant=False,
                gts_leave_one_out=False,
            ),
        )


def test_build_spec_from_args_allows_matching_explicit_study_scenario_override():
    study = ScowStudy.from_study_id("lower_austria_wind_2024")

    spec = build_spec_from_args(
        study,
        Namespace(
            preset=None,
            scenario="modern",
            method=None,
            num_runs=None,
            selection_criterion=None,
            seed=None,
            single=False,
            full_significant=False,
            gts_leave_one_out=False,
        ),
    )

    assert spec.scenario == "modern"
    assert spec.artifact_key == "lower_austria_wind_2024"
