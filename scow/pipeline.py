from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from scow.study_registry import REPO_ROOT, StudyDefinition, get_study


@dataclass(frozen=True)
class PipelineStage:
    """One numbered workflow stage for one explicit study."""

    number: int
    script_path: Path
    name: str
    purpose: str
    status: str
    outputs: tuple[str, ...]
    blockers: tuple[str, ...] = ()


@dataclass(frozen=True)
class StudyPipeline:
    """A study-specific view of the current five-stage SCOW workflow."""

    study: StudyDefinition
    stages: tuple[PipelineStage, ...]

    @property
    def executable(self) -> bool:
        return self.study.executable and all(stage.status == "implemented" for stage in self.stages)

    @property
    def blocked_stages(self) -> tuple[PipelineStage, ...]:
        return tuple(stage for stage in self.stages if stage.blockers)


def _script_path(name: str) -> Path:
    return REPO_ROOT / name


def _implemented_pipeline(study: StudyDefinition) -> StudyPipeline:
    return StudyPipeline(
        study=study,
        stages=(
            PipelineStage(
                number=1,
                script_path=_script_path("01_spatial_data.py"),
                name="spatial_data",
                purpose="Build the analysis atlas and export the design matrix.",
                status="implemented",
                outputs=(
                    "data/processed/dc_data_2014.csv",
                    "data/processed/LandscapeAtlas_AUT_pro.nc",
                    "data/processed/WindAtlas_AUT_pro.nc",
                ),
            ),
            PipelineStage(
                number=2,
                script_path=_script_path("02_estimate_wtp.py"),
                name="estimate_wtp",
                purpose="Estimate repeated discrete-choice models and select a parsimonious specification.",
                status="implemented",
                outputs=(
                    "data/results/spatialdc_coefs_2014_<run>.csv",
                    "data/results/spatialdc_loglik_2014_<run>.csv",
                    "parsimonious_name.txt",
                ),
            ),
            PipelineStage(
                number=3,
                script_path=_script_path("03_postprocess_wtp.py"),
                name="postprocess_wtp",
                purpose="Convert coefficients into a local social-cost raster.",
                status="implemented",
                outputs=("data/results/soco.nc",),
            ),
            PipelineStage(
                number=4,
                script_path=_script_path("04_optimal_siting.py"),
                name="optimal_siting",
                purpose="Solve the siting model for lcoe, loco, and soco objectives.",
                status="implemented",
                outputs=(
                    "data/results/opt_locations_2014_lcoe.csv",
                    "data/results/opt_locations_2014_loco.csv",
                    "data/results/opt_locations_2014_soco.csv",
                ),
            ),
            PipelineStage(
                number=5,
                script_path=_script_path("05_postprocess_sites.py"),
                name="postprocess_sites",
                purpose="Compare portfolios and generate maps, figures, and tables.",
                status="implemented",
                outputs=(
                    "figures/map_soco_2014.png",
                    "figures/map_turbines_soco.png",
                    "cost_comparison.tex",
                ),
            ),
        ),
    )


def _wind_2024_pipeline(study: StudyDefinition) -> StudyPipeline:
    stage1_blockers = (
        "The stage-1 payload now uses the official 2024 zoning geometry plus the 2024 zoning-hull and greenland-building sources.",
        "The full 2024 variable inventory is still incomplete, especially moor, RGZ, MLR, radio links, military radar areas, and recreation/tourism layers.",
        "The explicit 2024 namespace is now wired through stages 1 to 5, but the underlying covariate recipe still needs a 2024-specific scientific review.",
    )
    return StudyPipeline(
        study=study,
        stages=(
            PipelineStage(
                number=1,
                script_path=_script_path("01_spatial_data.py"),
                name="spatial_data",
                purpose="Build a 2024-specific wind atlas, zoning outcome, and design matrix.",
                status="implemented",
                outputs=("data/processed/dc_data_lower_austria_wind_2024.csv",),
                blockers=stage1_blockers,
            ),
            PipelineStage(
                number=2,
                script_path=_script_path("02_estimate_wtp.py"),
                name="estimate_wtp",
                purpose="Estimate 2024 wind preferences from the amended decision record.",
                status="implemented",
                outputs=("data/results/spatialdc_coefs_lower_austria_wind_2024_<run>.csv",),
                blockers=(
                    "Stage 2 now accepts the explicit lower_austria_wind_2024 namespace.",
                    "The active variable set is still the current baseline specification and still needs a 2024-specific review.",
                ),
            ),
            PipelineStage(
                number=3,
                script_path=_script_path("03_postprocess_wtp.py"),
                name="postprocess_wtp",
                purpose="Build a 2024 wind local social-cost raster.",
                status="implemented",
                outputs=("data/results/soco_lower_austria_wind_2024.nc",),
                blockers=(
                    "Stage 3 now accepts the explicit lower_austria_wind_2024 namespace.",
                    "Execution still depends on stage-2 outputs for the same explicit study key.",
                ),
            ),
            PipelineStage(
                number=4,
                script_path=_script_path("04_optimal_siting.py"),
                name="optimal_siting",
                purpose="Run the current wind siting model against the 2024 study artifacts.",
                status="implemented",
                outputs=(
                    "data/results/opt_locations_lower_austria_wind_2024_<objective>.csv",
                    "data/results/opt_cost_lower_austria_wind_2024_<objective>.csv",
                ),
                blockers=(
                    "Stage 4 now accepts the explicit lower_austria_wind_2024 namespace.",
                    "Execution still depends on the explicit 2024 stage-3 social-cost artifact.",
                ),
            ),
            PipelineStage(
                number=5,
                script_path=_script_path("05_postprocess_sites.py"),
                name="postprocess_sites",
                purpose="Generate 2024 reporting outputs and portfolio comparisons.",
                status="implemented",
                outputs=(
                    "figures/map_soco_lower_austria_wind_2024.png",
                    "figures/map_turbines_soco_lower_austria_wind_2024.png",
                    "cost_comparison_lower_austria_wind_2024.tex",
                ),
                blockers=(
                    "Stage 5 now accepts the explicit lower_austria_wind_2024 namespace.",
                    "Execution still depends on explicit 2024 stage-3 and stage-4 artifacts.",
                ),
            ),
        ),
    )


def _pv_2022_pipeline(study: StudyDefinition) -> StudyPipeline:
    stage1_blockers = (
        "There is no PV-specific stage-1 implementation path in the current codebase.",
        "The current workflow assumes wind resource, wind turbines, and wind-specific spacing semantics.",
        "Core PV inputs remain missing, especially APG data, military cadastre data, KnollConsult raster inputs, raw 2021 ELT/RGZ screening data, and tourism/soil supplements.",
    )
    return StudyPipeline(
        study=study,
        stages=(
            PipelineStage(
                number=1,
                script_path=_script_path("01_spatial_data.py"),
                name="spatial_data",
                purpose="Build a PV-specific atlas, decision outcome, and design matrix.",
                status="blocked",
                outputs=("data/processed/dc_data_lower_austria_pv_2022.csv",),
                blockers=stage1_blockers,
            ),
            PipelineStage(
                number=2,
                script_path=_script_path("02_estimate_wtp.py"),
                name="estimate_wtp",
                purpose="Estimate PV zoning preferences from the 2022 PV siting decision.",
                status="blocked",
                outputs=("data/results/spatialdc_coefs_lower_austria_pv_2022_<run>.csv",),
                blockers=("Blocked on a PV stage-1 recipe and PV-specific model payload.",),
            ),
            PipelineStage(
                number=3,
                script_path=_script_path("03_postprocess_wtp.py"),
                name="postprocess_wtp",
                purpose="Build a PV local social-cost raster.",
                status="blocked",
                outputs=("data/results/soco_lower_austria_pv_2022.nc",),
                blockers=("Blocked on stage 2 outputs and on PV-specific valuation semantics.",),
            ),
            PipelineStage(
                number=4,
                script_path=_script_path("04_optimal_siting.py"),
                name="optimal_siting",
                purpose="Run a PV-specific siting optimization stage.",
                status="blocked",
                outputs=("data/results/opt_locations_lower_austria_pv_2022_<objective>.csv",),
                blockers=("Blocked on stages 1 to 3 and on PV-specific production and spacing logic.",),
            ),
            PipelineStage(
                number=5,
                script_path=_script_path("05_postprocess_sites.py"),
                name="postprocess_sites",
                purpose="Generate PV reporting outputs and portfolio comparisons.",
                status="blocked",
                outputs=("figures/map_soco_lower_austria_pv_2022.png",),
                blockers=("Blocked on stages 3 and 4.",),
            ),
        ),
    )


def build_pipeline(study_id: str) -> StudyPipeline:
    """Build the current five-stage project pipeline for one explicit study id."""

    study = get_study(study_id)
    if study_id == "lower_austria_wind_2014":
        return _implemented_pipeline(study)
    if study_id == "lower_austria_wind_2024":
        return _wind_2024_pipeline(study)
    if study_id == "lower_austria_pv_2022":
        return _pv_2022_pipeline(study)
    raise KeyError(f"No pipeline builder rule for study_id {study_id!r}")
