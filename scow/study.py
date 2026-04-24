from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from social_cost_inference.optimization import OptimizationRequest
from social_cost_inference.reporting import ReportingRequest
from social_cost_inference.studies import SpatialDecisionStudy
from scow.study_registry import get_study

AREA_NAME = "Niederösterreich"

SUPPORTED_SCENARIOS = {"2014", "modern"}
SUPPORTED_EXPLICIT_STUDIES = {"lower_austria_wind_2014", "lower_austria_wind_2024"}

BASELINE_BINARY_VARIABLES = (
    "airports_buff",
    "alps_convention",
    "broadleaved",
    "coniferous",
    "important_bird_areas",
    "protected_areas:important_bird_areas",
    "pastures",
    "preservation",
    "protected_areas",
    "restricted_military_areas",
    "water_bodies",
)

BASELINE_FLOAT_VARIABLES = (
    "distance_buildings_in_greenland",
    "distance_existing_turbines",
    "distance_greenland_zonings",
    "distance_other_building_land",
    "distance_power_lines",
    "distance_residential_buildings",
    "distance_roads",
    "elevation",
    "min_lcoe",
    "overnight_stays",
    "slope",
    "tree_cover_density",
    "tree_cover_density:broadleaved",
    "tree_cover_density:coniferous",
)

BASELINE_RESTRICTED_DATA_VARIABLES = (
    "important_bird_areas",
    "overnight_stays",
    "protected_areas:important_bird_areas",
)

BASELINE_PROTECTED_VARIABLES = ("min_lcoe",)

TURBINES_2014 = (
    "Enercon.E40.500",
    "Enercon.E82.3000",
    "Enercon.E101.3050",
    "Enercon.E115.3000",
    "Vestas.V100.1800",
    "Vestas.V100.2000",
    "Vestas.V112.3075",
)
TURBINES_MODERN = (
    "Enercon.E138.3500",
    "Enercon.E160.5560",
    "Vestas.V150.4200",
)

LOSS_FACTOR = 0.71197
ECONOMICS = {
    "discount_rate": 0.05,
    "lifetime_a": 25,
    "om_fixed_eur_per_kw_a": 20.0,
    "om_variable_eur_per_kwh": 0.008,
    "bos_cost_share": 0.3,
    "grid_connect_cost_eur_per_kw": 0.0,
}

OBJECTIVE_COST_THRESHOLDS_EUR_PER_MWH = {
    "loco": 100,
    "lcoe": 75,
    "soco": 125,
}


@dataclass(frozen=True)
class Stage1DesignMatrixSpec:
    """Describe one stage-1 design-matrix export target."""

    name: str
    turbines: tuple[str, ...]


@dataclass(frozen=True)
class Stage1Payload:
    """Study-specific stage-1 vector sources and export targets."""

    zoning_source: Path
    municipal_zoning_source: Path
    buildings_in_greenland_source: Path
    existing_turbines_built_before_year: int
    design_matrix_specs: tuple[Stage1DesignMatrixSpec, ...]
    important_bird_areas_source: Path | None = None
    overnight_stays_source: Path | None = None


@dataclass(frozen=True)
class ScowStudy(SpatialDecisionStudy):
    """Project-local study definition for the current SCOW workflow."""

    scenario: str
    study_id: str | None = None
    spacing: int = 3
    num_turbines: int | str = "auto"
    num_runs: int = 500
    selection_criterion: str = "BIC"
    seed: int = 42
    area_name: str = AREA_NAME
    generation_target: float = 8.0
    slice_count: int = 6
    min_turbines: int = 1200
    max_turbines: int = 6000
    objective_thresholds_eur_per_mwh: dict[str, float] = field(
        default_factory=lambda: dict(OBJECTIVE_COST_THRESHOLDS_EUR_PER_MWH)
    )

    def __post_init__(self) -> None:
        if self.scenario not in SUPPORTED_SCENARIOS:
            raise ValueError(f"Unsupported scenario: {self.scenario}")
        if self.study_id is not None:
            if self.study_id not in SUPPORTED_EXPLICIT_STUDIES:
                raise ValueError(f"Unsupported explicit study id for ScowStudy: {self.study_id}")
            expected_scenario = {
                "lower_austria_wind_2014": "2014",
                "lower_austria_wind_2024": "modern",
            }[self.study_id]
            if self.scenario != expected_scenario:
                raise ValueError(
                    f"Explicit study {self.study_id} requires scenario {expected_scenario!r}, "
                    f"got {self.scenario!r}"
                )

    @classmethod
    def from_study_id(cls, study_id: str, **overrides) -> "ScowStudy":
        """Build a wind-study runtime configuration from one explicit study id."""

        if study_id == "lower_austria_wind_2014":
            return cls(scenario="2014", study_id=study_id, **overrides)
        if study_id == "lower_austria_wind_2024":
            return cls(scenario="modern", study_id=study_id, **overrides)
        raise NotImplementedError(
            f"ScowStudy currently supports explicit wind studies only, got {study_id!r}"
        )

    @property
    def turbines(self) -> tuple[str, ...]:
        if self.scenario == "2014":
            return TURBINES_2014
        return TURBINES_MODERN

    @property
    def resolved_study_id(self) -> str | None:
        """Return the explicit study id when one is selected."""

        return self.study_id

    @property
    def artifact_key(self) -> str:
        """Return the shared artifact namespace for this study."""

        if self.study_id is None:
            return self.scenario
        if self.study_id == "lower_austria_wind_2014":
            return "2014"
        return self.study_id

    @property
    def social_cost_key(self) -> str:
        """Compatibility alias for the shared artifact namespace."""

        return self.artifact_key

    @property
    def stage4_key(self) -> str:
        """Compatibility alias for the shared artifact namespace."""

        return self.artifact_key

    @property
    def loss_factor(self) -> float:
        return LOSS_FACTOR

    @property
    def economics(self) -> dict[str, float]:
        return dict(ECONOMICS)

    @property
    def objective_family(self) -> tuple[str, ...]:
        return tuple(self.objective_thresholds_eur_per_mwh.keys())

    def objective_threshold(self, objective: str) -> float:
        return self.objective_thresholds_eur_per_mwh[objective]

    def stage1_payload(self, repo_root: Path) -> Stage1Payload:
        """Return the stage-1 vector payload and export plan for this study."""

        site_root = repo_root / "data" / "site"
        if self.study_id == "lower_austria_wind_2024":
            study = get_study(self.study_id)
            assert study.source_documents_root is not None
            assert study.source_data_root is not None
            return Stage1Payload(
                zoning_source=study.source_documents_root / "official_wind_zones_2024_epsg4326.zip",
                municipal_zoning_source=study.source_data_root
                / "03_distance_settlement_inputs"
                / "noe_widmungshuellen_rru_wi_huelle_epsg31259.zip",
                buildings_in_greenland_source=study.source_data_root
                / "03_distance_settlement_inputs"
                / "noe_erhaltenswerte_gebaeude_rru_wi_geb_epsg31259.zip",
                important_bird_areas_source=None,
                overnight_stays_source=None,
                existing_turbines_built_before_year=study.decision_year + 1,
                design_matrix_specs=(
                    Stage1DesignMatrixSpec("lower_austria_wind_2024", self.turbines),
                ),
            )

        return Stage1Payload(
            zoning_source=site_root / "zones" / "Zonierung_noe.shp",
            municipal_zoning_source=site_root / "landuse" / "Widmungen_noe.shp",
            buildings_in_greenland_source=site_root / "landuse" / "Erhaltenswerte_gebaeude.shp",
            important_bird_areas_source=site_root / "iba" / "Europe_KBA.shp",
            overnight_stays_source=site_root / "tourism" / "Naechtigungen.shp",
            existing_turbines_built_before_year=2015,
            design_matrix_specs=(
                Stage1DesignMatrixSpec("2014", TURBINES_2014),
                Stage1DesignMatrixSpec("modern", TURBINES_MODERN),
            ),
        )

    def estimation_spec(self):
        from scow.estimation import EstimationSpec

        return EstimationSpec(
            scenario=self.scenario,
            method="mnl_logitr",
            search_mode="gts_leave_one_out",
            num_runs=self.num_runs,
            selection_criterion=self.selection_criterion,
            seed=self.seed,
            binary_variables=BASELINE_BINARY_VARIABLES,
            float_variables=BASELINE_FLOAT_VARIABLES,
            protected_variables=BASELINE_PROTECTED_VARIABLES,
            restricted_data_variables=BASELINE_RESTRICTED_DATA_VARIABLES,
            artifact_key=self.artifact_key,
            tag=None,
        )

    def optimization_request(self, objective: str) -> OptimizationRequest:
        """Return the current stage-4 request for one named objective."""

        return OptimizationRequest(
            scenario=self.scenario,
            artifact_key=self.artifact_key,
            objective=objective,
            spacing=self.spacing,
            num_turbines=self.num_turbines,
            slice_count=self.slice_count,
            min_turbines=self.min_turbines,
            max_turbines=self.max_turbines,
            cost_threshold=self.objective_threshold(objective),
        )

    def reporting_request(self) -> ReportingRequest:
        """Return the current stage-5 reporting request for this study."""

        return ReportingRequest(
            scenario=self.scenario,
            artifact_key=self.artifact_key,
            area_name=self.area_name,
            generation_target=self.generation_target,
        )
