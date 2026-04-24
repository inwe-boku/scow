from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class StudyDataState:
    """Summarize current local data readiness for one explicit study."""

    summary: str
    ready_count: int = 0
    partial_count: int = 0
    missing_count: int = 0
    restricted_inputs: tuple[str, ...] = ()
    major_missing_inputs: tuple[str, ...] = ()


@dataclass(frozen=True)
class StudyDefinition:
    """Tracked project-level definition for one case-study target."""

    study_id: str
    title: str
    technology: str
    region: str
    decision_year: int
    legal_reference: str
    decision_description: str
    status: str
    executable: bool
    source_documents_root: Path | None = None
    source_data_root: Path | None = None
    current_runtime_scenario: str | None = None
    paper_reference: str | None = None
    notes: tuple[str, ...] = ()
    data_state: StudyDataState = field(default_factory=lambda: StudyDataState(summary="unknown"))


_STUDIES = {
    "lower_austria_wind_2014": StudyDefinition(
        study_id="lower_austria_wind_2014",
        title="Lower Austria wind zoning 2014",
        technology="wind",
        region="Lower Austria",
        decision_year=2014,
        legal_reference="LGBl. 8001/1-0",
        decision_description="Baseline wind-zoning study underlying the published SCOW paper.",
        status="implemented",
        executable=True,
        source_documents_root=REPO_ROOT / ".local" / "zoning_documents",
        source_data_root=REPO_ROOT / "data",
        current_runtime_scenario="2014",
        paper_reference="Wehrle, Regner, Morawetz, Schmidt (2024), Energy Economics 139, 107865",
        notes=(
            "This is the only study currently implemented end-to-end in the numbered pipeline.",
            "The current executable baseline is a reduced-data local baseline, not a pristine full-data rebuild.",
            "The original 2014 ordinance bundle is not yet ingested as a dedicated local study bundle.",
        ),
        data_state=StudyDataState(
            summary="Executable baseline artifacts are present locally, with restricted-data substitutions still active.",
            restricted_inputs=(
                "important_bird_areas",
                "overnight_stays",
                "protected_areas:important_bird_areas",
            ),
        ),
    ),
    "lower_austria_pv_2022": StudyDefinition(
        study_id="lower_austria_pv_2022",
        title="Lower Austria PV zoning 2022",
        technology="pv",
        region="Lower Austria",
        decision_year=2022,
        legal_reference="LGBl. Nr. 94/2022",
        decision_description="Ground-mounted PV zoning program in greenland.",
        status="documents_ingested",
        executable=False,
        source_documents_root=REPO_ROOT / ".local" / "zoning_documents" / "noe_pv_2022",
        notes=(
            "Public core documents are ingested locally.",
            "There is no dedicated stage-1 data bundle yet.",
            "The current numbered pipeline has no PV-specific implementation path.",
        ),
        data_state=StudyDataState(
            summary="Documents are present, but study-specific GIS inputs remain incomplete and code is not implemented.",
            major_missing_inputs=(
                "APG grid data delivery used in the PV zoning process",
                "Military cadastre layers used in the process",
                "Bundesstrassenplanungen",
                "Leitungsinfrastrukturen",
                "KnollConsult green-space evaluation raster",
                "Raw 2021 ELT/RGZ screening data",
                "Tourism and excursion GIS layer",
                "AGES soil supplement for the nine eBod-missing municipalities",
                "Ad hoc municipal zoning-plan GIS inputs",
            ),
        ),
    ),
    "lower_austria_wind_2024": StudyDefinition(
        study_id="lower_austria_wind_2024",
        title="Lower Austria wind zoning amendment 2024",
        technology="wind",
        region="Lower Austria",
        decision_year=2024,
        legal_reference="LGBl. Nr. 47/2024",
        decision_description="First amendment to the Lower Austria wind zoning program.",
        status="stage5_namespace_implemented",
        executable=False,
        source_documents_root=REPO_ROOT / ".local" / "zoning_documents" / "noe_wind_2024",
        source_data_root=REPO_ROOT / ".local" / "zoning_data" / "noe_wind_2024",
        current_runtime_scenario="modern",
        notes=(
            "This is the best-prepared follow-on study.",
            "A dedicated 2024 raw-data bundle exists locally.",
            "The stage-1 payload is now wired into the live spatial-preprocessing path.",
            "Stage 2 now accepts the explicit lower_austria_wind_2024 namespace.",
            "Stage 3 now accepts the explicit lower_austria_wind_2024 namespace.",
            "Stages 4 and 5 now accept the explicit lower_austria_wind_2024 namespace.",
            "The study is still not marked executable because the 2024 covariate recipe and missing-input inventory are not yet closed out.",
        ),
        data_state=StudyDataState(
            summary="Documents and a dedicated raw-data bundle are present, and the explicit stage-1 payload is wired into the spatial-preprocessing code.",
            ready_count=35,
            partial_count=15,
            missing_count=12,
            major_missing_inputs=(
                "moor",
                "regional_green_zone",
                "multifunctional_landscape",
                "radio_links",
                "military_radar_areas",
                "recreation_and_health_tourism_centers",
            ),
        ),
    ),
}


def list_studies() -> tuple[StudyDefinition, ...]:
    """Return the explicit case studies tracked by the project."""

    return tuple(_STUDIES[study_id] for study_id in sorted(_STUDIES))


def get_study(study_id: str) -> StudyDefinition:
    """Return one tracked study definition by explicit study id."""

    try:
        return _STUDIES[study_id]
    except KeyError as exc:
        known = ", ".join(sorted(_STUDIES))
        raise KeyError(f"Unknown study_id {study_id!r}. Known studies: {known}") from exc


def study_ids() -> tuple[str, ...]:
    """Return all explicit study ids."""

    return tuple(sorted(_STUDIES))
