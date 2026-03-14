from __future__ import annotations

from pathlib import Path

from config import repo


def _scenario_suffix(scenario: str) -> str:
    return "" if scenario == "2014" else f"_{scenario}"


def _label_segment(output_label: str | None) -> str:
    return "" if output_label is None else f"__{output_label}"


def results_dir() -> Path:
    path = repo / "data" / "results"
    path.mkdir(parents=True, exist_ok=True)
    return path


def processed_dir() -> Path:
    return repo / "data" / "processed"


def figures_dir() -> Path:
    path = repo / "figures"
    path.mkdir(parents=True, exist_ok=True)
    return path


def opt_dir() -> Path:
    path = repo / "opt"
    path.mkdir(parents=True, exist_ok=True)
    return path


def design_matrix_path(scenario: str) -> Path:
    return processed_dir() / f"dc_data_{scenario}.csv"


def wind_atlas_path() -> Path:
    return processed_dir() / "WindAtlas_AUT_pro.nc"


def landscape_atlas_path() -> Path:
    return processed_dir() / "LandscapeAtlas_AUT_pro.nc"


def coefs_path(scenario: str, run_name: str, output_label: str | None = None) -> Path:
    label_segment = _label_segment(output_label)
    if label_segment:
        filename = f"spatialdc_coefs_{scenario}{label_segment}__{run_name}.csv"
    else:
        filename = f"spatialdc_coefs_{scenario}_{run_name}.csv"
    return results_dir() / filename


def loglik_path(scenario: str, run_name: str, output_label: str | None = None) -> Path:
    label_segment = _label_segment(output_label)
    if label_segment:
        filename = f"spatialdc_loglik_{scenario}{label_segment}__{run_name}.csv"
    else:
        filename = f"spatialdc_loglik_{scenario}_{run_name}.csv"
    return results_dir() / filename


def stage2_metadata_path(scenario: str, run_name: str, output_label: str | None = None) -> Path:
    label_segment = _label_segment(output_label)
    if label_segment:
        filename = f"spatialdc_meta_{scenario}{label_segment}__{run_name}.json"
    else:
        filename = f"spatialdc_meta_{scenario}_{run_name}.json"
    return results_dir() / filename


def parsimonious_name_path(scenario: str, output_label: str | None = None) -> Path:
    return repo / f"parsimonious_name{_scenario_suffix(scenario)}{_label_segment(output_label)}.txt"


def soco_path(scenario: str) -> Path:
    return results_dir() / f"soco{_scenario_suffix(scenario)}.nc"


def opt_locations_path(scenario: str, objective: str) -> Path:
    return results_dir() / f"opt_locations_{scenario}_{objective}.csv"


def opt_cost_path(scenario: str, objective: str) -> Path:
    return results_dir() / f"opt_cost_{scenario}_{objective}.csv"


def opt_metadata_path(scenario: str, objective: str) -> Path:
    return results_dir() / f"opt_meta_{scenario}_{objective}.json"


def chunk_locations_cache_path(scenario: str, objective: str, chunk: int) -> Path:
    return opt_dir() / f"chunk_locations_{scenario}_{objective}_{chunk}.parquet"


def chunk_locations_meta_path(scenario: str, objective: str, chunk: int) -> Path:
    return opt_dir() / f"chunk_locations_{scenario}_{objective}_{chunk}.json"


def latex_estimates_path(scenario: str, output_label: str | None = None) -> Path:
    return repo / f"dc_estimates{_scenario_suffix(scenario)}{_label_segment(output_label)}.tex"


def latex_valuation_path(scenario: str, output_label: str | None = None) -> Path:
    return repo / f"valuation{_scenario_suffix(scenario)}{_label_segment(output_label)}.tex"


def latex_summary_stats_path(scenario: str, output_label: str | None = None) -> Path:
    return repo / f"sum_stats{_scenario_suffix(scenario)}{_label_segment(output_label)}.tex"


def latex_cost_comparison_path() -> Path:
    return repo / "cost_comparison.tex"


def map_soco_figure_path(scenario: str) -> Path:
    return figures_dir() / f"map_soco_{scenario}.png"


def map_lcoe_figure_path(scenario: str) -> Path:
    return figures_dir() / f"map_lcoe_{scenario}.png"


def dist_soco_figure_path(scenario: str) -> Path:
    return figures_dir() / f"dist_soco_{scenario}.png"


def soco_facet_figure_path(scenario: str) -> Path:
    return figures_dir() / f"soco_facet{_scenario_suffix(scenario)}.png"


def map_turbines_figure_path(objective: str) -> Path:
    return figures_dir() / f"map_turbines_{objective}.png"


def cost_curve_figure_path() -> Path:
    return figures_dir() / "cost_curve.png"
