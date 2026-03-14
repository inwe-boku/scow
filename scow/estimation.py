from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from config import repo
from scow.paths import (
    coefs_path,
    design_matrix_path,
    latex_estimates_path,
    latex_summary_stats_path,
    latex_valuation_path,
    parsimonious_name_path,
)
from scow.tex_convert import (
    coeffs_to_latex,
    convert_names,
    stats_to_latex,
    summary_statistics_to_latex,
    units,
)
from scow.utils import leave_one_out, postprocess_spatialdc, subprocess_rscript

SUPPORTED_SCENARIOS = {"2014", "modern"}
SUPPORTED_METHODS = {"mnl_logitr"}
SUPPORTED_SEARCH_MODES = {"single", "full_significant", "gts_leave_one_out"}
SUPPORTED_SELECTION_CRITERIA = {"AIC", "BIC"}

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


def _dedupe_preserve_order(values: tuple[str, ...] | list[str]) -> tuple[str, ...]:
    ordered = dict.fromkeys(values)
    return tuple(ordered.keys())


def _artifact_run_name(run_name: str) -> str:
    return run_name.replace(":", "_")


@dataclass(frozen=True)
class EstimationSpec:
    """Stage-2 estimation configuration.

    The spec captures the canonical stage-2 choices: scenario, econometric method,
    model-selection mode, repeated-run count, variable lists, and reporting-related
    metadata. It is intentionally Python-first and supports immutable-style
    semantic overrides such as `.single()` or `.with_num_runs(...)`.
    """

    scenario: str
    method: str
    search_mode: str
    num_runs: int
    selection_criterion: str
    seed: int
    binary_variables: tuple[str, ...]
    float_variables: tuple[str, ...]
    protected_variables: tuple[str, ...]
    restricted_data_variables: tuple[str, ...]
    output_label: str | None = None
    tag: str | None = None

    def __post_init__(self) -> None:
        if self.scenario not in SUPPORTED_SCENARIOS:
            raise ValueError(f"Unsupported scenario: {self.scenario}")
        if self.method not in SUPPORTED_METHODS:
            raise NotImplementedError(
                f"Method {self.method!r} is not implemented for canonical stage-2 execution."
            )
        if self.search_mode not in SUPPORTED_SEARCH_MODES:
            raise ValueError(f"Unsupported search mode: {self.search_mode}")
        if self.selection_criterion not in SUPPORTED_SELECTION_CRITERIA:
            raise ValueError(f"Unsupported selection criterion: {self.selection_criterion}")
        if self.num_runs <= 0:
            raise ValueError("num_runs must be a positive integer")

        binary_variables = _dedupe_preserve_order(self.binary_variables)
        float_variables = _dedupe_preserve_order(self.float_variables)
        protected_variables = _dedupe_preserve_order(self.protected_variables)
        restricted_data_variables = _dedupe_preserve_order(self.restricted_data_variables)

        overlap = set(binary_variables) & set(float_variables)
        if overlap:
            raise ValueError(f"Binary and float variables must not overlap: {sorted(overlap)}")
        if len(binary_variables) != len(self.binary_variables):
            raise ValueError("binary_variables contain duplicates")
        if len(float_variables) != len(self.float_variables):
            raise ValueError("float_variables contain duplicates")
        if "min_lcoe" not in float_variables:
            raise ValueError("min_lcoe must remain part of float_variables")
        if "min_lcoe" not in protected_variables:
            raise ValueError("min_lcoe must remain a protected variable")
        missing_protected = set(protected_variables) - set(binary_variables) - set(float_variables)
        if missing_protected:
            raise ValueError(
                f"Protected variables are not present in the model variable set: {sorted(missing_protected)}"
            )

    @property
    def variables(self) -> tuple[str, ...]:
        return self.binary_variables + self.float_variables

    @property
    def integer_count(self) -> int:
        return len(self.binary_variables)

    @property
    def suffix(self) -> str:
        return "" if self.scenario == "2014" else f"_{self.scenario}"

    @property
    def parsimonious_name_file(self) -> Path:
        return parsimonious_name_path(self.scenario, output_label=self.output_label)

    @property
    def canonical_output_namespace(self) -> bool:
        return self.output_label is None

    def _updated(self, **changes: object) -> "EstimationSpec":
        values = {
            "scenario": self.scenario,
            "method": self.method,
            "search_mode": self.search_mode,
            "num_runs": self.num_runs,
            "selection_criterion": self.selection_criterion,
            "seed": self.seed,
            "binary_variables": self.binary_variables,
            "float_variables": self.float_variables,
            "protected_variables": self.protected_variables,
            "restricted_data_variables": self.restricted_data_variables,
            "output_label": self.output_label,
            "tag": self.tag,
        }
        values.update(changes)
        return EstimationSpec(**values)

    def for_scenario(self, scenario: str) -> "EstimationSpec":
        return self._updated(scenario=scenario)

    def using_method(self, method: str) -> "EstimationSpec":
        return self._updated(method=method)

    def single(self) -> "EstimationSpec":
        return self._updated(search_mode="single")

    def full_significant(self) -> "EstimationSpec":
        return self._updated(search_mode="full_significant")

    def gts_leave_one_out(self) -> "EstimationSpec":
        return self._updated(search_mode="gts_leave_one_out")

    def with_num_runs(self, num_runs: int) -> "EstimationSpec":
        return self._updated(num_runs=int(num_runs))

    def with_selection_criterion(self, selection_criterion: str) -> "EstimationSpec":
        return self._updated(selection_criterion=selection_criterion)

    def with_seed(self, seed: int) -> "EstimationSpec":
        return self._updated(seed=int(seed))

    def exclude_variables(self, variables: list[str] | tuple[str, ...]) -> "EstimationSpec":
        excluded = set(variables)
        protected_overlap = excluded & set(self.protected_variables)
        if protected_overlap:
            raise ValueError(
                f"Cannot exclude protected variables: {sorted(protected_overlap)}"
            )
        return self._updated(
            binary_variables=tuple(v for v in self.binary_variables if v not in excluded),
            float_variables=tuple(v for v in self.float_variables if v not in excluded),
        )

    def protect_variables(self, variables: list[str] | tuple[str, ...]) -> "EstimationSpec":
        return self._updated(
            protected_variables=_dedupe_preserve_order(self.protected_variables + tuple(variables))
        )

    def with_tag(self, tag: str | None) -> "EstimationSpec":
        return self._updated(tag=tag)

    def with_output_label(self, output_label: str | None) -> "EstimationSpec":
        return self._updated(output_label=output_label)


@dataclass(frozen=True)
class EstimationResult:
    """Result metadata for one stage-2 estimation run."""

    spec: EstimationSpec
    executed_run_names: tuple[str, ...]
    selected_run_name: str
    selected_criterion_value: float | None
    removed_zero_variance: tuple[str, ...]
    removed_restricted_zero_variance: tuple[str, ...]
    surviving_binary_variables: tuple[str, ...]
    surviving_float_variables: tuple[str, ...]
    generated_reports: tuple[str, ...]

    @property
    def run_names(self) -> tuple[str, ...]:
        return self.executed_run_names

    @property
    def parsimonious_name(self) -> str:
        return self.selected_run_name

    def coef_path(self, run_name: str) -> Path:
        return coefs_path(
            self.spec.scenario,
            _artifact_run_name(run_name),
            output_label=self.spec.output_label,
        )

    def loglik_path(self, run_name: str) -> Path:
        return loglik_path(
            self.spec.scenario,
            _artifact_run_name(run_name),
            output_label=self.spec.output_label,
        )

    @property
    def parsimonious_coef_path(self) -> Path:
        return self.coef_path(self.selected_run_name)

    @property
    def parsimonious_loglik_path(self) -> Path:
        return self.loglik_path(self.selected_run_name)

    def summary(self) -> str:
        return (
            f"scenario={self.spec.scenario}, method={self.spec.method}, "
            f"search_mode={self.spec.search_mode}, selected={self.selected_run_name}, "
            f"runs={len(self.executed_run_names)}"
        )


class _PresetNamespace:
    """Namespace of stage-2 preset factories."""

    def baseline(self) -> EstimationSpec:
        return EstimationSpec(
            scenario="2014",
            method="mnl_logitr",
            search_mode="gts_leave_one_out",
            num_runs=500,
            selection_criterion="BIC",
            seed=42,
            binary_variables=BASELINE_BINARY_VARIABLES,
            float_variables=BASELINE_FLOAT_VARIABLES,
            protected_variables=BASELINE_PROTECTED_VARIABLES,
            restricted_data_variables=BASELINE_RESTRICTED_DATA_VARIABLES,
            tag=None,
        )

    def single_check(self) -> EstimationSpec:
        return self.baseline().single()

    def fast_debug(self) -> EstimationSpec:
        return self.baseline().single().with_num_runs(25)

    def modern_baseline(self) -> EstimationSpec:
        return self.baseline().for_scenario("modern")


presets = _PresetNamespace()


def model_term_series(dataframe: pd.DataFrame, term: str) -> pd.Series | None:
    """Return the data series for a model term, including interaction terms."""
    if term in dataframe.columns:
        return dataframe[term]

    if ":" not in term:
        return None

    parts = term.split(":")
    if not all(part in dataframe.columns for part in parts):
        return None

    series = dataframe[parts[0]]
    for part in parts[1:]:
        series = series * dataframe[part]
    return series


def drop_zero_variance_variables(
    dataframe: pd.DataFrame,
    integer_variables: tuple[str, ...],
    float_variables: tuple[str, ...],
) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
    """Return surviving variables after removing terms with no variation."""
    surviving_integers: list[str] = []
    surviving_floats: list[str] = []
    removed: list[str] = []

    for term in integer_variables + float_variables:
        series = model_term_series(dataframe, term)
        if series is None:
            raise KeyError(f"Model term {term!r} is not available in the design matrix.")
        if series.dropna().nunique() <= 1:
            removed.append(term)
        elif term in integer_variables:
            surviving_integers.append(term)
        else:
            surviving_floats.append(term)

    return tuple(surviving_integers), tuple(surviving_floats), tuple(removed)


def eliminate_variable(
    variables_list: tuple[str, ...] | list[str],
    eliminate: tuple[str, ...] | list[str],
    integer_variables: tuple[str, ...] | list[str],
) -> tuple[tuple[str, ...], tuple[str, ...], int]:
    """Drop selected variables and return the surviving integer/float partitions."""
    eliminate_with_interactions = set(eliminate)
    for elim_var in eliminate:
        if ":" in elim_var:
            eliminate_with_interactions.update([var for var in variables_list if elim_var in var])
        else:
            eliminate_with_interactions.add(elim_var)

    surviving_integers = tuple(
        ivar
        for ivar in variables_list
        if ivar in integer_variables and ivar not in eliminate_with_interactions
    )
    surviving_floats = tuple(
        var
        for var in variables_list
        if var not in eliminate_with_interactions and var not in integer_variables
    )
    return surviving_integers, surviving_floats, len(surviving_integers)


def _estimate_spatial_dc(
    spec: EstimationSpec,
    run_name: str,
    integer_count: int,
    variables: tuple[str, ...] | list[str],
    compute: bool,
) -> tuple[pd.DataFrame, float, list[str], pd.Series]:
    """Run one spatial discrete-choice specification and postprocess it."""
    artifact_run_name = _artifact_run_name(run_name)
    if compute:
        subprocess_rscript(
            str(repo),
            spec.scenario,
            str(spec.num_runs),
            artifact_run_name,
            str(integer_count),
            list(variables),
            seed=spec.seed,
            output_label=spec.output_label,
        )

    coeff, loglik = postprocess_spatialdc(
        repo,
        spec.scenario,
        artifact_run_name,
        output_label=spec.output_label,
    )
    insignificant = list(coeff[coeff["pval"] > 0.1].index)
    criterion = loglik[spec.selection_criterion]
    if pd.isna(criterion):
        raise ValueError(
            f"{spec.selection_criterion} is NA for run {artifact_run_name}. "
            "Check whether zero-variance or restricted-data terms are still entering the model."
        )
    return coeff, float(criterion), insignificant, loglik


def _leave_one_out_run_names(
    variables: tuple[str, ...] | list[str],
    protected_variables: tuple[str, ...] | list[str],
    prefix: str,
) -> tuple[str, ...]:
    removable_variables = [var for var in variables if var not in set(protected_variables)]
    return tuple(_artifact_run_name(f"{prefix}_{var}") for var in removable_variables)


def _write_parsimonious_name(spec: EstimationSpec, parsimonious_name: str) -> None:
    spec.parsimonious_name_file.write_text(parsimonious_name)


def _write_estimates_report(
    spec: EstimationSpec,
    full_coefs: pd.DataFrame,
    full_stats: pd.Series,
    parsimonious_coeffs: pd.DataFrame,
    parsimonious_stats: pd.Series,
) -> Path:
    latex_coeffs = coeffs_to_latex(full_coefs, parsimonious_coeffs, 5)
    latex_stats = stats_to_latex(full_stats, parsimonious_stats)
    outpath = latex_estimates_path(spec.scenario, output_label=spec.output_label)
    outpath.write_text(latex_coeffs + "\n\\midrule\n" + latex_stats)
    return outpath


def _write_valuation_report(spec: EstimationSpec, parsimonious_name: str, parsimonious_coeffs: pd.DataFrame) -> Path:
    ests = pd.read_csv(
        coefs_path(spec.scenario, parsimonious_name, output_label=spec.output_label)
    )
    ests = ests.set_index("term")

    wtp_confidence = []
    for mod_nr in range(1, spec.num_runs + 1):
        run_estimates = ests[ests["mod"] == mod_nr]["estimate"]
        wtp_confidence.append(run_estimates / run_estimates.loc["min_lcoe"])

    wtp_confidence = pd.concat(wtp_confidence, axis=1)
    wtp = pd.DataFrame(index=wtp_confidence.index, columns=["Estimate", "ci95_lower", "ci95_upper"])
    wtp["ci95_lower"] = wtp_confidence.quantile(0.025, axis=1)
    wtp["ci95_upper"] = wtp_confidence.quantile(0.975, axis=1)
    wtp["Estimate"] = parsimonious_coeffs["estimate"] / parsimonious_coeffs.loc["min_lcoe", "estimate"]
    wtp["95% confidence interval"] = wtp.apply(
        lambda row: f"[{row['ci95_lower']:.3f}, {row['ci95_upper']:.3f}]",
        axis=1,
    )
    wtp = wtp.drop(columns=["ci95_lower", "ci95_upper"])

    units_wtp = {
        key: r"\euro/\si{\mega\watt\hour} / " + value
        if value not in ["--", r"\si{\percent}"]
        else r"\euro/\si{\mega\watt\hour}"
        for key, value in units.items()
    }
    wtp.loc["coniferous:tree_cover_density", "Unit"] = r"\euro/\si{\mega\watt\hour}"

    wtp = wtp.rename(index=convert_names)
    wtp["Unit"] = wtp.index.map(units_wtp)
    wtp.loc["Quasi-LCOE", "Unit"] = ""
    wtp.loc["Touristic overnights", "Unit"] = r"\euro/\si{\mega\watt\hour}" + " / 1000"
    wtp.loc[wtp["Unit"].isna(), "Unit"] = r"\euro/\si{\mega\watt\hour}"
    wtp = wtp[["Unit", "Estimate", "95% confidence interval"]]
    latex_wtp = wtp.to_latex(float_format="%.3f")
    lines = latex_wtp.split("\n")[5:-3]
    if lines[-1].endswith("\\\\"):
        lines[-1] = lines[-1][:-2]
    outpath = latex_valuation_path(spec.scenario, output_label=spec.output_label)
    outpath.write_text("\n".join(lines))
    return outpath


def run_estimation(spec: EstimationSpec, compute: bool = True, report: bool = True) -> EstimationResult:
    """Run canonical stage-2 estimation and optional reporting.

    This first implementation slice supports only the current canonical
    `mnl_logitr` path and preserves the existing SCOW artifact namespace.
    """

    if spec.method != "mnl_logitr":
        raise NotImplementedError(
            f"Canonical stage-2 execution is currently implemented only for {SUPPORTED_METHODS}."
        )
    if compute and (
        spec.method != "mnl_logitr"
        or spec.search_mode != "gts_leave_one_out"
        or spec.output_label is not None
    ):
        if spec.output_label is None and spec.search_mode != "gts_leave_one_out":
            raise ValueError(
                "Noncanonical stage-2 compute runs require an explicit output_label to avoid "
                "overwriting canonical artifacts."
            )

    results_dir = coefs_path(spec.scenario, "full", output_label=spec.output_label).parent
    results_dir.mkdir(parents=True, exist_ok=True)

    design_matrix = pd.read_csv(design_matrix_path(spec.scenario))
    generated_reports: list[str] = []

    if report:
        summary_statistics_to_latex(
            design_matrix,
            str(latex_summary_stats_path(spec.scenario, output_label=spec.output_label)),
        )
        generated_reports.append(
            latex_summary_stats_path(spec.scenario, output_label=spec.output_label).name
        )

    surviving_binary, surviving_float, removed_zero_variance = drop_zero_variance_variables(
        design_matrix,
        spec.binary_variables,
        spec.float_variables,
    )
    surviving_variables = surviving_binary + surviving_float
    removed_restricted = tuple(
        term for term in removed_zero_variance if term in set(spec.restricted_data_variables)
    )

    if removed_zero_variance:
        print(f"Excluding zero-variance model terms: {', '.join(removed_zero_variance)}")
    if removed_restricted:
        print(
            "NOTICE: restricted-data substitutions are active in this checkout. "
            f"The following restricted-data terms are excluded from estimation because they have zero variance: "
            f"{', '.join(removed_restricted)}"
        )
    if not surviving_variables:
        raise ValueError("All model variables were removed by the zero-variance guard.")

    executed_run_names: list[str] = []

    full_coefs, full_criterion, insignificant, full_stats = _estimate_spatial_dc(
        spec,
        "full",
        len(surviving_binary),
        surviving_variables,
        compute=compute,
    )
    executed_run_names.append("full")

    if spec.search_mode == "single":
        parsimonious_name = "full"
        parsimonious_coeffs = full_coefs
        parsimonious_stats = full_stats
        selected_criterion_value = full_criterion
    else:
        significant_binary, significant_float, significant_int_count = eliminate_variable(
            surviving_variables,
            tuple(insignificant),
            surviving_binary,
        )
        significant_variables = significant_binary + significant_float

        parsimonious_coeffs, selected_criterion_value, _, parsimonious_stats = _estimate_spatial_dc(
            spec,
            "full_significant",
            significant_int_count,
            significant_variables,
            compute=compute,
        )
        executed_run_names.append("full_significant")

        if spec.search_mode == "full_significant":
            parsimonious_name = "full_significant"
        else:
            gts_iteration = 1
            criterion = selected_criterion_value
            exclude = list(insignificant)

            while True:
                run_prefix = f"exclude_{gts_iteration}"
                vars_int, vars_float, _ = eliminate_variable(
                    significant_variables,
                    tuple(exclude),
                    significant_binary,
                )
                executed_run_names.extend(
                    _leave_one_out_run_names(
                        vars_int + vars_float,
                        spec.protected_variables,
                        run_prefix,
                    )
                )
                selection, _ = leave_one_out(
                    vars_int + vars_float,
                    vars_int,
                    run_prefix,
                    repo,
                    spec.scenario,
                    str(spec.num_runs),
                    compute=compute,
                    protected_variables=set(spec.protected_variables),
                    output_label=spec.output_label,
                )

                best_variable_to_exclude = selection.loc[
                    selection[spec.selection_criterion].idxmin(), "variable"
                ]
                best_criterion = float(selection[spec.selection_criterion].min())
                if best_criterion >= criterion:
                    break

                exclude.append(best_variable_to_exclude)
                criterion = best_criterion
                gts_iteration += 1

            if gts_iteration == 1:
                parsimonious_name = "full_significant"
                parsimonious_coeffs, parsimonious_stats = postprocess_spatialdc(
                    repo,
                    spec.scenario,
                    parsimonious_name,
                    output_label=spec.output_label,
                )
            else:
                parsimonious_name = f"exclude_{gts_iteration - 1}_{exclude[-1].replace(':', '_')}"
                parsimonious_coeffs, parsimonious_stats = postprocess_spatialdc(
                    repo,
                    spec.scenario,
                    parsimonious_name,
                    output_label=spec.output_label,
                )
            selected_criterion_value = criterion

    _write_parsimonious_name(spec, parsimonious_name)

    if report:
        generated_reports.append(
            _write_estimates_report(
                spec,
                full_coefs,
                full_stats,
                parsimonious_coeffs,
                parsimonious_stats,
            ).name
        )
        generated_reports.append(
            _write_valuation_report(spec, parsimonious_name, parsimonious_coeffs).name
        )

    return EstimationResult(
        spec=spec,
        executed_run_names=tuple(executed_run_names),
        selected_run_name=parsimonious_name,
        selected_criterion_value=selected_criterion_value,
        removed_zero_variance=removed_zero_variance,
        removed_restricted_zero_variance=removed_restricted,
        surviving_binary_variables=surviving_binary,
        surviving_float_variables=surviving_float,
        generated_reports=tuple(generated_reports),
    )
