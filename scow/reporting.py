from __future__ import annotations

import logging
from pathlib import Path

import contextily as cx
import geopandas as gpd
import matplotlib.pyplot as plt
import seaborn as sns

from scow.artifact_validators import validate_stage4_costs
from scow.logging_utils import configure_script_logging
from scow.paths import (
    cost_curve_figure_path,
    dist_soco_figure_path,
    latex_cost_comparison_path,
    map_lcoe_figure_path,
    map_soco_figure_path,
    map_turbines_figure_path,
    soco_facet_figure_path,
)
from scow.postprocess_sites_helpers import annual_cost, build_cost_comparison, sites_to_target
from scow.postprocess_wtp_helpers import cost_distribution_frame
from scow.social_cost import get_social_cost_atlas, open_social_cost_atlas
from scow.stage4_outputs import load_stage4_costs


def latex_cost_comparison_fragment(df):
    """Render the current LaTeX body-row fragment without relying on pandas line offsets."""

    lines = []
    for label, row in df.iterrows():
        cells = [str(int(row["Wind Turbines"]))]
        cells.extend(f"{float(row[col]):.1f}" for col in df.columns if col != "Wind Turbines")
        lines.append(f"{label} & " + " & ".join(cells) + r" \\")
    return "\n".join(lines)


def generate_reports(
    *,
    repo_root: Path,
    scenario: str,
    area_name: str,
    generation_target: float,
    logger: logging.Logger | logging.LoggerAdapter | None = None,
) -> None:
    """Generate the current unified stage-5 reporting outputs."""

    if logger is None:
        logger = configure_script_logging("stage5")

    logger.info("Loading stage-5 reporting inputs")
    state = gpd.read_file(repo_root / "data" / "site" / "vgd" / "vgd_oesterreich.shp")
    state = state[["BL", "geometry"]].dissolve(by="BL").reset_index()
    state = state[state.BL == area_name]

    social_cost_atlas = get_social_cost_atlas(scenario)
    social_cost = open_social_cost_atlas(social_cost_atlas, engine="netcdf4")
    opt_cost_soco = validate_stage4_costs(load_stage4_costs(scenario, "soco"), "soco")
    opt_cost_loco = validate_stage4_costs(load_stage4_costs(scenario, "loco"), "loco")
    opt_cost_lcoe = validate_stage4_costs(load_stage4_costs(scenario, "lcoe"), "lcoe")

    lcoe = social_cost["min_lcoe"].copy()
    lcoe.name = "Quasi-LCOE"
    lcoe.attrs["unit"] = "€/MWh"

    fig, ax = plt.subplots(figsize=(7, 5))
    social_cost["Local Social Cost"].plot(ax=ax, robust=True)
    ax.set_title("")
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(map_soco_figure_path(scenario), dpi=300)
    plt.close()

    fig, ax = plt.subplots(figsize=(7, 5))
    lcoe.plot(ax=ax, robust=True)
    ax.set_title("")
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(map_lcoe_figure_path(scenario), dpi=300)
    plt.close()
    logger.info("Maps of Local Social Cost and LCOE plotted")

    cost_by_turbine = cost_distribution_frame(social_cost, 0.999)
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.kdeplot(data=cost_by_turbine, x="Cost [€/MWh]", hue="Cost Component", fill=True)
    plt.grid(alpha=0.4)
    ax.set_axisbelow(True)
    plt.tight_layout()
    plt.savefig(dist_soco_figure_path(scenario), dpi=300)
    plt.close()
    logger.info("Cost distributions plotted")

    num_vars = len(social_cost.data_vars)
    num_cols = 4
    num_rows = (num_vars + num_cols - 1) // num_cols
    fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(26, num_rows * 5))
    axes = axes.flatten()
    for i, var in enumerate(social_cost.data_vars):
        social_cost[var].plot(ax=axes[i])
        axes[i].set_title(var)
        axes[i].tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    for j in range(i + 1, num_rows * num_cols):
        fig.delaxes(axes[j])
    plt.tight_layout()
    plt.savefig(soco_facet_figure_path(scenario), dpi=300)
    plt.close()
    logger.info("Facet of social cost variables plotted")

    opt_cost_soco["lcoe"] = social_cost["min_lcoe"].sel(
        x=opt_cost_soco.geometry.x.to_xarray(),
        y=opt_cost_soco.geometry.y.to_xarray(),
        method="nearest",
    ).to_pandas()
    opt_cost_soco["loco"] = opt_cost_soco["soco"] - opt_cost_soco["lcoe"]
    opt_cost_soco = annual_cost(opt_cost_soco)

    opt_cost_lcoe["loco"] = social_cost["Local Social Cost"].sel(
        x=opt_cost_lcoe.geometry.x.to_xarray(),
        y=opt_cost_lcoe.geometry.y.to_xarray(),
        method="nearest",
    ).to_pandas()
    opt_cost_lcoe["soco"] = opt_cost_lcoe["lcoe"] + opt_cost_lcoe["loco"]
    opt_cost_lcoe = annual_cost(opt_cost_lcoe)

    opt_cost_loco["lcoe"] = social_cost["min_lcoe"].sel(
        x=opt_cost_loco.geometry.x.to_xarray(),
        y=opt_cost_loco.geometry.y.to_xarray(),
        method="nearest",
    ).to_pandas()
    opt_cost_loco["soco"] = opt_cost_loco["lcoe"] + opt_cost_loco["loco"]
    opt_cost_loco = annual_cost(opt_cost_loco)

    cost_comparison = build_cost_comparison(opt_cost_lcoe, opt_cost_soco, generation_target)
    latex_cost_comparison_path().write_text(latex_cost_comparison_fragment(cost_comparison))
    logger.info("Wrote cost comparison table to %s", latex_cost_comparison_path())

    for objective, df in zip(
        ["lcoe", "loco", "soco"],
        [opt_cost_lcoe, opt_cost_loco, opt_cost_soco],
    ):
        fig, ax = plt.subplots(figsize=(8, 5))
        state.plot(ax=ax, alpha=0.3, edgecolor="k")
        df.loc[sites_to_target(df, generation_target)].plot(ax=ax, color="red", markersize=2)
        try:
            cx.add_basemap(ax, crs=state.crs, source=cx.providers.BasemapAT.terrain, zoom=10)
        except Exception as exc:
            logger.warning("Basemap fetch failed for %s: %s", objective, exc)
        ax.set_axis_off()
        plt.tight_layout()
        plt.savefig(map_turbines_figure_path(objective), dpi=300)
        plt.close()

    logger.info("Maps of optimal wind turbine sites saved")

    fig, ax = plt.subplots(figsize=(8, 5))
    plt.stackplot(
        opt_cost_soco["CumSumEnergy"],
        opt_cost_soco["lcoe"],
        opt_cost_soco["loco"],
        alpha=0.825,
    )
    ax.set_xlabel("Energy [TWh]")
    ax.set_ylabel("Social Cost [€/MWh]")
    ax.set_xlim([0, 60])
    plt.grid()
    ax.set_axisbelow(True)
    plt.legend(["Quasi-LCOE", "Local Social Cost"], loc="upper left")
    plt.tight_layout()
    plt.savefig(cost_curve_figure_path(), dpi=300)
    plt.close()

    logger.info("Social cost curve saved")
