"""Pure helpers for stage-5 postprocessing."""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import wkt


def to_geopandas(df):
    """Convert a CSV-loaded stage-5 table with WKT geometry into a GeoDataFrame."""
    df = df.copy()
    df["geometry"] = df["geometry"].apply(wkt.loads)
    cols = [c for c in list(df.columns) if c not in ["geometry", "Unnamed: 0"]]
    return gpd.GeoDataFrame(geometry=df["geometry"], data=df[cols])


def sites_to_target(df, target, cumulative_col="CumSumEnergy"):
    """Return a positional mask including the first row that reaches the target."""
    reached = df[cumulative_col].to_numpy() >= target
    if reached.any():
        cutoff_pos = int(reached.argmax())
        mask = np.zeros(len(df), dtype=bool)
        mask[:cutoff_pos + 1] = True
        return mask
    return np.ones(len(df), dtype=bool)


def annual_cost(df):
    """Add annualized social, local-social, and private cost columns in million EUR."""
    df = df.copy()
    df["Annual Social Cost"] = df["soco"] * df["Energy"] / 1000
    df["Annual Local Social Cost"] = df["loco"] * df["Energy"] / 1000
    df["Annual Private Cost"] = df["lcoe"] * df["Energy"] / 1000
    return df


def build_cost_comparison(private_df, social_df, target):
    """Return the stage-5 annual cost comparison table up to the energy target."""
    cost_comparison = pd.DataFrame(
        index=["Private Cost", "Social Cost"],
        columns=["Wind Turbines", "Local Social Cost", "Private Cost", "Total Local Social Cost"],
    )

    private_sites = sites_to_target(private_df, target)
    cost_comparison.loc["Private Cost", "Wind Turbines"] = private_sites.sum()
    cost_comparison.loc["Private Cost", "Private Cost"] = private_df.loc[private_sites, "Annual Private Cost"].sum()
    cost_comparison.loc["Private Cost", "Local Social Cost"] = private_df.loc[
        private_sites, "Annual Local Social Cost"
    ].sum()
    cost_comparison.loc["Private Cost", "Total Local Social Cost"] = private_df.loc[
        private_sites, "Annual Social Cost"
    ].sum()

    social_sites = sites_to_target(social_df, target)
    cost_comparison.loc["Social Cost", "Wind Turbines"] = social_sites.sum()
    cost_comparison.loc["Social Cost", "Private Cost"] = social_df.loc[social_sites, "Annual Private Cost"].sum()
    cost_comparison.loc["Social Cost", "Local Social Cost"] = social_df.loc[
        social_sites, "Annual Local Social Cost"
    ].sum()
    cost_comparison.loc["Social Cost", "Total Local Social Cost"] = social_df.loc[
        social_sites, "Annual Social Cost"
    ].sum()
    return cost_comparison
