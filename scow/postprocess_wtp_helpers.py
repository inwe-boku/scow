"""Pure helpers for stage-3 WTP and social-cost postprocessing."""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr

COST_UNIT_EUR_PER_MWH = "€/MWh"


def social_cost_display_name(name: str) -> str:
    """Return the canonical display label for one stage-3 cost layer."""

    if name == "min_lcoe":
        return "Quasi-LCOE"
    return name


def annotate_social_cost_dataset(soco: xr.Dataset) -> xr.Dataset:
    """Normalize stage-3 cost-layer metadata for plotting and NetCDF export."""

    normalized = soco.copy()
    for name, da in normalized.data_vars.items():
        attrs = {
            "long_name": social_cost_display_name(name),
            "units": COST_UNIT_EUR_PER_MWH,
            "unit": COST_UNIT_EUR_PER_MWH,
        }
        if name == "min_lcoe":
            attrs.update({key: value for key, value in da.attrs.items() if key.startswith("cleo:")})
        elif name != "Local Social Cost":
            attrs["source_variable"] = name
        da.attrs = attrs
    return normalized


def build_interaction_layers(landscape_data: xr.Dataset, terms) -> xr.Dataset:
    """Return a dataset with required interaction layers added for stage-3 terms."""
    enriched = landscape_data.copy()
    for term in terms:
        if ":" not in term:
            continue
        left, right = term.split(":")
        enriched[term] = enriched[left] * enriched[right]
    return enriched


def compute_social_cost_dataset(spatial_data: xr.Dataset, valuations: pd.Series) -> xr.Dataset:
    """Return component and local-social-cost rasters in €/MWh.

    `valuations` is expected to be indexed by the same variable names present in
    `spatial_data`. `"Local Social Cost"` includes every valued component except
    `min_lcoe`, which remains the private-cost anchor and is excluded from the
    local-social-cost sum.
    """
    soco = xr.Dataset()
    terms = list(valuations.index)
    local_terms = []

    for term in terms:
        soco[term] = spatial_data[term] * valuations.loc[term]
        if term != "min_lcoe":
            local_terms.append(soco[term])

    if local_terms:
        soco["Local Social Cost"] = sum(local_terms[1:], local_terms[0])
    else:
        template = soco["min_lcoe"] if "min_lcoe" in soco.data_vars else next(iter(soco.data_vars.values()))
        soco["Local Social Cost"] = xr.zeros_like(template)

    local_min = soco["Local Social Cost"].min()
    if float(local_min) < 0:
        soco["Local Social Cost"] = soco["Local Social Cost"] - local_min

    return annotate_social_cost_dataset(soco)


def cost_distribution_frame(soco: xr.Dataset, cut_off: float) -> pd.DataFrame:
    """Return the plotted cost-distribution frame for stage 3 in €/MWh."""

    def _flat_finite(da: xr.DataArray) -> np.ndarray:
        values = np.asarray(da).ravel()
        return values[np.isfinite(values)]

    loco_values = _flat_finite(soco["Local Social Cost"])
    cost_distribution = pd.DataFrame({"Cost [€/MWh]": loco_values})
    cost_distribution["target"] = 0
    cost_distribution["Cost Component"] = "Local Social Cost"

    lcoe_values = _flat_finite(soco["min_lcoe"])
    lcoe_frame = pd.DataFrame({"Cost [€/MWh]": lcoe_values})
    lcoe_frame = lcoe_frame.loc[
        lcoe_frame["Cost [€/MWh]"] <= lcoe_frame["Cost [€/MWh]"].quantile(cut_off)
    ].copy()
    lcoe_frame["target"] = 1
    lcoe_frame["Cost Component"] = "Quasi-LCOE"
    cost_distribution = pd.concat([cost_distribution, lcoe_frame], ignore_index=True)

    total_values = _flat_finite(soco["Local Social Cost"] + soco["min_lcoe"])
    total_frame = pd.DataFrame({"Cost [€/MWh]": total_values})
    total_frame = total_frame.loc[
        total_frame["Cost [€/MWh]"] <= total_frame["Cost [€/MWh]"].quantile(cut_off)
    ].copy()
    total_frame["target"] = 2
    total_frame["Cost Component"] = "Total Local Social Cost"
    return pd.concat([cost_distribution, total_frame], ignore_index=True)
