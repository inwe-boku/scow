import pandas as pd
import xarray as xr

from scow.postprocess_wtp_helpers import (
    build_interaction_layers,
    compute_social_cost_dataset,
    cost_distribution_frame,
)


def test_build_interaction_layers_adds_requested_products():
    landscape = xr.Dataset(
        {
            "a": xr.DataArray([[1.0, 2.0]]),
            "b": xr.DataArray([[3.0, 4.0]]),
        }
    )

    result = build_interaction_layers(landscape, ["a", "a:b"])

    assert "a:b" in result.data_vars
    assert result["a:b"].values.tolist() == [[3.0, 8.0]]


def test_compute_social_cost_dataset_excludes_min_lcoe_from_local_sum_and_normalizes():
    spatial = xr.Dataset(
        {
            "min_lcoe": xr.DataArray([[5.0, 6.0]]),
            "birds": xr.DataArray([[-2.0, 1.0]]),
            "roads": xr.DataArray([[1.0, 3.0]]),
        }
    )
    valuations = pd.Series({"min_lcoe": 1.0, "birds": 2.0, "roads": 1.0})

    soco = compute_social_cost_dataset(spatial, valuations)

    assert soco["min_lcoe"].values.tolist() == [[5.0, 6.0]]
    assert soco["Local Social Cost"].values.tolist() == [[0.0, 8.0]]
    assert soco["Local Social Cost"].attrs["unit"] == "€/MWh"


def test_cost_distribution_frame_labels_components():
    soco = xr.Dataset(
        {
            "min_lcoe": xr.DataArray([[1.0, 2.0], [3.0, 4.0]]),
            "Local Social Cost": xr.DataArray([[0.5, 1.5], [2.5, 3.5]]),
        }
    )

    result = cost_distribution_frame(soco, 1.0)

    assert set(result["Cost Component"]) == {
        "Local Social Cost",
        "Quasi-LCOE",
        "Total Local Social Cost",
    }
    assert "Cost [€/MWh]" in result.columns
