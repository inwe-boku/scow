import pytest
import xarray as xr

from scow.artifact_validators import validate_soco_dataset


def test_validate_soco_dataset_accepts_current_contract():
    dataset = xr.Dataset(
        {
            "min_lcoe": xr.DataArray([[1.0]]),
            "Local Social Cost": xr.DataArray([[2.0]], attrs={"unit": "€/MWh"}),
        }
    )

    validated = validate_soco_dataset(dataset)

    assert "min_lcoe" in validated


def test_validate_soco_dataset_rejects_missing_local_social_cost_unit():
    dataset = xr.Dataset(
        {
            "min_lcoe": xr.DataArray([[1.0]]),
            "Local Social Cost": xr.DataArray([[2.0]]),
        }
    )

    with pytest.raises(ValueError, match="unit '€/MWh'"):
        validate_soco_dataset(dataset)
