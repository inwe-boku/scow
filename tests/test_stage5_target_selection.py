import pandas as pd

from scow.postprocess_sites_helpers import annual_cost, build_cost_comparison, sites_to_target


def test_sites_to_target_includes_first_row_crossing_threshold():
    df = pd.DataFrame(
        {"CumSumEnergy": [1.2, 2.4, 3.7, 5.1]},
        index=[10, 20, 30, 40],
    )

    mask = sites_to_target(df, 3.0)

    assert mask.tolist() == [True, True, True, False]


def test_sites_to_target_returns_all_rows_when_target_not_reached():
    df = pd.DataFrame({"CumSumEnergy": [1.0, 1.5, 2.0]})

    mask = sites_to_target(df, 5.0)

    assert mask.tolist() == [True, True, True]


def test_annual_cost_adds_expected_columns():
    df = pd.DataFrame(
        {
            "soco": [100.0],
            "loco": [40.0],
            "lcoe": [60.0],
            "Energy": [2000.0],
        }
    )

    result = annual_cost(df)

    assert result.loc[0, "Annual Social Cost"] == 200.0
    assert result.loc[0, "Annual Local Social Cost"] == 80.0
    assert result.loc[0, "Annual Private Cost"] == 120.0


def test_build_cost_comparison_uses_target_masks_for_both_portfolios():
    private_df = pd.DataFrame(
        {
            "CumSumEnergy": [1.0, 3.0, 6.0],
            "Annual Private Cost": [10.0, 20.0, 30.0],
            "Annual Local Social Cost": [1.0, 2.0, 3.0],
            "Annual Social Cost": [11.0, 22.0, 33.0],
        }
    )
    social_df = pd.DataFrame(
        {
            "CumSumEnergy": [2.0, 5.5, 7.0],
            "Annual Private Cost": [12.0, 18.0, 24.0],
            "Annual Local Social Cost": [2.0, 1.0, 1.0],
            "Annual Social Cost": [14.0, 19.0, 25.0],
        }
    )

    result = build_cost_comparison(private_df, social_df, 5.0)

    assert result.loc["Private Cost", "Wind Turbines"] == 3
    assert result.loc["Private Cost", "Private Cost"] == 60.0
    assert result.loc["Social Cost", "Wind Turbines"] == 2
    assert result.loc["Social Cost", "Total Local Social Cost"] == 33.0
