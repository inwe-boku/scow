# %% imports
import logging
import cleo
import pandas as pd
import xarray as xr
import gams.transfer as gt
from scow.utils import sliced_location_optimization, locations_to_gdf
from config import repo, data_ver, spacing, gams_conf, nturbines

slice_count = 6
read_locations = False

# %%
atlas = cleo.Atlas(repo, "AUT", "EPSG:31287")
if atlas.region != "Niederösterreich":
    atlas.clip_to_nuts("Niederösterreich")
if not atlas.wind_turbines:
    atlas.wind_turbines = [
        "Enercon.E40.500", "Enercon.E82.3000", "Enercon.E101.3050", "Enercon.E115.3000",
        "Vestas.V100.1800", "Vestas.V100.2000", "Vestas.V112.3075"
    ]
if 'capacity_factors' not in atlas.wind.data.data_vars:
    atlas.wind.simulate_capacity_factors(bias_correction=0.71197)  # 0.71197
if 'lcoe' not in atlas.wind.data.data_vars:
    atlas.wind.compute_lcoe(turbine_cost_share=0.7)
if 'min_lcoe' not in atlas.wind.data.data_vars:
    atlas.wind.minimum_lcoe()
if 'optimal_power' not in atlas.wind.data.data_vars:
    atlas.wind.compute_optimal_power_energy()

energy = atlas.wind.data.optimal_energy.compute()
power = atlas.wind.data.optimal_power.compute()

social_cost = xr.open_dataset(repo / "data" / "results" / "soco.nc")

objectives = ["loco", "lcoe", "soco"]

for objective in objectives:
    if objective == "lcoe":
        cost_array = social_cost["min_lcoe"].compute()
        cost_threshold = 75
    elif objective == "loco":
        cost_array = social_cost["Local Social Cost"].compute()
        cost_threshold = 100
    elif objective == "soco":
        cost_array = (social_cost["Local Social Cost"] + social_cost["min_lcoe"]).compute()
        cost_threshold = 125

    if read_locations is False:
        gams_container = gt.Container()
        locations = sliced_location_optimization(gams_conf, gams_container, cost_array, num_slices=slice_count,
                                                 num_turbines=nturbines, space_px=spacing, min_turbines=1200,
                                                 max_turbines=6000, lcoe_thresh=cost_threshold,
                                                 gdx_out_string=objective, read_only=False)

        locations.to_csv(repo / "data" / "results" / f"opt_locations_{data_ver}_{objective}.csv")
    else:
        locations = pd.read_csv(repo / "data" / "results" / f"opt_locations_{data_ver}_{objective}.csv", index_col=[0])

    total_num_turbines = len(locations)

    # uncover total social cost at optimal locations
    cost_at_best_locations = [locations_to_gdf(cost_array, locations[(i * total_num_turbines):((i + 1) * total_num_turbines)],
                                               cost_name=objective, energy_array=energy, power_array=power) for i in range(0, slice_count)]
    cost_at_best_locations = pd.concat(cost_at_best_locations)
    cost_at_best_locations = cost_at_best_locations.dropna(how="any", axis=0)
    cost_at_best_locations = cost_at_best_locations.sort_values(by=objective)

    cost_at_best_locations["CumSumEnergy"] = cost_at_best_locations["Energy"].cumsum() / 1000
    cost_at_best_locations.to_csv(repo / "data" / "results" / f"opt_cost_{data_ver}_{objective}.csv")
    logging.info(f"Optimal cost for {objective} written to: "
                 f"{str(repo / 'data' / 'results' / f'opt_cost_{data_ver}_{objective}.csv')}")
