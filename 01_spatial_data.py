# %% imports
import tempfile
import time
from pathlib import Path
from contextlib import contextmanager

import cleo
import geopandas as gpd
import netCDF4
import numpy as np
import pandas as pd
import rioxarray as rxr
import xarray as xr
from fiona import drvsupport

from config import repo
from scow.site_data import (
    process_gip,
    process_protected_areas,
    process_water_bodies,
    process_wind_turbine_locations,
)
from scow.stages import AREA_NAME, ECONOMICS, LOSS_FACTOR, TURBINES_2014, TURBINES_MODERN

drvsupport.supported_drivers["LIBKML"] = "rw"

touch = True

ATLAS_ROOT = repo / "data" / "atlas"
SITE_ROOT = repo / "data" / "site"
PROCESSED_ROOT = repo / "data" / "processed"
PROCESSED_ROOT.mkdir(parents=True, exist_ok=True)


@contextmanager
def timed_step(name: str):
    start = time.perf_counter()
    print(f"[start] {name}", flush=True)
    try:
        yield
    finally:
        elapsed = time.perf_counter() - start
        print(f"[done] {name} ({elapsed:.1f}s)", flush=True)


def _write_temp_raster(da: xr.DataArray, tmpdir: Path, name: str) -> Path:
    raster = da.rename(name)
    if raster.rio.crs is None:
        raster = raster.rio.write_crs("EPSG:31287")
    path = tmpdir / f"{name}.tif"
    raster.rio.to_raster(path)
    return path


def _stage_temp_raster(atlas: cleo.Atlas, tmpdir: Path, da: xr.DataArray, name: str, if_exists: str = "replace") -> None:
    atlas.landscape.add(name, _write_temp_raster(da, tmpdir, name), if_exists=if_exists).materialize()


def _zero_on_valid_mask(valid_mask: xr.DataArray, name: str) -> xr.DataArray:
    data = xr.where(valid_mask, 0.0, np.nan).rename(name)
    data.attrs["units"] = "1"
    return data.rio.write_crs("EPSG:31287")


def _slope_degrees(elevation: xr.DataArray) -> xr.DataArray:
    elev = elevation.squeeze().load()
    data = elev.values.astype(np.float64)
    dy = float(abs(np.diff(elev["y"].values).mean()))
    dx = float(abs(np.diff(elev["x"].values).mean()))
    grad_y, grad_x = np.gradient(data, dy, dx)
    slope = np.degrees(np.arctan(np.sqrt(grad_x**2 + grad_y**2)))
    slope[~np.isfinite(data)] = np.nan
    return xr.DataArray(
        slope,
        coords=elev.coords,
        dims=elev.dims,
        name="slope",
        attrs={"units": "degree"},
    ).rio.write_crs("EPSG:31287")


def _prepare_tree_cover_density(atlas: cleo.Atlas) -> xr.DataArray:
    tree_cover_density = rxr.open_rasterio(
        SITE_ROOT / "forest" / "2015" / "TCD_2015_020m_at_31287_d02" / "TCD_2015_020m_at_31287_d02.tif",
        masked=True,
    ).squeeze()
    tree_cover_density = xr.where(tree_cover_density > 200, np.nan, tree_cover_density / 100)
    tree_cover_density.attrs["units"] = "%"
    template = atlas.landscape.data["elevation"]
    tree_cover_density = tree_cover_density.rio.write_crs("EPSG:31287").interp_like(template)
    return tree_cover_density.rename("tree_cover_density")


def _prepare_leaf_types(atlas: cleo.Atlas) -> tuple[xr.DataArray, xr.DataArray]:
    dominant_leaf_type = rxr.open_rasterio(
        SITE_ROOT / "forest" / "2015" / "DLT_2015_020m_at_31287_d02" / "DLT_2015_020m_at_31287_d02.tif",
        masked=True,
    ).squeeze()
    dominant_leaf_type = dominant_leaf_type.where(~(dominant_leaf_type == 255), np.nan)
    dominant_leaf_type = dominant_leaf_type.rio.write_crs("EPSG:31287").interp_like(atlas.landscape.data["elevation"])
    broadleaved = xr.where(dominant_leaf_type == 1, 1.0, 0.0).where(atlas.landscape.data["valid_mask"])
    coniferous = xr.where(dominant_leaf_type == 2, 1.0, 0.0).where(atlas.landscape.data["valid_mask"])
    broadleaved = broadleaved.rename("broadleaved").rio.write_crs("EPSG:31287")
    coniferous = coniferous.rename("coniferous").rio.write_crs("EPSG:31287")
    broadleaved.attrs["units"] = "1"
    coniferous.attrs["units"] = "1"
    return broadleaved, coniferous


def _compute_min_lcoe(atlas: cleo.Atlas, turbines: list[str]) -> xr.DataArray:
    with timed_step(f"select turbines ({len(turbines)})"):
        atlas.wind.select(turbines=turbines)
    try:
        with timed_step("materialize capacity_factors"):
            atlas.wind.compute("capacity_factors", loss_factor=LOSS_FACTOR).materialize(
                overwrite=True,
                allow_method_change=True,
            )
        with timed_step("materialize lcoe"):
            atlas.wind.compute(
                "lcoe",
                cf={"loss_factor": LOSS_FACTOR},
                economics=ECONOMICS,
            ).materialize(overwrite=True)
        with timed_step("materialize min_lcoe_turbine"):
            atlas.wind.compute(
                "min_lcoe_turbine",
                cf={"loss_factor": LOSS_FACTOR},
                economics=ECONOMICS,
            ).materialize(overwrite=True)
        with timed_step("materialize optimal_power"):
            atlas.wind.compute(
                "optimal_power",
                cf={"loss_factor": LOSS_FACTOR},
                economics=ECONOMICS,
            ).materialize(overwrite=True)
        with timed_step("materialize optimal_energy"):
            atlas.wind.compute(
                "optimal_energy",
                cf={"loss_factor": LOSS_FACTOR},
                economics=ECONOMICS,
            ).materialize(overwrite=True)
        with timed_step("derive min_lcoe from lcoe"):
            return atlas.wind.data["lcoe"].min(dim="turbine").rename("min_lcoe")
    finally:
        atlas.wind.clear_selection()


def _sanitize_netcdf_attrs(ds: xr.Dataset) -> xr.Dataset:
    def _coerce(value):
        if isinstance(value, (bool, np.bool_)):
            return int(value)
        if isinstance(value, (int, float, np.integer, np.floating, str)):
            return value
        return str(value)

    sanitized = ds.copy()
    sanitized.attrs = {key: _coerce(value) for key, value in sanitized.attrs.items()}
    for var in sanitized.data_vars:
        sanitized[var].attrs = {key: _coerce(value) for key, value in sanitized[var].attrs.items()}
    for coord in sanitized.coords:
        sanitized[coord].attrs = {key: _coerce(value) for key, value in sanitized[coord].attrs.items()}
    return sanitized


def _distance_unit_label(da: xr.DataArray):
    return da.attrs.get("units") or da.attrs.get("unit")


with timed_step("initialize atlas"):
    atlas = cleo.Atlas(ATLAS_ROOT, "AUT", "EPSG:31287")
    atlas.select(area=AREA_NAME, inplace=True)

# %% build missing landscape variables
with tempfile.TemporaryDirectory(prefix="scow_stage1_") as tmp:
    tmpdir = Path(tmp)
    landscape_vars = set(atlas.landscape.data.data_vars)

    # Known bootstrap exemptions: unavailable closed datasets become zero layers
    if "overnight_stays" not in landscape_vars:
        with timed_step("bootstrap exempt overnight_stays"):
            _stage_temp_raster(
                atlas,
                tmpdir,
                _zero_on_valid_mask(atlas.landscape.data["valid_mask"], "overnight_stays"),
                "overnight_stays",
            )
    if "important_bird_areas" not in landscape_vars:
        with timed_step("bootstrap exempt important_bird_areas"):
            _stage_temp_raster(
                atlas,
                tmpdir,
                _zero_on_valid_mask(atlas.landscape.data["valid_mask"], "important_bird_areas"),
                "important_bird_areas",
            )

    if "zoning" not in atlas.landscape.data.data_vars:
        with timed_step("rasterize zoning"):
            atlas.landscape.rasterize(
                SITE_ROOT / "zones" / "Zonierung_noe.shp",
                name="zoning",
                all_touched=touch,
                if_exists="noop",
            ).materialize()

    if "elevation" not in atlas.landscape.data.data_vars:
        with timed_step("add elevation"):
            atlas.landscape.add(
                "elevation",
                SITE_ROOT / "elevation" / "dhm_at_lamb_10m_2018.tif",
                if_exists="noop",
            ).materialize()

    if "slope" not in atlas.landscape.data.data_vars:
        with timed_step("compute slope"):
            _stage_temp_raster(atlas, tmpdir, _slope_degrees(atlas.landscape.data["elevation"]), "slope")

    if "water_bodies" not in atlas.landscape.data.data_vars:
        with timed_step("process water_bodies"):
            running_waters, standing_waters = process_water_bodies(SITE_ROOT / "water_bodies")
            atlas.landscape.rasterize(
                pd.concat([standing_waters, running_waters]),
                name="water_bodies",
                all_touched=touch,
                if_exists="noop",
            ).materialize()

    if "existing_turbines" not in atlas.landscape.data.data_vars:
        with timed_step("process existing_turbines"):
            turbines = process_wind_turbine_locations(atlas.landscape)
            atlas.landscape.rasterize(
                turbines.loc[turbines["Jahr"] < 2015, :],
                name="existing_turbines",
                all_touched=touch,
                if_exists="noop",
            ).materialize()

    if "preservation" not in atlas.landscape.data.data_vars:
        with timed_step("rasterize preservation"):
            atlas.landscape.rasterize(
                SITE_ROOT / "preservation" / "RRU_REGROP_ELTPolygon_0.shp",
                name="preservation",
                all_touched=touch,
                if_exists="noop",
            ).materialize()

    if "protected_areas" not in atlas.landscape.data.data_vars:
        with timed_step("process protected_areas"):
            protected_areas_path = tmpdir / "protected_areas.gpkg"
            process_protected_areas(repo, dissolve=False).to_file(protected_areas_path, driver="GPKG")
            atlas.landscape.rasterize(
                protected_areas_path,
                name="protected_areas",
                all_touched=touch,
                if_exists="replace",
            ).materialize()

    if "alps_convention" not in atlas.landscape.data.data_vars:
        with timed_step("rasterize alps_convention"):
            atlas.landscape.rasterize(
                SITE_ROOT / "alps" / "alpsconvention.shp",
                name="alps_convention",
                all_touched=touch,
                if_exists="noop",
            ).materialize()

    if "tree_cover_density" not in atlas.landscape.data.data_vars:
        with timed_step("process tree_cover_density"):
            _stage_temp_raster(atlas, tmpdir, _prepare_tree_cover_density(atlas), "tree_cover_density")

    if "broadleaved" not in atlas.landscape.data.data_vars or "coniferous" not in atlas.landscape.data.data_vars:
        with timed_step("process leaf types"):
            broadleaved, coniferous = _prepare_leaf_types(atlas)
            if "broadleaved" not in atlas.landscape.data.data_vars:
                _stage_temp_raster(atlas, tmpdir, broadleaved, "broadleaved")
            if "coniferous" not in atlas.landscape.data.data_vars:
                _stage_temp_raster(atlas, tmpdir, coniferous, "coniferous")

    if "pastures" not in atlas.landscape.data.data_vars:
        with timed_step("add pastures"):
            atlas.landscape.add_clc_category(231, name="pastures", if_exists="noop").materialize()

    if "airports_buff" not in atlas.landscape.data.data_vars:
        with timed_step("process airports_buff"):
            airports_clc = gpd.read_file(SITE_ROOT / "clc" / "CLC_2018_AT.shp")
            airports_clc["CODE_18"] = airports_clc["CODE_18"].astype(int)
            airports_clc = airports_clc.loc[airports_clc["CODE_18"] == 124].to_crs(atlas.crs)
            airports_clc = airports_clc.clip(atlas.get_nuts_area(AREA_NAME).geometry)

            airport_buffers = {
                "Vöslau": 2500,
                "Wr. Neustadt/West MIL": 5000,
                "Wr. Neustadt/Ost": 5000,
                "Wien-Schwechat": 6250,
                "Tulln MIL": 5000,
            }
            fname_most_recent_airspace = sorted((SITE_ROOT / "airspace").glob("*.kml"), reverse=True)[0]
            airports_ac = gpd.read_file(fname_most_recent_airspace, driver="LIBKML", layer="Airports").to_crs(atlas.crs)
            airports = airports_clc.sjoin(airports_ac, predicate="contains", how="left")
            for airport in list(airports["Name"]):
                airports.loc[airports["Name"] == airport, "geometry"] = airports[airports["Name"] == airport].buffer(
                    airport_buffers[airport]
                )
            atlas.landscape.rasterize(
                airports,
                name="airports_buff",
                all_touched=touch,
                if_exists="noop",
            ).materialize()

    if "restricted_military_areas" not in atlas.landscape.data.data_vars:
        with timed_step("process restricted_military_areas"):
            fname_most_recent_airspace = sorted((SITE_ROOT / "airspace").glob("*.kml"), reverse=True)[0]
            mil_areas = gpd.read_file(fname_most_recent_airspace, driver="LIBKML", layer="Restricted areas military (R)")
            atlas.landscape.rasterize(
                mil_areas,
                name="restricted_military_areas",
                all_touched=touch,
                if_exists="noop",
            ).materialize()

    if "power_lines" not in atlas.landscape.data.data_vars:
        with timed_step("process power_lines"):
            lines = gpd.read_file(SITE_ROOT / "buildings" / "BAU_2700_STROMLEITUNG_L_0.shp").to_crs(atlas.crs)
            lines = lines.clip(atlas.get_nuts_area(AREA_NAME))
            atlas.landscape.rasterize(lines, name="power_lines", all_touched=touch, if_exists="replace").materialize()

    if "roads" not in atlas.landscape.data.data_vars:
        with timed_step("process roads"):
            atlas.landscape.rasterize(
                process_gip(SITE_ROOT / "gip", dissolve=False),
                name="roads",
                all_touched=touch,
                if_exists="replace",
            ).materialize()

    if "buildings_in_greenland" not in atlas.landscape.data.data_vars:
        with timed_step("process buildings_in_greenland"):
            gruengeb = gpd.read_file(SITE_ROOT / "landuse" / "Erhaltenswerte_gebaeude.shp", encoding="utf-8").dissolve()
            atlas.landscape.rasterize(
                gruengeb,
                name="buildings_in_greenland",
                all_touched=True,
                if_exists="noop",
            ).materialize()

    if "residential_buildings" not in atlas.landscape.data.data_vars:
        with timed_step("process municipal zoning"):
            municipal_zoning = gpd.read_file(SITE_ROOT / "landuse" / "Widmungen_noe.shp", encoding="utf-8")

            residential_cats = [
                "Bauland Wohn- oder Mischnutzung",
                "Bauland Wohn- oder Mischnutzung nachhaltig",
            ]
            atlas.landscape.rasterize(
                municipal_zoning.loc[municipal_zoning["WIDMUNG"].str.contains("|".join(residential_cats)), :],
                name="residential_buildings",
                all_touched=touch,
                if_exists="replace",
            ).materialize()

            other_building_cats = [
                "Bauland Betriebsnutzung",
                "Bauland-Sondergebiet",
                "Bauland-Kerngebiet Handelseinrichtungen",
            ]
            atlas.landscape.rasterize(
                municipal_zoning.loc[municipal_zoning["WIDMUNG"].str.contains("|".join(other_building_cats)), :],
                name="other_building_land",
                all_touched=touch,
                if_exists="replace",
            ).materialize()

            greenland_zoning_cats = [
                "Grünland-Campingplatz",
                "Grünland-Kleingarten",
                "Grünland-land- und forstwirtschaftliche Hofstelle",
            ]
            atlas.landscape.rasterize(
                municipal_zoning.loc[municipal_zoning["WIDMUNG"].str.contains("|".join(greenland_zoning_cats)), :],
                name="greenland_zonings",
                all_touched=True,
                if_exists="replace",
            ).materialize()

    distance_sources = [
        "roads",
        "power_lines",
        "existing_turbines",
        "buildings_in_greenland",
        "residential_buildings",
        "other_building_land",
        "greenland_zonings",
    ]
    distance_names = [f"distance_{source}" for source in distance_sources]
    missing_distance_sources = [
        source
        for source, name in zip(distance_sources, distance_names, strict=True)
        if name not in atlas.landscape.data.data_vars
    ]
    if missing_distance_sources:
        with timed_step("compute distances"):
            atlas.landscape.compute("distance", source=missing_distance_sources, if_exists="noop").materialize()

    distance_names_needing_conversion = [
        name for name in distance_names if _distance_unit_label(atlas.landscape.data[name]) != "km"
    ]
    if distance_names_needing_conversion:
        with timed_step("convert distances to km"):
            for name in distance_names_needing_conversion:
                converted = atlas.landscape.convert_units(name, "km")
                _stage_temp_raster(atlas, tmpdir, converted, name)

    # %% save non-scenario-specific exports to disk
    with timed_step("export processed atlases"):
        _sanitize_netcdf_attrs(atlas.wind.data).to_netcdf(path=PROCESSED_ROOT / "WindAtlas_AUT_pro.nc")
        _sanitize_netcdf_attrs(atlas.landscape.data).to_netcdf(path=PROCESSED_ROOT / "LandscapeAtlas_AUT_pro.nc")

    # %% write data at time of decision to csv
    with timed_step("compute 2014 min_lcoe"):
        min_lcoe_2014 = _compute_min_lcoe(atlas, TURBINES_2014)
    with timed_step("export 2014 design matrix"):
        _stage_temp_raster(atlas, tmpdir, min_lcoe_2014, "min_lcoe")
        df = atlas.flatten(domain="landscape")
        df = df.dropna(how="any", axis=0)
        df.to_csv(PROCESSED_ROOT / "dc_data_2014.csv")

    with timed_step("compute modern min_lcoe"):
        min_lcoe_modern = _compute_min_lcoe(atlas, TURBINES_MODERN)
    with timed_step("export modern design matrix"):
        _stage_temp_raster(atlas, tmpdir, min_lcoe_modern, "min_lcoe")
        df = atlas.flatten(domain="landscape")
        df = df.dropna(how="any", axis=0)
        df.to_csv(PROCESSED_ROOT / "dc_data_modern.csv")
