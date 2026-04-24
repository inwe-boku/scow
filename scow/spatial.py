from __future__ import annotations

import time
from contextlib import contextmanager
from pathlib import Path

import cleo
import geopandas as gpd
import numpy as np
import pandas as pd
import rioxarray as rxr
import xarray as xr
from fiona import drvsupport

from scow.site_data import (
    process_gip,
    process_protected_areas,
    process_water_bodies,
    process_wind_turbine_locations,
)
from scow.study import ECONOMICS, LOSS_FACTOR, ScowStudy, Stage1DesignMatrixSpec

drvsupport.supported_drivers["LIBKML"] = "rw"

TOUCH = True


@contextmanager
def timed_step(name: str):
    """Print a coarse stage-1 timing line for one build step."""
    start = time.perf_counter()
    print(f"[start] {name}", flush=True)
    try:
        yield
    finally:
        elapsed = time.perf_counter() - start
        print(f"[done] {name} ({elapsed:.1f}s)", flush=True)


def _site_root(repo_root: Path) -> Path:
    return repo_root / "data" / "site"


def _atlas_root(repo_root: Path) -> Path:
    return repo_root / "data" / "atlas"


def _processed_root(repo_root: Path) -> Path:
    path = repo_root / "data" / "processed"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _stage_computed_raster(
    atlas: cleo.Atlas,
    da: xr.DataArray,
    name: str,
    if_exists: str = "replace",
) -> None:
    raster = da.rename(name)
    if raster.rio.crs is None:
        raster = raster.rio.write_crs(atlas.crs)
    atlas.landscape.add_dataarray(name, raster, if_exists=if_exists).materialize()


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


def _prepare_tree_cover_density(atlas: cleo.Atlas, repo_root: Path) -> xr.DataArray:
    site_root = _site_root(repo_root)
    tree_cover_density = rxr.open_rasterio(
        site_root / "forest" / "2015" / "TCD_2015_020m_at_31287_d02" / "TCD_2015_020m_at_31287_d02.tif",
        masked=True,
    ).squeeze()
    tree_cover_density = xr.where(tree_cover_density > 200, np.nan, tree_cover_density / 100)
    tree_cover_density.attrs["units"] = "%"
    template = atlas.landscape.data["elevation"]
    tree_cover_density = tree_cover_density.rio.write_crs("EPSG:31287").interp_like(template)
    return tree_cover_density.rename("tree_cover_density")


def _prepare_leaf_types(atlas: cleo.Atlas, repo_root: Path) -> tuple[xr.DataArray, xr.DataArray]:
    site_root = _site_root(repo_root)
    dominant_leaf_type = rxr.open_rasterio(
        site_root / "forest" / "2015" / "DLT_2015_020m_at_31287_d02" / "DLT_2015_020m_at_31287_d02.tif",
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


def _compute_min_lcoe(
    atlas: cleo.Atlas,
    turbines: tuple[str, ...] | list[str],
    *,
    loss_factor: float,
    economics: dict[str, float],
) -> xr.DataArray:
    with timed_step(f"select turbines ({len(turbines)})"):
        atlas.wind.select(turbines=list(turbines))
    try:
        with timed_step("materialize capacity_factors"):
            atlas.wind.compute("capacity_factors", loss_factor=loss_factor).materialize(
                overwrite=True,
                allow_method_change=True,
            )
        with timed_step("materialize lcoe"):
            atlas.wind.compute(
                "lcoe",
                cf={"loss_factor": loss_factor},
                economics=economics,
            ).materialize(overwrite=True)
        with timed_step("materialize min_lcoe_turbine"):
            atlas.wind.compute(
                "min_lcoe_turbine",
                cf={"loss_factor": loss_factor},
                economics=economics,
            ).materialize(overwrite=True)
        with timed_step("materialize optimal_power"):
            atlas.wind.compute(
                "optimal_power",
                cf={"loss_factor": loss_factor},
                economics=economics,
            ).materialize(overwrite=True)
        with timed_step("materialize optimal_energy"):
            atlas.wind.compute(
                "optimal_energy",
                cf={"loss_factor": loss_factor},
                economics=economics,
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


def _load_vector_source(source: Path, *, target_crs: str, encoding: str | None = None) -> gpd.GeoDataFrame:
    """Load one vector source path and reproject it to the target CRS."""

    frame = gpd.read_file(source, encoding=encoding)
    if frame.crs is None:
        raise ValueError(f"Vector source has no CRS: {source}")
    return frame.to_crs(target_crs)


def build_spatial_atlas(study: ScowStudy, repo_root: Path) -> cleo.Atlas:
    """Initialize and populate the SCOW atlas stores used by stage 1."""
    with timed_step("initialize atlas"):
        atlas = cleo.Atlas(_atlas_root(repo_root), "AUT", "EPSG:31287")
        atlas.select(area=study.area_name, inplace=True)

    ensure_scow_landscape_layers(atlas, study, repo_root)
    ensure_scow_distance_layers(atlas)
    return atlas


def ensure_scow_landscape_layers(atlas: cleo.Atlas, study: ScowStudy, repo_root: Path) -> None:
    """Ensure the current SCOW landscape recipe is present in the atlas."""
    site_root = _site_root(repo_root)
    payload = study.stage1_payload(repo_root)
    landscape_vars = set(atlas.landscape.data.data_vars)

    if "overnight_stays" not in landscape_vars:
        with timed_step("bootstrap exempt overnight_stays"):
            _stage_computed_raster(
                atlas,
                _zero_on_valid_mask(atlas.landscape.data["valid_mask"], "overnight_stays"),
                "overnight_stays",
            )
    if "important_bird_areas" not in landscape_vars:
        with timed_step("bootstrap exempt important_bird_areas"):
            _stage_computed_raster(
                atlas,
                _zero_on_valid_mask(atlas.landscape.data["valid_mask"], "important_bird_areas"),
                "important_bird_areas",
            )

    if "zoning" not in atlas.landscape.data.data_vars:
        with timed_step("rasterize zoning"):
            atlas.landscape.rasterize(
                _load_vector_source(payload.zoning_source, target_crs=atlas.crs),
                name="zoning",
                all_touched=TOUCH,
                if_exists="noop",
            ).materialize()

    if "elevation" not in atlas.landscape.data.data_vars:
        with timed_step("add elevation"):
            atlas.landscape.add(
                "elevation",
                site_root / "elevation" / "dhm_at_lamb_10m_2018.tif",
                if_exists="noop",
            ).materialize()

    if "slope" not in atlas.landscape.data.data_vars:
        with timed_step("compute slope"):
            _stage_computed_raster(atlas, _slope_degrees(atlas.landscape.data["elevation"]), "slope")

    if "water_bodies" not in atlas.landscape.data.data_vars:
        with timed_step("process water_bodies"):
            running_waters, standing_waters = process_water_bodies(site_root / "water_bodies")
            atlas.landscape.rasterize(
                pd.concat([standing_waters, running_waters]),
                name="water_bodies",
                all_touched=TOUCH,
                if_exists="noop",
            ).materialize()

    if "existing_turbines" not in atlas.landscape.data.data_vars:
        with timed_step("process existing_turbines"):
            turbines = process_wind_turbine_locations(atlas.landscape)
            atlas.landscape.rasterize(
                turbines.loc[turbines["Jahr"] < payload.existing_turbines_built_before_year, :],
                name="existing_turbines",
                all_touched=TOUCH,
                if_exists="noop",
            ).materialize()

    if "preservation" not in atlas.landscape.data.data_vars:
        with timed_step("rasterize preservation"):
            atlas.landscape.rasterize(
                site_root / "preservation" / "RRU_REGROP_ELTPolygon_0.shp",
                name="preservation",
                all_touched=TOUCH,
                if_exists="noop",
            ).materialize()

    if "protected_areas" not in atlas.landscape.data.data_vars:
        with timed_step("process protected_areas"):
            atlas.landscape.rasterize(
                process_protected_areas(repo_root, dissolve=False),
                name="protected_areas",
                all_touched=TOUCH,
                if_exists="replace",
            ).materialize()

    if "alps_convention" not in atlas.landscape.data.data_vars:
        with timed_step("rasterize alps_convention"):
            atlas.landscape.rasterize(
                site_root / "alps" / "alpsconvention.shp",
                name="alps_convention",
                all_touched=TOUCH,
                if_exists="noop",
            ).materialize()

    if "tree_cover_density" not in atlas.landscape.data.data_vars:
        with timed_step("process tree_cover_density"):
            _stage_computed_raster(atlas, _prepare_tree_cover_density(atlas, repo_root), "tree_cover_density")

    if "broadleaved" not in atlas.landscape.data.data_vars or "coniferous" not in atlas.landscape.data.data_vars:
        with timed_step("process leaf types"):
            broadleaved, coniferous = _prepare_leaf_types(atlas, repo_root)
            if "broadleaved" not in atlas.landscape.data.data_vars:
                _stage_computed_raster(atlas, broadleaved, "broadleaved")
            if "coniferous" not in atlas.landscape.data.data_vars:
                _stage_computed_raster(atlas, coniferous, "coniferous")

    if "pastures" not in atlas.landscape.data.data_vars:
        with timed_step("add pastures"):
            atlas.landscape.add_clc_category(231, name="pastures", if_exists="noop").materialize()

    if "airports_buff" not in atlas.landscape.data.data_vars:
        with timed_step("process airports_buff"):
            airports_clc = gpd.read_file(site_root / "clc" / "CLC_2018_AT.shp")
            airports_clc["CODE_18"] = airports_clc["CODE_18"].astype(int)
            airports_clc = airports_clc.loc[airports_clc["CODE_18"] == 124].to_crs(atlas.crs)
            airports_clc = airports_clc.clip(atlas.get_nuts_area(study.area_name).geometry)

            airport_buffers = {
                "Vöslau": 2500,
                "Wr. Neustadt/West MIL": 5000,
                "Wr. Neustadt/Ost": 5000,
                "Wien-Schwechat": 6250,
                "Tulln MIL": 5000,
            }
            fname_most_recent_airspace = sorted((site_root / "airspace").glob("*.kml"), reverse=True)[0]
            airports_ac = gpd.read_file(fname_most_recent_airspace, driver="LIBKML", layer="Airports").to_crs(atlas.crs)
            airports = airports_clc.sjoin(airports_ac, predicate="contains", how="left")
            for airport in list(airports["Name"]):
                airports.loc[airports["Name"] == airport, "geometry"] = airports[airports["Name"] == airport].buffer(
                    airport_buffers[airport]
                )
            atlas.landscape.rasterize(
                airports,
                name="airports_buff",
                all_touched=TOUCH,
                if_exists="noop",
            ).materialize()

    if "restricted_military_areas" not in atlas.landscape.data.data_vars:
        with timed_step("process restricted_military_areas"):
            fname_most_recent_airspace = sorted((site_root / "airspace").glob("*.kml"), reverse=True)[0]
            mil_areas = gpd.read_file(fname_most_recent_airspace, driver="LIBKML", layer="Restricted areas military (R)")
            atlas.landscape.rasterize(
                mil_areas,
                name="restricted_military_areas",
                all_touched=TOUCH,
                if_exists="noop",
            ).materialize()

    if "power_lines" not in atlas.landscape.data.data_vars:
        with timed_step("process power_lines"):
            lines = gpd.read_file(site_root / "buildings" / "BAU_2700_STROMLEITUNG_L_0.shp").to_crs(atlas.crs)
            lines = lines.clip(atlas.get_nuts_area(study.area_name))
            atlas.landscape.rasterize(lines, name="power_lines", all_touched=TOUCH, if_exists="replace").materialize()

    if "roads" not in atlas.landscape.data.data_vars:
        with timed_step("process roads"):
            atlas.landscape.rasterize(
                process_gip(site_root / "gip", dissolve=False),
                name="roads",
                all_touched=TOUCH,
                if_exists="replace",
            ).materialize()

    if "buildings_in_greenland" not in atlas.landscape.data.data_vars:
        with timed_step("process buildings_in_greenland"):
            gruengeb = _load_vector_source(
                payload.buildings_in_greenland_source,
                target_crs=atlas.crs,
                encoding="utf-8",
            ).dissolve()
            atlas.landscape.rasterize(
                gruengeb,
                name="buildings_in_greenland",
                all_touched=True,
                if_exists="noop",
            ).materialize()

    if "residential_buildings" not in atlas.landscape.data.data_vars:
        with timed_step("process municipal zoning"):
            municipal_zoning = _load_vector_source(
                payload.municipal_zoning_source,
                target_crs=atlas.crs,
                encoding="utf-8",
            )

            residential_cats = [
                "Bauland Wohn- oder Mischnutzung",
                "Bauland Wohn- oder Mischnutzung nachhaltig",
            ]
            atlas.landscape.rasterize(
                municipal_zoning.loc[municipal_zoning["WIDMUNG"].str.contains("|".join(residential_cats)), :],
                name="residential_buildings",
                all_touched=TOUCH,
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
                all_touched=TOUCH,
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


def ensure_scow_distance_layers(atlas: cleo.Atlas) -> None:
    """Ensure the SCOW distance layers exist and are stored in kilometers."""
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
                _stage_computed_raster(atlas, converted, name)


def export_processed_atlases(atlas: cleo.Atlas, processed_root: Path) -> None:
    """Write the current processed atlas exports to the canonical NetCDF paths."""
    with timed_step("export processed atlases"):
        _sanitize_netcdf_attrs(atlas.wind.data).to_netcdf(path=processed_root / "WindAtlas_AUT_pro.nc")
        _sanitize_netcdf_attrs(atlas.landscape.data).to_netcdf(path=processed_root / "LandscapeAtlas_AUT_pro.nc")


def export_design_matrix_for_turbines(
    atlas: cleo.Atlas,
    spec: Stage1DesignMatrixSpec,
    processed_root: Path,
) -> None:
    """Compute scenario-specific min_lcoe and write the flattened design matrix."""
    with timed_step(f"compute {spec.name} min_lcoe"):
        min_lcoe = _compute_min_lcoe(
            atlas,
            spec.turbines,
            loss_factor=LOSS_FACTOR,
            economics=ECONOMICS,
        )
    with timed_step(f"export {spec.name} design matrix"):
        _stage_computed_raster(atlas, min_lcoe, "min_lcoe")
        design_matrix = atlas.flatten(domain="landscape")
        design_matrix = design_matrix.dropna(how="any", axis=0)
        design_matrix.to_csv(processed_root / f"dc_data_{spec.name}.csv")


def prepare_spatial_data(study: ScowStudy, repo_root: Path) -> None:
    """Run the current SCOW stage-1 build and export workflow."""
    processed_root = _processed_root(repo_root)
    atlas = build_spatial_atlas(study, repo_root)
    export_processed_atlases(atlas, processed_root)
    payload = study.stage1_payload(repo_root)
    for spec in payload.design_matrix_specs:
        export_design_matrix_for_turbines(atlas, spec, processed_root)
