# %% imports
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import rioxarray as rxr
import netCDF4
import cleo
from fiona import drvsupport
from xrspatial import slope as xrs_slope
from config import repo
from scow.site_data import (process_tourism, process_water_bodies, process_gip, process_protected_areas,
                            generate_data_dict, process_wind_turbine_locations)

drvsupport.supported_drivers['LIBKML'] = 'rw'

touch = True

# %% initialize Atlas for Austria
atlas = cleo.Atlas(str(repo), "AUT", "epsg:31287")
atlas.clip_to_nuts("Niederösterreich", inplace=True)
if not atlas.wind_turbines:
    atlas.wind_turbines = [
        "Enercon.E40.500", "Enercon.E82.3000", "Enercon.E101.3050", "Enercon.E115.3000",
        "Vestas.V100.1800", "Vestas.V100.2000", "Vestas.V112.3075"]  # wind turbines as of 2014

atlas.wind.compute_wind_shear_coefficient()
atlas.wind.compute_air_density_correction()
atlas.wind.compute_weibull_pdf()

# %% download site data sets
dict = generate_data_dict(atlas.path / "data" / "site", wdpa_date="Dec2024")
atlas.landscape.load_and_extract_from_dict(dict)

# # # # # # PROPRIETARY DATA # # # # #
# %% tourism
# request "Ankünfte und Nächtigungen nach Gemeinden, 2018" from Statistik Austria: info@statistik.gv.at
stays = process_tourism(atlas.path / 'data' / 'site')
stays["overnight_stays"] = stays["overnight_stays"] / 1000
atlas.landscape.rasterize(stays, column='overnight_stays', name="overnight_stays", all_touched=touch)
del stays

# %% important bird areas
# request important bird areas GIS data via: https://datazone.birdlife.org/site/requestgis
iba = gpd.read_file(atlas.path / 'data' / 'site' / 'iba' / 'EUROPE_KBA.shp')
# BUG: potential incompatibility of GDAL, geopandas and fiona may result in fiona.ogrext.FeatureBuilder.build ValueError
# however, error does not affect spatial data
iba = iba.loc[iba['ISO3'] == 'AUT', :]
atlas.landscape.rasterize(iba, name='important_bird_areas', all_touched=touch)
del iba

# # # # # OPEN DATA # # # # #
# %% wind power zoning
if "zoning" not in atlas.landscape.data.data_vars:
    atlas.landscape.rasterize(str(atlas.path / 'data' / 'site' / 'zones' / 'Zonierung_noe.shp'), name='zoning',
                              all_touched=touch)

# %% elevation and slope
if "elevation" not in atlas.landscape.data.data_vars:
    elevation = rxr.open_rasterio(atlas.path / 'data' / 'site' / 'elevation' / 'dhm_at_lamb_10m_2018.tif', masked=True)
    elevation = elevation.squeeze()
    # clip elevation to nuts region
    clips = atlas.get_nuts_region(atlas.region)
    # clips = clips.to_crs(atlas.crs)
    elevation = elevation.rio.reproject(atlas.crs)
    elevation = elevation.rio.clip(clips.geometry)
    # compute slope
    slope = xrs_slope(elevation)
    elevation = elevation.interp_like(atlas.landscape.data["template"])
    slope = slope.interp_like(atlas.landscape.data["template"])
    # add data
    atlas.landscape.add(elevation, name='elevation')
    atlas.landscape.add(slope, name='slope')
    del elevation, slope

# %% water bodies
if "water_bodies" not in atlas.landscape.data.data_vars:
    running_waters, standing_waters = process_water_bodies(atlas.path / 'data' / 'site' / 'water_bodies')
    water_bodies = pd.concat([standing_waters, running_waters])
    atlas.landscape.rasterize(water_bodies, name='water_bodies', all_touched=touch)
    del running_waters, standing_waters, water_bodies

# %% pre-existing wind turbine locations
if "existing_turbines" not in atlas.landscape.data.data_vars:
    turbines = process_wind_turbine_locations(atlas.landscape)
    atlas.landscape.rasterize(turbines.loc[turbines['Jahr'] < 2015, :], name='existing_turbines', all_touched=touch)

# %% landscapes worthy of preservation - "Erhaltenswerte Landschaftsteile Niederösterreich"
if "preservation" not in atlas.landscape.data.data_vars:
    atlas.landscape.rasterize(str(atlas.path / 'data' / 'site' / 'preservation' / 'RRU_REGROP_ELTPolygon_0.shp'),
                              name='preservation', all_touched=touch)

# %% protected areas
if "protected_areas" not in atlas.landscape.data.data_vars:
    protected_areas = process_protected_areas(atlas.path)
    atlas.landscape.rasterize(protected_areas, name='protected_areas', all_touched=touch)
    del protected_areas

# %% alpine convention
if "alps_conventions" not in atlas.landscape.data.data_vars:
    atlas.landscape.rasterize(str(atlas.path / 'data' / 'site' / 'alps' / 'alpsconvention.shp'), name='alps_convention')

# %% tree cover density
if "tree_cover_density" not in atlas.landscape.data.data_vars:
    tree_cover_density = rxr.open_rasterio(
        atlas.path / 'data' / 'site' / 'forest' / '2015' / 'TCD_2015_020m_at_31287_d02/TCD_2015_020m_at_31287_d02.tif').squeeze()
    tree_cover_density = xr.where(tree_cover_density > 200, np.nan, tree_cover_density / 100).rio.write_crs(
        "epsg:31287")
    tree_cover_density.attrs['unit'] = '%'
    atlas.landscape.add(tree_cover_density, name='tree_cover_density')
    del tree_cover_density

# %% broadleaved and coniferous forests
# leaf type: 0 - not tree covered, 1: broadleaved trees, 2: coniferous trees, >2: unclassifiable
if "broadleaved" not in atlas.landscape.data.data_vars:
    dominant_leaf_type = rxr.open_rasterio(
        atlas.path / 'data' / 'site' / 'forest' / '2015' / 'DLT_2015_020m_at_31287_d02/DLT_2015_020m_at_31287_d02.tif').squeeze()
    dominant_leaf_type = dominant_leaf_type.where(~(dominant_leaf_type == 255), np.nan)
    dominant_leaf_type = dominant_leaf_type.rio.clip(atlas.get_nuts_region(atlas.region).geometry).interp_like(
        atlas.landscape.data.template)

    dominant_broadleaved = xr.where(dominant_leaf_type == 1, 1, atlas.landscape.data.template).rio.write_crs(
        "epsg:31287")
    atlas.landscape.add(dominant_broadleaved, name='broadleaved')

    dominant_coniferous = xr.where(dominant_leaf_type == 2, 1, atlas.landscape.data.template).rio.write_crs(
        "epsg:31287")
    atlas.landscape.add(dominant_coniferous, name='coniferous')

    del dominant_leaf_type, dominant_broadleaved, dominant_coniferous

# %% Pastures
# clc: pastures (231) & crop areas (211, 212, 213, 221, 222, 223, 241, 242, 243)
# 211, 212, 213: arable land
# 221, 222, 223: permanent crops
# 241, 242, 243: heterogeneous agricultural areas -- check whether 242 and 243 should be excluded
# TODO: allow aggregates of classes like crop areas:
#  crops = clc[clc['CODE_18'].isin([211, 212, 213, 221, 222, 223, 241, 242, 243])]

atlas.landscape.add_corine_land_cover(231)

# %% airports
# airport area as GeoDataFrame from corine land cover
if "airport_buff" not in atlas.landscape.data.data_vars:
    airports_clc = gpd.read_file(atlas.path / 'data' / 'site' / 'clc' / 'CLC_2018_AT.shp')
    airports_clc['CODE_18'] = airports_clc['CODE_18'].astype('int')
    airports_clc = airports_clc.loc[airports_clc['CODE_18'] == 124].to_crs(atlas.crs)
    airports_clc = airports_clc.clip(atlas.get_nuts_region(atlas.region).geometry)

    airport_buffers = {
        'Vöslau': 2500,
        'Wr. Neustadt/West MIL': 5000,
        'Wr. Neustadt/Ost': 5000,
        'Wien-Schwechat': 6250,
        'Tulln MIL': 5000,
    }
    # airports - austro control
    fname_most_recent_airspace = sorted((atlas.path / 'data' / 'site' / 'airspace').glob('*.kml'), reverse=True)[0]
    airports_ac = gpd.read_file(fname_most_recent_airspace, driver='LIBKML', layer='Airports')
    airports_ac = airports_ac.to_crs(atlas.crs)

    # merge airports and apply buffers
    airports = airports_clc.sjoin(airports_ac, predicate="contains", how="left")
    for airport in list(airports['Name']):
        airports.loc[airports['Name'] == airport, 'geometry'] = \
            airports[airports['Name'] == airport].buffer(airport_buffers[airport])

    atlas.landscape.rasterize(airports, name='airports_buff', all_touched=touch)
    del airport, airports, airports_ac, airports_clc

# %% restricted military areas
if "restriceted_military_areas" not in atlas.landscape.data.data_vars:
    fname_most_recent_airspace = sorted((atlas.path / 'data' / 'site' / 'airspace').glob('*.kml'), reverse=True)[0]
    mil_areas = gpd.read_file(fname_most_recent_airspace, driver='LIBKML', layer='Restricted areas military (R)')
    atlas.landscape.rasterize(mil_areas, name='restricted_military_areas', all_touched=touch)
    del mil_areas, fname_most_recent_airspace

# %% grid infrastructure
if "power_lines" not in atlas.landscape.data.data_vars:
    lines = gpd.read_file(atlas.path / 'data' / 'site' / 'buildings' / 'BAU_2700_STROMLEITUNG_L_0.shp')
    lines = lines.to_crs(atlas.crs)
    lines = lines.clip(atlas.get_nuts_region(atlas.region))
    lines = lines.dissolve()
    atlas.landscape.rasterize(lines, name='power_lines', all_touched=touch)
    del lines

# %% high-ranking roads
if "roads" not in atlas.landscape.data.data_vars:
    roads = process_gip(atlas.path / 'data' / 'site' / 'gip')
    atlas.landscape.rasterize(roads, name='roads', all_touched=touch)
    del roads

# %% buildings in the greenland worth preserving
# load RRU_WI_GEBPoint.shp
if "buildings_in_greenland" not in atlas.landscape.data.data_vars:
    gruengeb = gpd.read_file(atlas.path / 'data' / 'site' / 'landuse' / 'Erhaltenswerte_gebaeude.shp', encoding='utf-8')
    gruengeb = gruengeb.dissolve()
    atlas.landscape.rasterize(gruengeb, name='buildings_in_greenland', all_touched=True)
    del gruengeb

# %% municipal zoning
if "residential_buildings" not in atlas.landscape.data.data_vars:
    # load RRU_WI_HUELLEPolygon_0.shp
    municipal_zoning = gpd.read_file(atlas.path / 'data' / 'site' / 'landuse' / 'Widmungen_noe.shp', encoding='utf-8')
    municipal_zoning = municipal_zoning.dissolve(by="WIDMUNG")
    municipal_zoning.reset_index(inplace=True)

    # residential buildings
    residential_cats = [
        'Bauland Wohn- oder Mischnutzung',
        'Bauland Wohn- oder Mischnutzung nachhaltig',
    ]
    residential_zoning = municipal_zoning.loc[municipal_zoning['WIDMUNG'].str.contains('|'.join(residential_cats)), :]
    atlas.landscape.rasterize(residential_zoning, name='residential_buildings', all_touched=touch)

    # other building land
    other_building_cats = [
        'Bauland Betriebsnutzung', 'Bauland-Sondergebiet', 'Bauland-Kerngebiet Handelseinrichtungen'
    ]
    # 'Bauland Wohn- oder Mischnutzung',
    bau = municipal_zoning.loc[municipal_zoning['WIDMUNG'].str.contains('|'.join(other_building_cats)), :]
    atlas.landscape.rasterize(bau, name='other_building_land', all_touched=touch)

    # greenland zoning
    greenland_zoning_cats = [
        'Grünland-Campingplatz', 'Grünland-Kleingarten', 'Grünland-land- und forstwirtschaftliche Hofstelle'
    ]
    greenland_zoning = municipal_zoning.loc[municipal_zoning['WIDMUNG'].str.contains('|'.join(greenland_zoning_cats)),
                       :]
    atlas.landscape.rasterize(greenland_zoning, name='greenland_zonings', all_touched=True)
    del municipal_zoning, greenland_zoning_cats, other_building_cats, residential_cats
    del greenland_zoning, residential_zoning, bau

# %% compute distances
if "distance_to_roads" not in atlas.landscape.data.data_vars:
    distances_to_compute = [
        "roads",
        "power_lines",
        "existing_turbines",
        "buildings_in_greenland",
        "residential_buildings",
        "other_building_land",
        "greenland_zonings"
    ]
    atlas.landscape.compute_distance(distances_to_compute, inplace=True)

# %% convert all distances to kilometers
distance_vars = [var for var in atlas.landscape.data.data_vars if "distance" in var]
atlas.landscape.convert(distance_vars, to_unit="km", inplace=True)

# %% save to disk
atlas.wind.data.to_netcdf(path=atlas.path/'data'/'processed'/'WindAtlas_AUT_pro.nc')
atlas.landscape.data.to_netcdf(path=atlas.path/'data'/'processed'/'LandscapeAtlas_AUT_pro.nc')

# %% wind turbines as of 2014
atlas.wind.simulate_capacity_factors(bias_correction=0.71197)
atlas.wind.compute_lcoe(turbine_cost_share=0.7)
atlas.wind.minimum_lcoe()
atlas.wind.compute_optimal_power_energy()
atlas.save(scenario='2014')

# write data at time of decision to csv
atlas.landscape.add(atlas.wind.data.min_lcoe, name='min_lcoe')
# atlas.landscape.add(atlas.data.optimal_power)
# atlas.landscape.add(atlas.data.optimal_energy)

df = atlas.landscape.flatten()
df = df.dropna(how="any", axis=0)
df.to_csv(atlas.path / "data" / "processed" / "dc_data_2014.csv")
del df

# %% modern wind turbines
atlas.wind_turbines = [
    "Enercon.E138.3500", "Enercon.E160.5560", "Vestas.V150.4200",
]

atlas.wind.data = atlas.wind.data.drop_vars(["capacity_factors", "lcoe", "min_lcoe"])
atlas.wind.simulate_capacity_factors(bias_correction=0.71197)
atlas.wind.compute_lcoe(turbine_cost_share=0.7)
atlas.wind.minimum_lcoe()

atlas.landscape.data = atlas.landscape.data.drop_vars("min_lcoe")
atlas.landscape.add(atlas.wind.data.min_lcoe, name='min_lcoe')
atlas.save(scenario='modern')
df = atlas.landscape.flatten()
df = df.dropna(how="any", axis=0)
df.to_csv(atlas.path / "data" / "processed" / "dc_data_modern.csv")
