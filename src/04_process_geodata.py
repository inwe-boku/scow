import sys;

print('Python %s on %s' % (sys.version, sys.platform))
sys.path.extend(['C:\\git_repos\\vala', 'C:\\git_repos\\vala', 'C:/git_repos/vala'])

# %% imports
import logging
import numpy as np
import pandas as pd
import geopandas as gpd
import datetime as dt
import xarray as xr
import rioxarray as rxr

from fiona import drvsupport
from shapely import wkt
from config import ROOTDIR, TURBINE_YEARS, country
from src.vala import clip_array_to_polygon, rasterize_shp2xr, dataset_to_pandas, raster_prox, combine_shapefiles
from logging_config import setup_logging

setup_logging()

# %% settings
INTERPOLATE = False
BL = 'Niederösterreich'
touch = True
qual = 'cost'

zones = {
    'Niederösterreich': ROOTDIR / 'data/zones/Zonierung_noe.shp',
}

port_buffs = {
    'Vöslau': 2500,
    'Wr. Neustadt/West MIL': 5000,
    'Wr. Neustadt/Ost': 5000,
    'Wien-Schwechat': 6250,
    'Tulln MIL': 5000
}

# %% functions

# %% read data
logging.info('Starting to process geodata')

land = gpd.read_file(ROOTDIR / 'data/vgd/vgd_oesterreich.shp')
land = land[['BL', 'geometry']].dissolve(by='BL')
land.reset_index(inplace=True)
land = land.loc[land.BL == BL]

xrtemplate = xr.open_dataarray(ROOTDIR / f'data/preprocessed/gwa_roughness_{country}.nc')
xrtemplate = xrtemplate.squeeze()
xrtemplate = clip_array_to_polygon(xrtemplate, land, 'template')

# %% wind power zoning
zoning = gpd.read_file(zones[BL])
zoning = zoning.to_crs(land.crs)
zoning.reset_index(inplace=True)
zoning_land = rasterize_shp2xr(xrtemplate, land, dummyshape=zoning, crs=land.crs, name='zoning', all_touched=touch)

geodat = xr.merge([zoning_land])
del zoning

logging.info('Zoning processed: [zoning_land]')

# %% existing turbines
turbs = gpd.read_file(ROOTDIR / 'data/AT_turbines/turbines.shp')
# turbine locations in 2014
turbines_land = rasterize_shp2xr(xrtemplate, land, dummyshape=turbs.loc[turbs['Jahr'] < 2015, :],
                                 crs=land.crs, name='existing_turbines', all_touched=touch)

turbines_prox_land = raster_prox(turbines_land, xrtemplate, 'prx_turbines')
turbines_prox_land = turbines_prox_land / 1000
turbines_prox_land.attrs['unit'] = 'km'

geodat = xr.merge([geodat, turbines_land, turbines_prox_land])
del turbs, turbines_land, turbines_prox_land

logging.info('Wind turbine data processed: [existing_turbines_land, prx_turbines]')

# %% landscapes worthy of preservation - "Erhaltenswerte Landschaftsteile Niederösterreich"
erh = gpd.read_file(ROOTDIR / 'data/preservation/RRU_REGROP_ELTPolygon_0.shp')
erh = rasterize_shp2xr(xrtemplate, land, dummyshape=erh, crs=land.crs, name='preserve', all_touched=touch)

geodat = xr.merge([geodat, erh])
del erh

logging.info('Preservation areas processed: [preserve]')

# %% protected area polygons
wdpa_date = dt.datetime.today().strftime('%b%Y')
wdpa = combine_shapefiles(ROOTDIR / 'data/schutzgebiete', f'WDPA_WDOECM_{wdpa_date}_Public_AUT_shp-polygons', [1, 2, 3])
wdpa = wdpa.to_crs(land.crs)
protareas_land = rasterize_shp2xr(xrtemplate, land, dummyshape=wdpa, crs=land.crs, name='prt_areas', all_touched=touch)

"""
# split up protected areas by IUCN categories:
iucn_cats = ['Ia', 'Ib', 'II', 'III', 'IV', 'V', 'VI']
wdpa_subcats = ['Birds', 'Habitats', 'Ramsar']
protection_categories = wdpa_categories(wdpa, iucn_cats, wdpa_subcats)
"""

geodat = xr.merge([geodat, protareas_land])
del wdpa

logging.info('Protected areas processed: [prt_areas]')

# %% Important Bird Areas - BirdLife
# aiba = gpd.read_file(ROOTDIR / 'data/iba/IBA bounbdaries Austria 6 9 2021.shp')  # Source: BirdLife Austria
iba = gpd.read_file(ROOTDIR / 'data/iba/EUROPE_KBA.shp')  # Source: BirdLife International
iba = iba.loc[iba['ISO3'] == 'AUT', :]
iba_land = rasterize_shp2xr(xrtemplate, land, dummyshape=iba, crs=land.crs, name='brd_areas', all_touched=touch)

geodat = xr.merge([geodat, iba_land])
del iba, iba_land

logging.info('Important Bird Areas processed: [brd_areas]')

# %% alpine convention
alpconv = gpd.read_file(ROOTDIR / 'data/alps/alpsconvention.shp')
alps = rasterize_shp2xr(xrtemplate, land, dummyshape=alpconv, crs=land.crs, name='alpconv', all_touched=touch)

geodat = xr.merge([geodat, alps])
del alpconv, alps

# %% forest areas
trcvdnsty = rxr.open_rasterio(ROOTDIR / 'data/forest/2015/TCD_2015_020m_at_31287_d02/TCD_2015_020m_at_31287_d02.tif')
trcvdnsty = trcvdnsty / 100
trcvdnsty = clip_array_to_polygon(trcvdnsty, land, 'tree_cover_density', interp_array=xrtemplate)
trcvdnsty = trcvdnsty.where(~xrtemplate.isnull(), np.nan)
trcvdnsty = xr.where(trcvdnsty > 1, 1, trcvdnsty)
trcvdnsty.attrs['unit'] = '%'

# leaf type: 0 - not tree covered, 1 - broadleaved trees, 2 - coniferous trees, >2 - unclassifiable
lftyp = rxr.open_rasterio(ROOTDIR / 'data/forest/2015/DLT_2015_020m_at_31287_d02/DLT_2015_020m_at_31287_d02.tif')
lftyp = lftyp.where(~(lftyp == 255), np.nan)
lftyp = clip_array_to_polygon(lftyp, land, 'leaf_type', interp_array=xrtemplate)

dlt_broadleaf = xr.where(lftyp == 1, 1, 0)
dlt_broadleaf = xr.where(xrtemplate == np.nan, np.nan, dlt_broadleaf)
dlt_broadleaf.name = 'broadleaved'

dlt_coniferous = xr.where(lftyp == 2, 1, 0)
dlt_coniferous = xr.where(xrtemplate == np.nan, np.nan, dlt_coniferous)
dlt_coniferous.name = 'coniferous'

geodat = xr.merge([geodat, trcvdnsty, lftyp, dlt_broadleaf, dlt_coniferous])
del trcvdnsty, lftyp, dlt_broadleaf, dlt_coniferous

logging.info('Forest data processed: [tree_cover_density, leaf_type, broadleaved, coniferous]')

# %% Corine Land Cover - pastures and crop area
clc = gpd.read_file(ROOTDIR / 'data/clc/CLC_2018_AT.shp')
clc['CODE_18'] = clc['CODE_18'].astype('int')
clc = clc.to_crs(land.crs)

# pastures
pastures = clc[clc['CODE_18'] == 231]
pastures_land = rasterize_shp2xr(xrtemplate, land, dummyshape=pastures, crs=land.crs, name='pastures',
                                 all_touched=touch)
# crop areas
crops = clc[clc['CODE_18'].isin([211, 212, 213, 221, 222, 223, 241, 242, 243])]
crops_land = rasterize_shp2xr(xrtemplate, land, dummyshape=crops, crs=land.crs, name='crops', all_touched=touch)

geodat = xr.merge([geodat, pastures_land, crops_land])
del pastures, pastures_land, crops, crops_land

logging.info('Corine Land Cover processed: [pastures, crops]')

# %% Airspace - airports and restricted military areas
drvsupport.supported_drivers['LIBKML'] = 'rw'
airspace = ROOTDIR / 'data/airspace/20231228LuftraumAT.kml'

# airports - corine land cover
airports_clc = clc[clc['CODE_18'] == 124]
airports_clc = airports_clc.clip(land)

# airports - austro control
airports_ac = gpd.read_file(airspace, driver='LIBKML', layer='Airports')
airports_ac = airports_ac.to_crs(land.crs)
# merge airports and apply buffers
airports = airports_clc.sjoin(airports_ac, predicate="contains", how="left")
for port in list(airports['Name']):
    airports.loc[airports['Name'] == port, 'geometry'] = \
        airports[airports['Name'] == port].buffer(port_buffs[port])

airports_buff = rasterize_shp2xr(xrtemplate, land, dummyshape=airports, crs=land.crs, name='airports_buff',
                                  all_touched=touch)

# restricted military areas
mil_areas = gpd.read_file(airspace, driver='LIBKML', layer='Restricted areas military (R)')
mil_areas = mil_areas.to_crs(land.crs)
mil_areas = gpd.clip(mil_areas, land)
res_mil_areas = rasterize_shp2xr(xrtemplate, land, dummyshape=mil_areas, crs=land.crs, name='mil_areas',
                                 all_touched=touch)

geodat = xr.merge([geodat, airports_buff, res_mil_areas])
del airports, airports_buff, mil_areas, res_mil_areas

logging.info('Airspace data processed: [airports, mil_areas]')

# %% municipal zoning ("widmungen")
umh = gpd.read_file(ROOTDIR / 'data/landuse/Widmungen_noe.shp', encoding='utf-8')  # RRU_WI_HUELLEPolygon_0.shp
# wohnwidmung
wohnwidmung = umh.loc[umh['WIDMUNG'].str.contains('Bauland Wohn- oder Mischnutzung'), :]
wohnwidmung = rasterize_shp2xr(xrtemplate, land, dummyshape=wohnwidmung, crs=land.crs,
                               name='wohnwidmung', all_touched=touch)
wohnwidmung_prox = raster_prox(wohnwidmung, xrtemplate, 'prx_wohnwidmung')
wohnwidmung_prox = wohnwidmung_prox / 1000
wohnwidmung_prox.attrs['unit'] = 'km'
# bauland-widmungen
bauwidmung = ['Bauland Betriebsnutzung', 'Bauland-Sondergebiet', 'Bauland-Kerngebiet Handelseinrichtungen']
# 'Bauland Wohn- oder Mischnutzung',
bau = umh.loc[umh['WIDMUNG'].str.contains('|'.join(bauwidmung)), :]
bauland = rasterize_shp2xr(xrtemplate, land, dummyshape=bau, crs=land.crs, name='bauland', all_touched=touch)
bauland_prox = raster_prox(bauland, xrtemplate, 'prx_bauland')
bauland_prox = bauland_prox / 1000
bauland_prox.attrs['unit'] = 'km'

# greenland zonings ("grünland-widmungen")
gruenwidmungen = ['Grünland-Campingplatz', 'Grünland-Kleingarten', 'Grünland-land- und forstwirtschaftliche Hofstelle']
gruenwid = umh.loc[umh['WIDMUNG'].str.contains('|'.join(gruenwidmungen)), :]
greenzone = rasterize_shp2xr(xrtemplate, land, dummyshape=gruenwid, crs=land.crs, name='greenzoning', all_touched=True)

greenzone_prox = raster_prox(greenzone, xrtemplate, 'prx_greenzoning')
greenzone_prox = greenzone_prox / 1000
greenzone_prox.attrs['unit'] = 'km'

# buildings in greenland worthy of protection ("erhaltenswerte gebäude im grünland")
gruengeb = gpd.read_file(ROOTDIR / 'data/landuse/Erhaltenswerte_gebaeude.shp', encoding='utf-8')  # RRU_WI_GEBPoint.shp
greenbuild = rasterize_shp2xr(xrtemplate, land, dummyshape=gruengeb, crs=land.crs, name='greenbuildings', all_touched=True)
greenbuild_prox = raster_prox(greenbuild, xrtemplate, 'prx_greenbuildings')
greenbuild_prox = greenbuild_prox / 1000
greenbuild_prox.attrs['unit'] = 'km'

geodat = xr.merge([geodat, bauland, bauland_prox, greenzone, greenzone_prox, greenbuild, greenbuild_prox,
                   wohnwidmung, wohnwidmung_prox])
del umh, bauwidmung, bau, bauland, bauland_prox, gruenwidmungen, gruenwid, gruengeb, greenzone, greenzone_prox
del greenbuild, greenbuild_prox, wohnwidmung, wohnwidmung_prox

logging.info('Municipal zoning processed: [bauland, prx_bauland, greenland, prx_greenland, wohnwidmung, prx_wohnwidmung]')

# %% tourism
overnights = gpd.read_file(ROOTDIR / 'data/tourism/Naechtigungen.shp')
overnights = overnights.to_crs(land.crs)
overnights_land = rasterize_shp2xr(xrtemplate, land, dummyshape=overnights, crs=land.crs, column='Naechtigun',
                                   name='overnights')
overnights_land = overnights_land / 1000
overnights_land = overnights_land.interp_like(xrtemplate)
overnights_land = xr.where(overnights_land < 0.0001, 0, overnights_land)

geodat = xr.merge([geodat, overnights_land])
del overnights, overnights_land

logging.info('Touristic overnights processed: [overnights]')

# %% roads
roads = gpd.read_file(ROOTDIR / 'data/gip/hrng_streets.shp')
roads_land = rasterize_shp2xr(xrtemplate, land, dummyshape=roads, crs=land.crs,
                              name='roads', all_touched=touch)

roads_prox_land = raster_prox(roads_land, xrtemplate, 'prx_roads')
# convert distances to km
roads_prox_land = roads_prox_land / 1000
roads_prox_land.attrs['unit'] = 'km'

geodat = xr.merge([geodat, roads_land, roads_prox_land])
del roads, roads_land, roads_prox_land

logging.info('Road network processed: [roads, prx_roads]')

# %% water bodies
waters = pd.concat([gpd.read_file(ROOTDIR / 'data/water_bodies/main_standing_waters.shp'),
                    gpd.read_file(ROOTDIR / 'data/water_bodies/main_running_waters.shp')])
waters_land = rasterize_shp2xr(xrtemplate, land, dummyshape=waters, crs=land.crs,
                               name='waters', all_touched=touch)

geodat = xr.merge([geodat, waters_land])
del waters, waters_land

logging.info('Water bodies processed: [waters, prx_waters]')

# %% distance to high-voltage grid
lines = gpd.read_file(ROOTDIR  / 'data/buildings/BAU_2700_STROMLEITUNG_L_0.shp')
lines = lines.to_crs(land.crs)
lines_noe = lines.clip(land)
lines_land = rasterize_shp2xr(xrtemplate, land, dummyshape=lines_noe, crs=land.crs, name='lines', all_touched=touch)
lines_prox_land = raster_prox(lines_land, xrtemplate, 'prx_grid')

# convert distances to km
lines_prox_land = lines_prox_land / 1000
lines_prox_land.attrs['unit'] = 'km'

geodat = xr.merge([geodat, lines_land, lines_prox_land])
del lines, lines_noe, lines_prox_land

logging.info('High voltage grid processed: [lines, prx_grid]')

# %% terrain slope
slope = xr.open_dataset(ROOTDIR / 'data/elevation/slope_31287.nc')
slope = slope['slope'].rio.write_crs('epsg:31287')
slope_land = clip_array_to_polygon(slope, land, 'slope')
slope_land = slope_land.interp_like(xrtemplate)

elev = xr.open_dataset(ROOTDIR / 'data/elevation/elevation_31287.nc')
elev = elev['elevation'].rio.write_crs('epsg:31287')
elev_land = clip_array_to_polygon(elev, land, 'elevation')
elev_land = elev_land.interp_like(xrtemplate)
elev_land.attrs['unit'] = 'm'

geodat = xr.merge([geodat, slope_land, elev_land])
del slope, slope_land, elev, elev_land

logging.info('Terrain slope processed: [slope, elevation]')

# %% LCOE according to turbine year
for turbine_year in TURBINE_YEARS:
    geoloop = geodat

    lcoe = xr.open_dataarray(ROOTDIR / f'data/preprocessed/lcoe_{country}_{turbine_year}.nc', mask_and_scale=True)
    lcoe = lcoe.rio.reproject(xrtemplate.rio.crs)
    lcoe_land = clip_array_to_polygon(lcoe, land, 'lcoe')
    # lcoe_land = lcoe_land.sel(turbine_models='min_lcoe')  # process only min_lcoe!
    lcoe_land = lcoe_land.interp_like(xrtemplate)

    geoloop = xr.merge([geoloop, lcoe_land])
    del lcoe

    logging.info('Wind resource data processed: [lcoe_land]')

    # save dataset
    if touch:
        name_ext = f'touch_{qual}_{turbine_year}'
    else:
        name_ext = f'notouch_{qual}_{turbine_year}'

    fname = ROOTDIR / f'data/preprocessed/geodat_{BL}_{name_ext}'.replace('ö', 'oe')
    geoloop.to_netcdf(path=ROOTDIR / f'{fname}.nc'.replace('ö', 'oe'), mode='w', format='NETCDF4', engine='netcdf4')
    # geoloop = xr.open_dataset(ROOTDIR / f'{fname}.nc'.replace('ö', 'oe'))

    # %% export tidy geodata
    # fname = ROOTDIR / f'data/preprocessed/geodat_{BL}_{name_ext}.csv'.replace('ö', 'oe')
    tidyvars = dataset_to_pandas(geoloop, digits=4, drop_labels=['band', 'spatial_ref', 'turbine_models'])
    tidyvars = tidyvars.dropna(how='any', axis=0)  # drops 1.7% of raster cells
    tidyvars.to_csv(f'{fname}.csv')

    logging.info(f'Geodata for {BL} exported to {fname}')

logging.info("Geodata processing completed")
