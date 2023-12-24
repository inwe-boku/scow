# %% imports
import os
import logging
import numpy as np
import pandas as pd
import geopandas as gpd
import rioxarray as rxr
import subprocess
from xrspatial import slope
from scipy.special import gamma
from config import ROOTDIR, turbines, country
from src.vala import process_power_curves
from logging_config import setup_logging

setup_logging()

# %% read data
logging.info('Preprocessing started')
# power curves
rnj_power_curves = pd.read_csv(ROOTDIR / 'data/power_curves/rnj_power_curve_000-smooth.csv')
oep_power_curves = pd.read_csv(ROOTDIR / 'data/power_curves/supply__wind_turbine_library.csv')
own_power_curves = pd.read_csv(ROOTDIR / 'data/power_curves/powercurves_research.csv', sep=';', decimal=',',
                               skiprows=[1, 2])
own_power_curves = own_power_curves.set_index('speed')
# air density and roughness
elevation = rxr.open_rasterio(ROOTDIR / f'data/gwa3/{country}_elevation_w_bathymetry.tif')
a50 = rxr.open_rasterio(ROOTDIR / f'data/gwa3/{country}_combined-Weibull-A_50.tif')
k50 = rxr.open_rasterio(ROOTDIR / f'data/gwa3/{country}_combined-Weibull-k_50.tif')
a100 = rxr.open_rasterio(ROOTDIR / f'data/gwa3/{country}_combined-Weibull-A_100.tif')
k100 = rxr.open_rasterio(ROOTDIR / f'data/gwa3/{country}_combined-Weibull-k_100.tif')
# tourism
austria = gpd.read_file(ROOTDIR / 'data/vgd/vgd_oesterreich.shp')
austria = austria[['GKZ', 'PG', 'BL', 'geometry']].dissolve(by='GKZ')
austria.reset_index(inplace=True)
austria['GKZ'] = austria['GKZ'].astype('int')

tour = pd.read_excel(ROOTDIR / 'data/tourism/Tabelle 30 GEH 2018.xlsx', header=[0, 1, 2, 3], skiprows=[0, 5, 6],
                     na_values=['GEH'])
tour = tour.loc[~tour.iloc[:, 0].isna(), :]
tour = tour.loc[[not isinstance(e, str) for e in tour.iloc[:, 0]], :]

# %% turbine power curves
oep_power_curves['type_string'] = oep_power_curves['manufacturer'] + '.' + oep_power_curves['turbine_type'].replace({'/': '.', '-': ''}, regex=True)
oep_power_curves.dropna(subset=['power_curve_values'], inplace=True)

powercurves = process_power_curves(turbines, rnj_power_curves, oep_power_curves)
powercurves = pd.concat([powercurves, own_power_curves], axis=1)

if not os.path.exists(ROOTDIR / 'data/preprocessed'):
    os.mkdir(ROOTDIR / 'data/preprocessed')
powercurves.to_csv(ROOTDIR / 'data/preprocessed/powercurves.csv', sep=';', decimal=',')
logging.info('Powercurves preprocessed')

# %% compute air density correction factor from elevation data
# compute air density correction - see https://wind-data.ch/tools/luftdichte.php
rho_correction = 1.247015 * np.exp(-0.000104 * elevation) / 1.225
rho_correction.to_netcdf(path=ROOTDIR / f'data/preprocessed/gwa_air_density_{country}.nc')

logging.info('Air density correction preprocessed')

# %% compute roughness factor alpha
u_mean_50 = a50 * gamma(1 / k50 + 1)
u_mean_100 = a100 * gamma(1/k100 + 1)
alpha = (np.log(u_mean_100) - np.log(u_mean_50)) / (np.log(100) - np.log(50))
alpha.to_netcdf(path=ROOTDIR / f'data/preprocessed/gwa_roughness_{country}.nc',
                mode='w', format='NETCDF4', engine='netcdf4')

logging.info('Roughness preprocessed')

# %% process tourism data
col = pd.DataFrame.from_records(data=tour.columns.to_flat_index())
col = col.replace({"-\n": "", "\n": " ", " ": ""}, regex=True)
for i in range(0, 4):
    col.loc[col[i].str.contains('Unnamed'), i] = ''
col = [col.loc[row, :].str.cat(sep=' ').strip() for row in col.index]
tour.columns = col
tour['GEM.KENNZIFFER'] = tour['GEM.KENNZIFFER'].astype('int')
tour['WINTERSAISON2017/2018 ÜBERNACHTUNGEN INSGESAMT'] = tour[['WINTERSAISON2017/2018 ÜBERNACHTUNGEN INSGESAMT']].replace('-', 0).astype('float')
tour['SOMMERSAISON2018 ÜBERNACHTUNGEN INSGESAMT'] = tour[['SOMMERSAISON2018 ÜBERNACHTUNGEN INSGESAMT']].astype('float').fillna(0)

tour['Naechtigungen'] = tour[['WINTERSAISON2017/2018 ÜBERNACHTUNGEN INSGESAMT', 'SOMMERSAISON2018 ÜBERNACHTUNGEN INSGESAMT']].sum(axis=1)

stays = austria.merge(tour[['GEM.KENNZIFFER', 'Naechtigungen']], left_on='GKZ', right_on='GEM.KENNZIFFER', how='left')
stays['Naechtigungen'] = stays['Naechtigungen'].fillna(0)
stays.to_file(ROOTDIR / 'data/tourism/Naechtigungen.shp')

logging.info('Touristic overnight stays preprocessed')

# %% process GIP geopackage
# documented in https://www.gip.gv.at/assets/downloads/GIP_Standardbeschreibung_2.3.2_FINAL.pdf
# street categories are explained in 2.3.2.6 Abschnittskategorie (GIP.EDGE.EDGECATEGORY)
# convert to shapefiles
os.chdir(ROOTDIR / 'data/gip/')
subprocess.run("""ogr2ogr -f "ESRI shapefile" shp gip_network_ogd.gpkg -sql "select cast(EDGECAT as character(255)) from EDGE_OGD" -overwrite -dialect ogrsql""")

# read data
gip_shp = ROOTDIR / 'data/gip/shp/EDGE_OGD.shp'
# abbreviations:
# A: Autobahn / highway
# S: Schnellstraße / highway
# B: Landesstraße B (ehem. Bundesstraße)
# L: Landesstraße L (ehem. Landeshauptstraße)
hochrangig = ['A', 'S', 'B', 'L']
hrng = gpd.GeoDataFrame()

# total number of rows: 1 532 485
for n in range(1, 63):
    gip = gpd.read_file(gip_shp, rows=slice((n-1)*25000, n*25000))
    gip = gip.loc[gip['EDGECAT'].isin(hochrangig)]
    hrng = pd.concat([hrng, gip])
    gip = []

hrng.to_file(ROOTDIR / 'data/gip/hrng_streets.shp')

logging.info('High-level road grid preprocessed')

# %% process water bodies
main_water_bodies = ['100 km² Gewässer', '1000 km² Gewässer', '10000 km² Gewässer', '500 km² Gewässer', '4000 km² Gewässer']
wbd = gpd.read_file(ROOTDIR / 'data/water_bodies/Fliessgewaesser.shp')
wbd = wbd.loc[wbd['GEW_KAT'].isin(main_water_bodies)]
wbd = wbd.dissolve()
wbd.to_file(ROOTDIR / 'data/water_bodies/main_running_waters.shp')

lks = gpd.read_file(ROOTDIR / 'data/water_bodies/stehendeGewaesser.shp')
lks = lks.loc[lks['FLAECHEKM2'] >= 0.03125, :]
lks.to_file(ROOTDIR / 'data/water_bodies/main_standing_waters.shp')

logging.info('Water bodies preprocessed')

# %% compute terrain slope
#states = austria[['BL', 'geometry']].dissolve(by='BL')
#states.reset_index(inplace=True)

#A100 = rxr.open_rasterio(ROOTDIR / f'data/gwa3/{country}_combined-Weibull-A_100.tif')
#A100 = A100.squeeze()
#A100 = A100.rio.reproject(states.crs)

elevation = rxr.open_rasterio(ROOTDIR / 'data/elevation/dhm_at_lamb_10m_2018.tif')
elevation = elevation.squeeze()
elevation = elevation.rio.reproject(austria.crs)
elevation.name = 'elevation'
elevation.to_netcdf(path=ROOTDIR / 'data/elevation/elevation_31287.nc')

slp = slope(elevation)
slp = slp.rio.reproject(austria.crs)
#slp = slp.interp_like(A100)
#slp = slp.where(~A100.isnull(), np.nan)
slp.to_netcdf(path=ROOTDIR / 'data/elevation/slope_31287.nc')

#elevation = elevation.interp_like(A100)
#elevation = elevation.where(~A100.isnull(), np.nan)


logging.info('Terrain elevation and slope preprocessed')

# %% merge wind turbine location data
bev = gpd.read_file(ROOTDIR / 'data/buildings/BAU_2200_BETRIEB_P_0.shp')
bev = bev.loc[bev.F_CODE == 2202, :]
bev = bev.to_crs(austria.crs)
bev['points'] = bev['geometry']
bev['geometry'] = bev['geometry'].buffer(20)

igw = pd.read_csv(ROOTDIR / 'data/AT_turbines/igwturbines.csv', sep=';', decimal=',')
igw = gpd.GeoDataFrame(data=igw, geometry=gpd.points_from_xy(igw.lon, igw.lat), crs='epsg:4326')
igw = igw.to_crs(austria.crs)

turbine_locations = bev.sjoin(igw, how='left', predicate='contains')
turbine_locations = turbine_locations.reset_index(drop=True)
turbine_locations.geometry = turbine_locations.points
turbine_locations['ERSTELLDAT'] = pd.to_datetime(turbine_locations['ERSTELLDAT'])
turbine_locations.loc[turbine_locations['Jahr'].isna(), 'Jahr'] = turbine_locations.loc[turbine_locations['Jahr'].isna(), 'ERSTELLDAT'].dt.year
turbine_locations['Jahr'] = turbine_locations['Jahr'].astype('int')
turbine_locations['ERSTELLDAT'] = turbine_locations['ERSTELLDAT'].astype('str')
turbine_locations = turbine_locations[['NAME', 'Name', 'BETREIBER', 'Betreiber1', 'Betreiber2', 'HOEHE', 'Nabenhöhe',
               'Rotordurchmesser', 'kW', 'Hersteller', 'Type', 'Jahr', 'ERSTELLDAT', 'n_Anlagen', 'url',
               'geometry']]
turbine_locations.to_file(ROOTDIR / 'data/AT_turbines/turbines.shp')

logging.info('Turbine locations preprocessed')

# %% create results directory
if not os.path.exists(ROOTDIR / 'data/results'):
    os.mkdir(ROOTDIR / 'data/results')
    logging.info('Directory for results created')

logging.info('preprocessing complete')
