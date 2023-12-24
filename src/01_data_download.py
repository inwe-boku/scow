# %% imports
import os
import json
import zipfile
import numpy as np
import pandas as pd
import logging
import datetime as dt

from config import ROOTDIR, country
from src.vala import download_file
from logging_config import setup_logging

setup_logging()

# %% proxy settings
wdpa_date = dt.datetime.today().strftime('%b%Y')

proxy = None  # 'https://proxy.wstw.energy-it.net:8080'
proxy_user = None  # 'DWSTW/netw4j'
proxy_pass = None  # 'Sch0lterzuck?'


# %% global wind atlas settings
url_gwa = 'https://globalwindatlas.info/api/gis/country'
layer = ['air-density', 'combined-Weibull-A', 'combined-Weibull-k']
ground = ['elevation_w_bathymetry']
height = ['50', '100', '150']

# %% file, directory, url dict
data_dict = {
    'bev_bauten.zip': [ROOTDIR / 'data/buildings', 'https://data.bev.gv.at/download/DLM/DLM_20230125/DLM_2000_BAUTEN_20230125.zip'],
    'rnj_power_curve_000-smooth.csv': [ROOTDIR / 'data/power_curves', 'https://raw.githubusercontent.com/renewables-ninja/vwf/master/power_curves/Wind%20Turbine%20Power%20Curves%20%7E%205%20(0.01ms%20with%200.00%20w%20smoother).csv'],
    'oep_wind_turbine_library.zip': [ROOTDIR / 'data/power_curves', """ "https://openenergy-platform.org/api/v0/schema/supply/tables/wind_turbine_library/rows/?form=datapackage" """],
    'vgd_oesterreich.zip': [ROOTDIR / 'data/vgd', """ "https://data.bev.gv.at/download/Verwaltungsgrenzen/shp/20211001/VGD_Oesterreich_gst_20211001.zip" """],
    'CLC_2018_AT.zip': [ROOTDIR / 'data/clc', 'https://docs.umweltbundesamt.at/s/beBw8fmwyCMA2ga/download/CLC_2018_AT_clip.zip'],
    'gridkit_europe.zip': [ROOTDIR / 'data/grid', 'https://zenodo.org/record/47317/files/gridkit_euorpe.zip?download=1'],
    'Zonierung_noe.zip': [ROOTDIR / 'data/zones', """ "https://sdi.noe.gv.at/at.gv.noe.geoserver/OGD/wfs?request=GetFeature&version=1.1.0&typeName=OGD:RRU_WIND_ZONEN_P19&srsName=EPSG:31259&outputFormat=shape-zip&format_options=CHARSET:UTF-8" """],
    'B_gip_network_ogd.zip': [ROOTDIR / 'data/gip', 'https://open.gip.gv.at/ogd/B_gip_network_ogd.zip'],
    'Fliessgewaesser.zip': [ROOTDIR / 'data/water_bodies', 'https://docs.umweltbundesamt.at/s/YkgTDiDs9DPstCJ/download/Routen_v16.zip'],
    'StehendeGewaesser.zip': [ROOTDIR / 'data/water_bodies', 'https://docs.umweltbundesamt.at/s/t4jHoXmrwrsjnea/download/stehendeGewaesser_v16.zip',],
    'Adressregister.zip': [ROOTDIR / 'data/gwr', 'https://data.bev.gv.at/download/Adressregister/Archiv_Adressregister/Adresse_Relationale_Tabellen_Stichtagsdaten_20220403.zip'],
    'elevation.zip': [ROOTDIR / 'data/elevation', 'https://gis.ktn.gv.at/OGD/Geographie_Planung/ogd-10m-at.zip'],
    'Widmungen_noe.zip': [ROOTDIR / 'data/landuse', 'https://sdi.noe.gv.at/at.gv.noe.geoserver/OGD/wfs?request=GetFeature&version=1.1.0&typeName=OGD:RRU_WI_HUELLE&srsName=EPSG:31259&outputFormat=shape-zip&format_options=CHARSET:UTF-8'],
    'wald_aut.zip': [ROOTDIR / 'data/forest', 'https://docs.umweltbundesamt.at/s/m46gNNrsMBTgCt9/download/2_FOR_HRL_forest_2015.zip'],
    'Erhaltenswerte_gebaeude.zip': [ROOTDIR / 'data/landuse', 'https://sdi.noe.gv.at/at.gv.noe.geoserver/OGD/wfs?request=GetFeature&version=1.1.0&typeName=OGD:RRU_WI_GEB&srsName=EPSG:31259&outputFormat=shape-zip&format_options=CHARSET:UTF-8'],
    'LuftraumAT.kmz': [ROOTDIR / 'data/airspace', 'https://www.austrocontrol.at/jart/prj3/ac/data/dokumente/20231228LuftraumAT_2023-10-31_0810655.kmz'],
    # find latest download-link at: https://www.austrocontrol.at/piloten/vor_dem_flug/aim_produkte/luftraumstruktur
    'alpsconvention.zip': [ROOTDIR / 'data/alps', 'https://www.atlas.alpconv.org/geoserver/ows?service=WFS&version=1.0.0&request=GetFeature&typename=geonode%3AAlpine_Convention_Perimeter_2018_v2&outputFormat=SHAPE-ZIP&srs=EPSG%3A3034&format_options=charset%3AUTF-8'],
    f'wdpa_{country}.zip': [ROOTDIR / 'data/schutzgebiete', f'https://d1gam3xoknrgr2.cloudfront.net/current/WDPA_WDOECM_{wdpa_date}_Public_{country}_shp.zip'],
    # find latest download link at: https://www.protectedplanet.net/country/AUT
    'Regionale_Raumordnung.zip': [ROOTDIR / 'data/preservation', 'https://sdi.noe.gv.at/at.gv.noe.geoserver/OGD/wfs?request=GetFeature&version=1.1.0&typeName=OGD:RRU_REGROP_ELT&srsName=EPSG:31259&outputFormat=shape-zip&format_options=CHARSET:UTF-8']
}

# %% download global wind atlas
if not os.path.exists(ROOTDIR / 'data/gwa3'):
    os.mkdir(ROOTDIR / 'data/gwa3')

for c in [country]:
    for l in layer:
        for h in height:
            fname = f'{c}_{l}_{h}.tif'
            download_file(f'{url_gwa}/{c}/{l}/{h}', ROOTDIR / 'data/gwa3' / fname,
                          proxy=proxy, proxy_user=proxy_user, proxy_pass=proxy_pass)
            logging.info(f'Download of {fname} complete')
    for g in ground:
        fname = f'{c}_{g}.tif'
        download_file(f'{url_gwa}/{c}/{g}', ROOTDIR / 'data/gwa3' / fname,
                      proxy=proxy, proxy_user=proxy_user, proxy_pass=proxy_pass)
        logging.info(f'Download of {fname} complete')

# %% download data_dict data
for file, addr in data_dict.items():
    if not os.path.exists(addr[0]):
        os.mkdir(addr[0])
    os.chdir(addr[0])
    dnld = download_file(addr[1], addr[0] / file,
                         proxy=proxy, proxy_user=proxy_user, proxy_pass=proxy_pass)
    logging.info(f'Download of {file} complete')
    if dnld:
        if (file[-3:] == 'zip') | (file[-3:] == 'kmz'):
            with zipfile.ZipFile(addr[0] / file) as zip_ref:
                zipinfos = zip_ref.infolist()
                zipext = [zifo.filename[-3:] for zifo in zipinfos]
                if 'shp' in zipext:
                    for zipinf in zipinfos:
                        zipinf.filename = f'{file[0:-3]}{zipinf.filename[-3:]}'
                        zip_ref.extract(zipinf)
                else:
                    zip_ref.extractall(addr[0])
        n = 0
        for fname in os.listdir(addr[0]):
            if fname.endswith('.zip'):
                with zipfile.ZipFile(addr[0] / fname) as nested_zip:
                    nzipinfos = nested_zip.infolist()
                    nzipext = [nzifo.filename[-3:] for nzifo in nzipinfos]
                    if 'shp' in nzipext:
                        for nzipinf in nzipinfos:
                            nzipinf.filename = f'{nzipinf.filename[0:-4]}_{n}.{nzipinf.filename[-3:]}'
                            nested_zip.extract(nzipinf)
                n += 1

# %% settings for downloading and processing wind turbine data from IG Windkraft
igwurl = 'https://www.igwindkraft.at/src_project/external/maps/generated/gmaps_daten.js'

streptyp = {
    'E-40': 'E40',
    'E40/5.40': 'E40 5.40',
    'E40 5.4': 'E40 5.40',
    'E66 18.7': 'E66 18.70',
    'E66/18.70': 'E66 18.70',
    'E66.18': 'E66 18.70',
    'E66 20.7': 'E66 20.70',
    'E70/E4': 'E70 E4',
    'E70/20.71': 'E70 E4',
    'E70': 'E70 E4',
    'E-101': 'E101',
    'E 101': 'E101',
    'E115/3.000': 'E115',
    '3.XM': '3XM',
    'V126/3450': 'V126',
}

strepher = {
    'ENERCON': 'Enercon',
    'DeWind': 'Dewind',
}

# %% retrieve and process turbine data

if not os.path.exists(ROOTDIR / 'data/AT_turbines'):
    os.mkdir(ROOTDIR / 'data/AT_turbines')
download_file(igwurl, ROOTDIR / 'data/AT_turbines/igwind.js',
              proxy=proxy, proxy_user=proxy_user, proxy_pass=proxy_pass)

with open(ROOTDIR / 'data/AT_turbines/igwind.js', 'rt') as f:
    with open(ROOTDIR / 'data/AT_turbines/turbines.json', 'wt') as g:
        for line in f:
            g.write(line.replace('var officeLayer = ', ''))
f.close()
g.close()

with open(ROOTDIR / 'data/AT_turbines/turbines.json', 'rt') as k:
    turbson = json.load(k)
k.close()

tlst = []
for i in range(0, len(turbson[1]['places'])):
    tlst.append(turbson[1]['places'][i]['data'])

igw = pd.DataFrame(tlst, columns=['Name', 'Betreiber1', 'Betreiber2', 'n_Anlagen', 'kW', 'Type', 'Jahr', 'x', 'lat',
                                  'lon', 'url', 'Hersteller', 'Nabenhöhe', 'Rotordurchmesser'])

igw['Type'] = igw['Type'].replace(streptyp)
igw['Hersteller'] = igw['Hersteller'].replace(strepher)

# clean Types
igw.loc[(igw['Type'] == 'E40') & (igw['kW'] == 500), 'Type'] = 'E40 5.40'
igw.loc[(igw['Type'] == 'E40') & (igw['kW'] == 600), 'Type'] = 'E40 6.44'
igw.loc[(igw['Type'] == 'E66') & (igw['kW'] == 1800), 'Type'] = 'E66 18.70'
igw.loc[(igw['Type'] == 'E82') & (igw['kW'] == 2300), 'Type'] = 'E82 E2'
igw.loc[(igw['Type'] == 'E115') & (igw['kW'] == 3200), 'Type'] = 'E115 E2'
igw.loc[(igw['Type'] == 'M114') & (igw['kW'] == 3170), 'Type'] = '3.2M114'

# Add detail for Oberwaltersdorf -
# source: https://www.ris.bka.gv.at/Dokumente/Bvwg/BVWGT_20150313_W102_2008321_1_00/BVWGT_20150313_W102_2008321_1_00.html
igw.loc[igw['Name'].str.contains('Oberwaltersdorf'), 'Type'] = 'V112'
igw.loc[igw['Name'].str.contains('Oberwaltersdorf'), 'Nabenhöhe'] = '140'
igw.loc[igw['Name'].str.contains('Oberwaltersdorf'), 'Rotordurchmesser'] = '112'

# Add detail for Pretul -
# source: https://www.bundesforste.at/fileadmin/erneuerbare_energie/Folder_Windpark-Pretul_FINAL_screen.pdf
igw.loc[igw['Name'].str.contains('Pretul'), 'Type'] = 'E82 E4'
igw.loc[igw['Name'].str.contains('Pretul'), 'Betreiber1'] = 'Österreichische Bundesforste'

igw.loc[igw['Nabenhöhe'] == '', 'Nabenhöhe'] = np.nan
igw['Nabenhöhe'] = igw['Nabenhöhe'].astype('float')

igw.loc[igw['Rotordurchmesser'] == '', 'Rotordurchmesser'] = np.nan
igw['Rotordurchmesser'] = igw['Rotordurchmesser'].astype('float')

igw.to_csv(ROOTDIR / 'data/AT_turbines/igwturbines.csv', sep=';', decimal=',', encoding='utf8')
tmod = igw[['Hersteller', 'Type']].drop_duplicates().sort_values(['Hersteller', 'Type'])
logging.info(f'Download of wind turbine data complete')
logging.info('Data download completed. Please add required additional data manually.')
logging.info('Required data: tourism/Tabelle 30 GEH 2018.xlsx requested from Statistik Austria; '
             'Important bird areas EUROPE_KBA.shp request from Bird Life')
