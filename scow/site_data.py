# %% imports
import os
import json
import logging
import subprocess
import numpy as np
import pandas as pd
import geopandas as gpd

from pathlib import Path
# from cleo.utils import download_file
from datetime import datetime

# TODO: gdal-bin needs to be installed for processing GIP
# sudo apt update
# sudo apt install gdal-bin


def process_gip(gip_path):
    # documented in https://www.gip.gv.at/assets/downloads/GIP_Standardbeschreibung_2.3.2_FINAL.pdf
    # street categories are explained in 2.3.2.6 Abschnittskategorie (GIP.EDGE.EDGECATEGORY)
    # convert to shapefiles
    os.chdir(gip_path)
    subprocess.run(
        """ogr2ogr -f "ESRI shapefile" shp gip_network_ogd.gpkg -sql "select cast(EDGECAT as character(255)) from EDGE_OGD" -overwrite -dialect ogrsql""",
        shell=True)

    # read data
    gip_shp = gip_path / 'shp' / 'EDGE_OGD.shp'

    # abbreviations:
    # A: Autobahn / highway
    # S: Schnellstraße / highway
    # B: Landesstraße B (ehem. Bundesstraße)
    # L: Landesstraße L (ehem. Landeshauptstraße)
    hochrangig = ['A', 'S', 'B', 'L']
    hrng = gpd.GeoDataFrame()

    # total number of rows: 1 532 485
    for n in range(1, 63):
        gip = gpd.read_file(gip_shp, rows=slice((n - 1) * 25000, n * 25000))
        gip = gip.loc[gip['EDGECAT'].isin(hochrangig)]
        hrng = pd.concat([hrng, gip])

    hrng = hrng.dissolve()
    # hrng.to_file(ROOTDIR / 'data/gip/hrng_streets.shp')
    logging.info('High-level road grid preprocessed')
    return hrng


def process_water_bodies(water_bodies_path):
    # process water bodies
    main_water_bodies = ['100 km² Gewässer', '1000 km² Gewässer', '10000 km² Gewässer', '500 km² Gewässer',
                         '4000 km² Gewässer']
    running = gpd.read_file(water_bodies_path / 'Fliessgewaesser.shp')
    running = running.loc[running['GEW_KAT'].isin(main_water_bodies)]
    running = running.dissolve()
    # wbd.to_file(ROOTDIR / 'data/water_bodies/main_running_waters.shp')

    standing = gpd.read_file(water_bodies_path / 'StehendeGewaesser.shp')
    standing = standing.loc[standing['FLAECHEKM2'] >= 0.03125, :]
    standing = standing.dissolve()
    # lks.to_file(ROOTDIR / 'data/water_bodies/main_standing_waters.shp')

    logging.info('Water bodies preprocessed')
    return running, standing


def process_tourism(data_path):

    # read geodata on austrian municipalities
    austria = gpd.read_file(data_path / 'vgd' / 'vgd_oesterreich.shp')
    austria = austria[['GKZ', 'PG', 'BL', 'geometry']].dissolve(by='GKZ')
    austria.reset_index(inplace=True)
    austria['GKZ'] = austria['GKZ'].astype('int')

    # process tourism data
    stays = pd.read_excel(data_path / 'tourism' / 'Tabelle 30 GEH 2018.xlsx',
                          header=[0, 1, 2, 3], skiprows=[0, 5, 6], na_values=['GEH'])
    stays = stays.loc[~stays.iloc[:, 0].isna(), :]
    stays = stays.loc[[not isinstance(e, str) for e in stays.iloc[:, 0]], :]

    col = pd.DataFrame.from_records(data=stays.columns.to_flat_index())
    col = col.replace({"-\n": "", "\n": " ", " ": ""}, regex=True)
    for i in range(0, 4):
        col.loc[col[i].str.contains('Unnamed'), i] = ''

    col = [col.loc[row, :].str.cat(sep=' ').strip() for row in col.index]

    winter = 'WINTERSAISON2017/2018 ÜBERNACHTUNGEN INSGESAMT'
    summer = 'SOMMERSAISON2018 ÜBERNACHTUNGEN INSGESAMT'

    stays.columns = col
    stays['GEM.KENNZIFFER'] = stays['GEM.KENNZIFFER'].astype('int')
    stays[winter] = stays[[winter]].replace('-', 0).astype('float')
    stays[summer] = stays[[summer]].astype('float').fillna(0)
    stays['overnight_stays'] = stays[[winter, summer]].sum(axis=1)

    # merge overnight stays into municipality-geodata
    geostays = austria.merge(stays[['GEM.KENNZIFFER', 'overnight_stays']], left_on='GKZ', right_on='GEM.KENNZIFFER',
                             how='left')
    geostays['overnight_stays'] = geostays['overnight_stays'].fillna(0)
    # geostays.to_file(ROOTDIR / 'data/tourism/overnight_stays.shp')

    logging.info('Touristic overnight stays preprocessed')
    return geostays

"""
def process_windturbines():

    # Settings for downloading and processing wind turbine data from IG Windkraft
    igwurl = 'https://www.igwindkraft.at/src_project/external/maps/generated/gmaps_daten.js'
    data_dir = Path('ROOTDIR') / 'data' / 'AT_turbines'
    igw_file = data_dir / 'igwind.js'
    turbines_file = data_dir / 'turbines.json'

    # Mapping dictionaries
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

    # Download and process turbine data
    download_file(igwurl, igw_file)

    with open(igw_file, 'r') as f:
        data = f.read().replace('var officeLayer = ', '')

    with open(turbines_file, 'w') as f:
        f.write(data)

    turbson = json.load(open(turbines_file, 'r'))
    tlst = [place['data'] for place in turbson[1]['places']]
    igw = pd.DataFrame(tlst, columns=['Name', 'Betreiber1', 'Betreiber2', 'n_Anlagen', 'kW', 'Type', 'Jahr', 'x', 'lat',
                                      'lon', 'url', 'Hersteller', 'Nabenhöhe', 'Rotordurchmesser'])

    igw['Type'].replace(streptyp, inplace=True)
    igw['Hersteller'].replace(strepher, inplace=True)

    # Clean Types
    type_kW_mapping = {
        ('E40', 500): 'E40 5.40',
        ('E40', 600): 'E40 6.44',
        ('E66', 1800): 'E66 18.70',
        ('E82', 2300): 'E82 E2',
        ('E115', 3200): 'E115 E2',
        ('M114', 3170): '3.2M114',
    }

    for (type_, kW), new_type in type_kW_mapping.items():
        igw.loc[(igw['Type'] == type_) & (igw['kW'] == kW), 'Type'] = new_type

    # Add details for specific turbine locations
    location_details = {
        'Oberwaltersdorf': {'Type': 'V112', 'Nabenhöhe': '140', 'Rotordurchmesser': '112'},
        'Pretul': {'Type': 'E82 E4', 'Betreiber1': 'Österreichische Bundesforste'}
    }

    for location, details in location_details.items():
        igw.loc[igw['Name'].str.contains(location), list(details.keys())] = list(details.values())

    # Convert columns to appropriate data types
    igw[['Nabenhöhe', 'Rotordurchmesser']] = igw[['Nabenhöhe', 'Rotordurchmesser']].apply(pd.to_numeric,
                                                                                          errors='coerce')

    # Save processed turbine data to CSV
    igw.to_csv(data_dir / 'igwturbines.csv', sep=';', decimal=',', encoding='utf-8', index=False)
    print('Download of wind turbine data complete')

    # Merge wind turbine location data
    bev = gpd.read_file('ROOTDIR/data/buildings/BAU_2200_BETRIEB_P_0.shp')
    bev = bev.loc[bev['F_CODE'] == 2202].to_crs('epsg:4326')
    bev['geometry'] = bev['geometry'].buffer(20)

    igw = gpd.GeoDataFrame(pd.read_csv(data_dir / 'igwturbines.csv', sep=';', decimal=','),
                           geometry=gpd.points_from_xy(igw['lon'], igw['lat']), crs='epsg:4326')

    turbine_locations = gpd.sjoin(bev, igw, how='left', op='contains')
    turbine_locations['ERSTELLDAT'] = pd.to_datetime(turbine_locations['ERSTELLDAT'])
    turbine_locations['Jahr'].fillna(turbine_locations['ERSTELLDAT'].dt.year, inplace=True)
    turbine_locations['Jahr'] = turbine_locations['Jahr'].astype(int)
    turbine_locations['ERSTELLDAT'] = turbine_locations['ERSTELLDAT'].astype(str)
    turbine_locations = turbine_locations[
        ['NAME', 'Name', 'BETREIBER', 'Betreiber1', 'Betreiber2', 'HOEHE', 'Nabenhöhe',
         'Rotordurchmesser', 'kW', 'Hersteller', 'Type', 'Jahr', 'ERSTELLDAT',
         'n_Anlagen', 'url', 'geometry']]
    turbine_locations.to_file(data_dir / 'turbines.shp')

    print('Turbine locations preprocessed')
"""

def process_protected_areas(path):
    protected_areas = gpd.GeoDataFrame()
    directory = path / "data" / "site" / "schutzgebiete"
    recent_files = {1: None, 2: None, 3: None}

    for file_path in directory.glob('*.shp'):
        parts = file_path.stem.split('_')
        if 'shp-polygons' in parts[-2]:
            date = datetime.strptime(parts[2], '%b%Y')  # Assuming the date part is at index 2
            number = int(parts[-1])  # Assuming the number part is the second last part
            current_date, _ = recent_files[number] if recent_files[number] else (None, None)
            if not current_date or date > current_date:
                recent_files[number] = (date, file_path)

    for key, datepath in recent_files.items():
        # fname_protected_areas = f"WDPA_WDOECM_Mar2024_Public_AUT_shp-polygons_{i}.shp"
        shp = gpd.read_file(datepath[1])
        protected_areas = pd.concat([protected_areas, shp], axis=0, join='outer')

    protected_areas = protected_areas.dissolve()
    return protected_areas

"""
def process_wind_turbine_locations(self):
    igw_url = 'https://www.igwindkraft.at/src_project/external/maps/generated/gmaps_daten.js'
    turbines_dir = self.parent.path / "data" / "site" / f"{self.parent.country}_turbines"
    turbines_dir.mkdir(parents=True, exist_ok=True)

    download_file(igw_url, turbines_dir / 'ig_wind.js')

    with open(turbines_dir / 'ig_wind.js', 'r') as f:
        with open(turbines_dir / 'turbines.json', 'w') as g:
            g.writelines(line.replace("var officeLayer = ", "") for line in f)

    turbine_list = [entry["data"] for entry in json.load((turbines_dir / "turbines.json").open("r"))[1]["places"]]
    igw = pd.DataFrame(turbine_list,
                       columns=['Name', 'Betreiber1', 'Betreiber2', 'n_Anlagen', 'kW', 'Type', 'Jahr', 'x', 'lat',
                                'lon', 'url', 'Hersteller', 'Nabenhöhe', 'Rotordurchmesser'])

    strep_typ = {
        'E-40': 'E40', 'E40/5.40': 'E40 5.40', 'E40 5.4': 'E40 5.40', 'E66 18.7': 'E66 18.70',
        'E66/18.70': 'E66 18.70', 'E66.18': 'E66 18.70', 'E66 20.7': 'E66 20.70', 'E70/E4': 'E70 E4',
        'E70/20.71': 'E70 E4', 'E70': 'E70 E4', 'E-101': 'E101', 'E 101': 'E101', 'E115/3.000': 'E115',
        '3.XM': '3XM', 'V126/3450': 'V126',
    }

    strep_her = {'ENERCON': 'Enercon', 'DeWind': 'Dewind'}

    igw["Type"] = igw["Type"].replace(strep_typ)
    igw["Hersteller"] = igw["Hersteller"].replace(strep_her)

    # clean Types
    igw.loc[(igw['Type'] == 'E40') & (igw['kW'] == 500), 'Type'] = 'E40 5.40'
    igw.loc[(igw['Type'] == 'E40') & (igw['kW'] == 600), 'Type'] = 'E40 6.44'
    igw.loc[(igw['Type'] == 'E66') & (igw['kW'] == 1800), 'Type'] = 'E66 18.70'
    igw.loc[(igw['Type'] == 'E82') & (igw['kW'] == 2300), 'Type'] = 'E82 E2'
    igw.loc[(igw['Type'] == 'E115') & (igw['kW'] == 3200), 'Type'] = 'E115 E2'
    igw.loc[(igw['Type'] == 'M114') & (igw['kW'] == 3170), 'Type'] = '3.2M114'

    # source: https://www.ris.bka.gv.at/Dokumente/Bvwg/BVWGT_20150313_W102_2008321_1_00/BVWGT_20150313_W102_2008321_1_00.html
    oberwaltersdorf = igw['Name'].str.contains('Oberwaltersdorf')
    igw.loc[oberwaltersdorf, ['Type', 'Nabenhöhe', 'Rotordurchmesser']] = ['V112', '140', '112']

    # source: https://www.bundesforste.at/fileadmin/erneuerbare_energie/Folder_Windpark-Pretul_FINAL_screen.pdf
    pretul = igw['Name'].str.contains('Pretul')
    igw.loc[pretul, ['Type', 'Betreiber1']] = ['E82 E4', 'Österreichische Bundesforste']

    igw[['Nabenhöhe', 'Rotordurchmesser']] = igw[['Nabenhöhe', 'Rotordurchmesser']].replace('', np.nan).astype(float)
    igw.to_csv(turbines_dir / 'igwturbines.csv', sep=';', decimal=',', encoding='utf8')
    logging.info(f'Download of wind turbine data complete')

    #  merge wind turbine location data from IGW and BEV
    bev = gpd.read_file(self.parent.path / 'data' / "site" / 'buildings' / 'BAU_2200_BETRIEB_P_0.shp').loc[
        lambda df: df.F_CODE == 2202].to_crs(self.parent.crs)
    bev['points'] = bev['geometry']
    bev['geometry'] = bev['geometry'].buffer(20)

    igw = gpd.GeoDataFrame(pd.read_csv(turbines_dir / 'igwturbines.csv', sep=';', decimal=','),
                           geometry=gpd.points_from_xy(igw['lon'], igw['lat']), crs='epsg:4326').to_crs(self.parent.crs)

    turbine_locations = bev.sjoin(igw, how='left', predicate='contains').reset_index(drop=True)
    turbine_locations['ERSTELLDAT'] = turbine_locations['ERSTELLDAT'].apply(pd.to_datetime)
    turbine_locations['Jahr'] = turbine_locations['Jahr'].fillna(turbine_locations['ERSTELLDAT'].dt.year)
    turbine_locations['Jahr'] = turbine_locations['Jahr'].astype(int)
    turbine_locations['ERSTELLDAT'] = turbine_locations['ERSTELLDAT'].astype(str)
    turbine_locations = turbine_locations[['NAME', 'Name', 'BETREIBER', 'Betreiber1', 'Betreiber2', 'HOEHE',
                                           'Nabenhöhe', 'Rotordurchmesser', 'kW', 'Hersteller', 'Type', 'Jahr',
                                           'ERSTELLDAT', 'n_Anlagen', 'url', 'geometry']]

    turbine_locations.to_file(turbines_dir / 'turbines.shp')
    logging.info('Turbine locations preprocessed')
    return turbine_locations
"""

def process_wind_turbine_locations(self):
    import json
    import logging
    import re
    from html import unescape as html_unescape

    import numpy as np
    import pandas as pd
    import requests
    import geopandas as gpd
    from bs4 import BeautifulSoup

    turbines_dir = self.parent.path / "data" / "site" / f"{self.parent.country}_turbines"
    turbines_dir.mkdir(parents=True, exist_ok=True)

    igw_page_url = "https://www.igwindkraft.at/aktuelles/windrad-karte"

    def _to_float_de(s):
        if s is None:
            return np.nan
        s = str(s).strip()
        if not s:
            return np.nan
        s = s.replace(" ", "")
        m = re.search(r"[-+]?\d+(?:[.,]\d+)?", s)
        if not m:
            return np.nan
        token = m.group(0).replace(",", ".")
        try:
            return float(token)
        except ValueError:
            return np.nan

    def _extract_wheels_json(html: str) -> dict:
        # DOM attempt: any node with a wheels attribute
        soup = BeautifulSoup(html, "html.parser")
        node = soup.select_one("[wheels]")
        if node is not None:
            wheels_raw = node.get("wheels")
            if wheels_raw:
                return json.loads(wheels_raw)

        # Regex fallback
        m = re.search(r"wheels\s*=\s*'(.+?)'", html, flags=re.DOTALL)
        if not m:
            m = re.search(r'wheels\s*=\s*"(.+?)"', html, flags=re.DOTALL)
        if not m:
            raise RuntimeError("Could not locate wheels=... in IGW HTML (windrad-karte).")

        wheels_raw = html_unescape(m.group(1))
        return json.loads(wheels_raw)

    def _parse_typeinfo(typeinfo_html: str):
        """
        Example:
          "Type: Seewind, S 20/110<br>Nabenhöhe: 36, Rotordurchmesser: 20"
        Returns:
          (manufacturer, model, hub_height_m, rotor_diameter_m)
        """
        raw = typeinfo_html or ""
        norm = raw.replace("<br>", "\n").replace("<br/>", "\n").replace("<br />", "\n")
        norm = BeautifulSoup(norm, "html.parser").get_text(separator="\n", strip=True)

        manufacturer = ""
        model = ""
        hub_height = np.nan
        rotor_diameter = np.nan

        m = re.search(r"(?im)^\s*Type\s*:\s*(.+?)\s*(?:,\s*(.+?))?\s*$", norm)
        if m:
            manufacturer = (m.group(1) or "").strip()
            model = (m.group(2) or "").strip()
        else:
            first_line = norm.splitlines()[0].strip() if norm.strip() else ""
            m2 = re.search(r"(?i)type\s*:\s*(.+)", first_line)
            payload = m2.group(1).strip() if m2 else first_line
            if "," in payload:
                manufacturer, model = [p.strip() for p in payload.split(",", 1)]
            else:
                manufacturer = payload.strip()

        mh = re.search(r"(?i)\bNabenhöhe\s*:\s*([0-9]+(?:[.,][0-9]+)?)", norm)
        if mh:
            hub_height = _to_float_de(mh.group(1))

        mr = re.search(r"(?i)\bRotordurchmesser\s*:\s*([0-9]+(?:[.,][0-9]+)?)", norm)
        if mr:
            rotor_diameter = _to_float_de(mr.group(1))

        return manufacturer, model, hub_height, rotor_diameter

    def _parse_facilityinfo(facilityinfo: str):
        """
        Example:
          "1 Anlage(n), je 0.11 MW, errichtet: 1994"
        Returns:
          (n_anlagen, mw_each, built_year)
        """
        if not facilityinfo:
            return np.nan, np.nan, np.nan
        s = str(facilityinfo)

        m_n = re.search(r"(?i)\b(\d+)\s*Anlage\(n\)\b", s)
        n_anlagen = int(m_n.group(1)) if m_n else np.nan

        m_mw = re.search(r"(?i)\bje\s*([0-9]+(?:[.,][0-9]+)?)\s*MW\b", s)
        mw_each = _to_float_de(m_mw.group(1)) if m_mw else np.nan

        m_y = re.search(r"(?i)\berrichtet\s*:\s*(\d{4})\b", s)
        built_year = int(m_y.group(1)) if m_y else np.nan

        return n_anlagen, mw_each, built_year

    # --- scrape current IGW data (windrad-karte) and emit the legacy CSV schema ---

    headers = {
        "User-Agent": (
            "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
            "AppleWebKit/605.1.15 (KHTML, like Gecko) Version/17.0 Safari/605.1.15"
        ),
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
        "Accept-Language": "de-AT,de;q=0.9,en;q=0.8",
    }
    r = requests.get(igw_page_url, headers=headers, timeout=60)
    r.raise_for_status()
    wheels = _extract_wheels_json(r.text)

    rows = []
    for group_id, items in wheels.items():
        for idx, item in enumerate(items):
            # Only keep meaningful keys (drop numeric-string nulls "1".."28")
            d = {k: v for k, v in dict(item).items() if not str(k).isdigit()}

            project = d.get("project", "")
            lon = d.get("x", np.nan)
            lat = d.get("y", np.nan)

            facilityInfo = d.get("facilityInfo", "")
            typeInfo = d.get("typeInfo", "")

            n_anlagen, mw_each, built_year = _parse_facilityinfo(facilityInfo)
            manufacturer, model, hub_h, rotor_d = _parse_typeinfo(typeInfo)

            # mweinzeln is already MW per unit (matches facilityInfo "je ... MW")
            mw_each2 = d.get("mweinzeln", mw_each)
            mw_each2 = _to_float_de(mw_each2)

            kW_each = np.nan
            if np.isfinite(mw_each2):
                kW_each = int(round(mw_each2 * 1000.0))

            # Legacy output schema (same columns as your old code)
            rows.append({
                "Name": project,
                "Betreiber1": "",          # not provided by current dataset
                "Betreiber2": "",          # not provided by current dataset
                "n_Anlagen": n_anlagen,
                "kW": kW_each,
                "Type": model,
                "Jahr": built_year,
                "x": lon,                  # legacy column existed; keep lon here
                "lat": lat,
                "lon": lon,
                "url": "",                 # not provided by current dataset
                "Hersteller": manufacturer,
                "Nabenhöhe": hub_h,
                "Rotordurchmesser": rotor_d,
            })

    igw = pd.DataFrame(
        rows,
        columns=[
            "Name", "Betreiber1", "Betreiber2", "n_Anlagen", "kW", "Type", "Jahr", "x", "lat",
            "lon", "url", "Hersteller", "Nabenhöhe", "Rotordurchmesser",
        ],
    )

    # --- your legacy normalization/cleaning (unchanged) ---

    strep_typ = {
        "E-40": "E40", "E40/5.40": "E40 5.40", "E40 5.4": "E40 5.40", "E66 18.7": "E66 18.70",
        "E66/18.70": "E66 18.70", "E66.18": "E66 18.70", "E66 20.7": "E66 20.70", "E70/E4": "E70 E4",
        "E70/20.71": "E70 E4", "E70": "E70 E4", "E-101": "E101", "E 101": "E101", "E115/3.000": "E115",
        "3.XM": "3XM", "V126/3450": "V126",
    }
    strep_her = {"ENERCON": "Enercon", "DeWind": "Dewind"}

    igw["Type"] = igw["Type"].replace(strep_typ)
    igw["Hersteller"] = igw["Hersteller"].replace(strep_her)

    # clean Types (same rules; depends on kW)
    igw.loc[(igw["Type"] == "E40") & (igw["kW"] == 500), "Type"] = "E40 5.40"
    igw.loc[(igw["Type"] == "E40") & (igw["kW"] == 600), "Type"] = "E40 6.44"
    igw.loc[(igw["Type"] == "E66") & (igw["kW"] == 1800), "Type"] = "E66 18.70"
    igw.loc[(igw["Type"] == "E82") & (igw["kW"] == 2300), "Type"] = "E82 E2"
    igw.loc[(igw["Type"] == "E115") & (igw["kW"] == 3200), "Type"] = "E115 E2"
    igw.loc[(igw["Type"] == "M114") & (igw["kW"] == 3170), "Type"] = "3.2M114"

    # source: https://www.ris.bka.gv.at/Dokumente/Bvwg/BVWGT_20150313_W102_2008321_1_00/BVWGT_20150313_W102_2008321_1_00.html
    oberwaltersdorf = igw["Name"].astype(str).str.contains("Oberwaltersdorf", na=False)
    igw.loc[oberwaltersdorf, ["Type", "Nabenhöhe", "Rotordurchmesser"]] = ["V112", 140, 112]

    # source: https://www.bundesforste.at/fileadmin/erneuerbare_energie/Folder_Windpark-Pretul_FINAL_screen.pdf
    pretul = igw["Name"].astype(str).str.contains("Pretul", na=False)
    igw.loc[pretul, ["Type", "Betreiber1"]] = ["E82 E4", "Österreichische Bundesforste"]

    igw[["Nabenhöhe", "Rotordurchmesser"]] = (
        igw[["Nabenhöhe", "Rotordurchmesser"]].replace("", np.nan).astype(float)
    )

    igw.to_csv(turbines_dir / "igwturbines.csv", sep=";", decimal=",", encoding="utf8", index=False)
    logging.info("Download of wind turbine data complete")

    # --- merge wind turbine location data from IGW and BEV (unchanged except one bugfix) ---

    bev = (
        gpd.read_file(self.parent.path / "data" / "site" / "buildings" / "BAU_2200_BETRIEB_P_0.shp")
        .loc[lambda df: df.F_CODE == 2202]
        .to_crs(self.parent.crs)
    )
    bev["points"] = bev["geometry"]
    bev["geometry"] = bev["geometry"].buffer(20)

    # BUGFIX: use the dataframe we just read for points_from_xy (not a stale name)
    igw_df = pd.read_csv(turbines_dir / "igwturbines.csv", sep=";", decimal=",")
    igw = gpd.GeoDataFrame(
        igw_df,
        geometry=gpd.points_from_xy(igw_df["lon"], igw_df["lat"]),
        crs="epsg:4326",
    ).to_crs(self.parent.crs)

    turbine_locations = bev.sjoin(igw, how="left", predicate="contains").reset_index(drop=True)
    turbine_locations["ERSTELLDAT"] = turbine_locations["ERSTELLDAT"].apply(pd.to_datetime)

    turbine_locations["Jahr"] = turbine_locations["Jahr"].fillna(turbine_locations["ERSTELLDAT"].dt.year)
    turbine_locations["Jahr"] = turbine_locations["Jahr"].astype(int)

    turbine_locations["ERSTELLDAT"] = turbine_locations["ERSTELLDAT"].astype(str)

    turbine_locations = turbine_locations[
        [
            "NAME", "Name", "BETREIBER", "Betreiber1", "Betreiber2", "HOEHE",
            "Nabenhöhe", "Rotordurchmesser", "kW", "Hersteller", "Type", "Jahr",
            "ERSTELLDAT", "n_Anlagen", "url", "geometry",
        ]
    ]

    turbine_locations.to_file(turbines_dir / "turbines.shp")
    logging.info("Turbine locations preprocessed")
    return turbine_locations



# %% data dict for Austrian data
def generate_data_dict(data_path, wdpa_date="Mar2024", country="AUT"):
    data_dict = {
        'bev_bauten.zip':
            [data_path / 'buildings',
             'https://data.bev.gv.at/download/DLM/DLM_20230125/DLM_2000_BAUTEN_20230125.zip'],
        'rnj_power_curve_000-smooth.csv':
            [data_path / 'power_curves',
             'https://raw.githubusercontent.com/renewables-ninja/vwf/master/power_curves/Wind%20Turbine%20Power%20Curves%20%7E%205%20(0.01ms%20with%200.00%20w%20smoother).csv'],
        'oep_wind_turbine_library.zip':
            [data_path / 'power_curves',
             "https://openenergy-platform.org/api/v0/schema/supply/tables/wind_turbine_library/rows/?form=datapackage"],
        'vgd_oesterreich.zip':
            [data_path / 'vgd',
             "https://data.bev.gv.at/download/Verwaltungsgrenzen/shp/20211001/VGD_Oesterreich_gst_20211001.zip"],
        'CLC_2018_AT.zip':
            [data_path / 'clc',
             'https://docs.umweltbundesamt.at/s/beBw8fmwyCMA2ga/download/CLC_2018_AT_clip.zip'],
        'gridkit_europe.zip':
            [data_path / 'grid',
             'https://zenodo.org/record/47317/files/gridkit_euorpe.zip?download=1'],
        'Zonierung_noe.zip':
            [data_path / 'zones',
             # "https://sdi.noe.gv.at/at.gv.noe.geoserver/OGD/wfs?request=GetFeature&version=1.1.0&typeName=OGD:RRU_WIND_ZONEN_P19&srsName=EPSG:31259&outputFormat=shape-zip&format_options=CHARSET:UTF-8"],  # Zoning as of 2024; no longer online
             "https://sdi.noe.gv.at/at.gv.noe.geoserver/OGD/wfs?request=GetFeature&version=1.1.0&typeName=OGD:RRU_WIND_ZONEN_P20_ROG14&srsName=EPSG:31259&outputFormat=shape-zip&format_options=CHARSET:UTF-8"],  # Zoning as of 2024
        'B_gip_network_ogd.zip':
            [data_path / 'gip',
             'https://open.gip.gv.at/ogd/B_gip_network_ogd.zip'],
        'Fliessgewaesser.zip':
            [data_path / 'water_bodies',
             'https://docs.umweltbundesamt.at/s/YkgTDiDs9DPstCJ/download/Routen_v16.zip'],
        'StehendeGewaesser.zip':
            [data_path / 'water_bodies',
             'https://docs.umweltbundesamt.at/s/t4jHoXmrwrsjnea/download/stehendeGewaesser_v16.zip', ],
        'Adressregister.zip':
            [data_path / 'gwr',
             'https://data.bev.gv.at/download/Adressregister/Archiv_Adressregister/Adresse_Relationale_Tabellen_Stichtagsdaten_20220403.zip'],
        'elevation.zip':
            [data_path / 'elevation',
             'https://gis.ktn.gv.at/OGD/Geographie_Planung/ogd-10m-at.zip'],
        'Widmungen_noe.zip':
            [data_path / 'landuse',
             'https://sdi.noe.gv.at/at.gv.noe.geoserver/OGD/wfs?request=GetFeature&version=1.1.0&typeName=OGD:RRU_WI_HUELLE&srsName=EPSG:31259&outputFormat=shape-zip&format_options=CHARSET:UTF-8'],
        'wald_aut.zip':
            [data_path / 'forest',
             'https://docs.umweltbundesamt.at/s/m46gNNrsMBTgCt9/download/2_FOR_HRL_forest_2015.zip'],
        'Erhaltenswerte_gebaeude.zip':
            [data_path / 'landuse',
             'https://sdi.noe.gv.at/at.gv.noe.geoserver/OGD/wfs?request=GetFeature&version=1.1.0&typeName=OGD:RRU_WI_GEB&srsName=EPSG:31259&outputFormat=shape-zip&format_options=CHARSET:UTF-8'],
        # find latest download-link at: https://www.austrocontrol.at/piloten/vor_dem_flug/aim_produkte/luftraumstruktur
        'LuftraumAT.kmz':
            [data_path / 'airspace',
             'https://www.austrocontrol.at/jart/prj3/ac/data/dokumente/20260319LuftraumAT_2026-01-23_080116.kmz'],
        'alpsconvention.zip':
            [data_path / 'alps',
             'https://maps.eurac.edu/geoserver/ows?format_options=charset%3AUTF-8&outputFormat=SHAPE-ZIP&request=GetFeature&service=WFS&srs=EPSG%3A32632&typename=geonode%3Aalpineconvention_2f19551311fb14f3075b85124055c19d&version=1.0.0'],
             # 'https://www.atlas.alpconv.org/geoserver/ows?service=WFS&version=1.0.0&request=GetFeature&typename=geonode%3AAlpine_Convention_Perimeter_2018_v2&outputFormat=SHAPE-ZIP&srs=EPSG%3A3034&format_options=charset%3AUTF-8'],
        # find latest download link at: https://www.protectedplanet.net/country/AUT
        f'wdpa_{country}.zip':
            [data_path / 'schutzgebiete',
             f'https://d1gam3xoknrgr2.cloudfront.net/current/WDPA_WDOECM_{wdpa_date}_Public_{country}_shp.zip'],
        'Regionale_Raumordnung.zip':
            [data_path / 'preservation',
             'https://sdi.noe.gv.at/at.gv.noe.geoserver/OGD/wfs?request=GetFeature&version=1.1.0&typeName=OGD:RRU_REGROP_ELT&srsName=EPSG:31259&outputFormat=shape-zip&format_options=CHARSET:UTF-8']
    }

    return data_dict
