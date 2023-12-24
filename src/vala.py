# %% imports
import os
import ssl

import urllib3
import certifi
import shutil
import logging
import numpy as np
import pandas as pd
import geopandas as gpd
import itertools
import subprocess
import rioxarray as rxr
import gams.transfer as gt
from xrspatial import proximity
from pathlib import Path
from operator import itemgetter
from scipy.spatial import KDTree
from scipy.interpolate import interp1d
from ast import literal_eval
from shapely.ops import substring
from shapely.geometry import LineString

import xarray as xr


def download_file(url, save_to, proxy=None, proxy_user=None, proxy_pass=None, overwrite=False):
    """
    downloads a file from a specified url to disk
    :param proxy_pass: proxy password
    :param proxy_user: proxy username
    :param proxy: proxy url:port
    :param url: url-string
    :param save_to: destination file name (string)
    :return:
    """
    dld = False
    if (not Path(save_to).is_file()) or (overwrite is True):
        if proxy is not None:
            default_headers = urllib3.make_headers(proxy_basic_auth=f'{proxy_user}:{proxy_pass}')
            http = urllib3.ProxyManager(proxy, proxy_headers=default_headers, ca_certs=certifi.where())
        else:
            http = urllib3.PoolManager(ca_certs=certifi.where())
        try:
            with http.request('GET', url.replace('"', '').replace(' ', ''),
                              preload_content=False) as r, open(save_to, 'wb') as out_file:
                shutil.copyfileobj(r, out_file)
        except:
            try:
                http = urllib3.PoolManager(cert_reqs='CERT_NONE')
                with http.request('GET', url.replace('"', '').replace(' ', ''),
                                  preload_content=False) as r, open(save_to, 'wb') as out_file:
                    shutil.copyfileobj(r, out_file)
            except:
                raise Exception
        dld = True
    return dld


def process_power_curves(turbines, renewables_ninja_curves, open_energy_curves):
    """
    convert power curves from renewables ninja and open energy to common format
    :param turbines: turbine types
    :param renewables_ninja_curves: power curves from renewables.ninja 
    :param open_energy_curves: power curves from open energy
    :param thresh_year: cut-off year for turbines
    :return: power curves in common format, pandas DataFrame
    """
    u_pwrcrv = np.linspace(0.5, 30, num=60)
    powercurves = pd.DataFrame(index=u_pwrcrv)

    missing_turbs = []
    for turbine, props in turbines.items():
        if open_energy_curves['type_string'].str.contains(turbine).any():
            sel = open_energy_curves['type_string'].str.contains(turbine)
            n = open_energy_curves[sel].index[0]
            f_itpl = interp1d(literal_eval(open_energy_curves.loc[n, 'power_curve_wind_speeds']),
                              literal_eval(open_energy_curves.loc[n, 'power_curve_values']) / open_energy_curves.loc[
                                  n, 'nominal_power'], kind='linear', bounds_error=False, fill_value=0)
            powercurves[open_energy_curves.loc[n, 'type_string']] = f_itpl(u_pwrcrv)
        elif turbine in renewables_ninja_curves.columns:
            f_itpl = interp1d(renewables_ninja_curves['speed'], renewables_ninja_curves[turbine], kind='linear')
            powercurves[turbine] = f_itpl(u_pwrcrv)
        else:
            missing_turbs.append(turbine)
    powercurves.index.name = 'speed'
    return powercurves


# %% capacity factor functions
def weibull_probability_density(u_power_curve, k, A):
    """
    Calculates probability density at points in u_power_curve given Weibull parameters in k and A
    :param u_power_curve:
    :param k:
    :param A:
    :return:
    """
    uar = np.asarray(u_power_curve)
    prb = [(k / A * (z / A) ** (k - 1)) * (np.exp(-(z / A) ** k)) for z in uar]
    pdf = xr.concat(prb, dim='wind_speed')
    pdf = pdf.assign_coords({'wind_speed': u_power_curve})
    pdf = pdf.squeeze()
    return pdf


def capacity_factor(pdf, alpha, u_power_curve, p_power_curve, h_turbine, h_reference=100, correction_factor=1):
    """
    calculates wind turbine capacity factors given Weibull probability density pdf, roughness factor alpha, wind turbine
    power curve data in u_power_curve and p_power_curve, turbine height h_turbine and reference height of wind speed
    modelling h_reference
    :param pdf: probability density function from weibull_probability_density()
    :param alpha: roughness coefficient
    :param u_power_curve:
    :param p_power_curve:
    :param h_turbine:
    :param h_reference:
    :return:
    """
    power_curve = xr.DataArray(data=p_power_curve, coords={'wind_speed': u_power_curve})
    u_adjusted = xr.DataArray(data=u_power_curve, coords={'wind_speed': u_power_curve}) @ (h_turbine/h_reference)**alpha
    cap_factor_values = np.trapz(pdf * power_curve, u_adjusted, axis=0)
    cap_factor = alpha.copy()
    cap_factor.values = cap_factor_values * correction_factor
    return cap_factor


def simulate_capacity_factors(powercurves, turbines, rootdir, country, bias_correction):
    # read data
    # read A and k parameters of windspeed Weibull distribution from Austrian wind atlas
    A100 = rxr.open_rasterio(rootdir / f'data/gwa3/{country}_combined-Weibull-A_100.tif')
    A100 = A100.squeeze()
    k100 = rxr.open_rasterio(rootdir / f'data/gwa3/{country}_combined-Weibull-k_100.tif')
    k100 = k100.squeeze()

    # read preprocessed data
    alpha = xr.open_dataarray(rootdir / f'data/preprocessed/gwa_roughness_{country}.nc')
    alpha = alpha.squeeze()
    rho = xr.open_dataarray(rootdir / f'data/preprocessed/gwa_air_density_{country}.nc')
    rho = rho.squeeze()

    u_pwrcrv = powercurves.index.values

    # calculate weibull wind speed probability density
    p = weibull_probability_density(u_pwrcrv, k100, A100)
    logging.info('Weibull probability density computed')

    # fold wind speed probability density and wind turbine power curve
    cf_arr = []
    for turbine_type in turbines.keys():
        cf = capacity_factor(p, alpha, u_pwrcrv, powercurves[turbine_type].values,
                             h_turbine=turbines[turbine_type][1], correction_factor=bias_correction) * rho
        cf_arr.append(cf)
        logging.info(f'Capacity factor for {turbine_type} computed')

    cap_factors = xr.concat(cf_arr, dim='turbine_models', fill_value=np.nan)
    cap_factors = cap_factors.assign_coords({'turbine_models': list(turbines.keys())})
    cap_factors = cap_factors.rio.write_nodata(np.nan)
    cap_factors = cap_factors.rio.reproject('epsg:3416')
    # cap_factors.close()
    return cap_factors


# %% geodata processing
def combine_shapefiles(directory, fname, iterator):
    gdf = gpd.GeoDataFrame()
    for i in iterator:
        shp = gpd.read_file(directory / f'{fname}_{i}.shp')
        gdf = pd.concat([gdf, shp], axis=0, join='outer')
    return gdf


def raster_prox(xrraster, xrtemplate, name):
    if len(xrraster.dims) == 2:
        prx = proximity(xr.where(xrraster > 0, 1, 0), x='x', y='y')
        prx = prx.interp_like(xrtemplate)
        prx = xr.where(xrtemplate.isnull(), np.nan, prx)
        prx = xr.where(prx < 0, 0, prx)
        prx = prx.rio.write_crs(xrtemplate.rio.crs)
    elif len(xrraster.dims) == 3:
        iterdim = [d for d in list(xrraster.dims) if d not in ['x', 'y']]
        prx = []
        for dim in xrraster[iterdim[0]]:
            dimraster = xrraster[xrraster[iterdim[0]] == dim]
            dimraster = dimraster.squeeze()
            dimprx = proximity(xr.where(dimraster > 0, 1, 0), x='x', y='y')
            dimprx = dimprx.interp_like(xrtemplate)
            dimprx = xr.where(xrtemplate[xrtemplate[iterdim[0]] == dim].isnull(), np.nan, dimprx)
            dimprx = xr.where(dimprx < 0, 0, dimprx)
            dimprx = dimprx.rio.write_crs(xrtemplate.rio.crs)
            # dimprx.name = name
            # dimprx.expand_dims({iterdim[0]: dim})
            prx.append(dimprx)

        prx = xr.concat(prx, dim=iterdim[0])
    else:
        raise ValueError('more than 3 dimensional arrays not supported')
    prx.name = name
    return prx



def clip_array_to_polygon(dataarray, polygon, name, interp_array=None, interp_method=None, all_touched=True):
    dataarray = dataarray.rio.reproject(polygon.crs)
    if interp_array is not None:
        if interp_method is None:
            dataarray = dataarray.interp_like(interp_array)
        else:
            dataarray = dataarray.interp_like(interp_array, method=interp_method)
    dataarray_clipped = dataarray.rio.clip(polygon.geometry.values, polygon.crs, all_touched=all_touched)
    dataarray_clipped = dataarray_clipped.squeeze()
    if hasattr(dataarray_clipped, '_FillValue'):
        if dataarray_clipped._FillValue is not None:
            dataarray_clipped = xr.where(dataarray_clipped == dataarray_clipped._FillValue, np.nan, dataarray_clipped)
    dataarray_clipped.name = name
    return dataarray_clipped


def rasterize_shp2xr(rasterarray, clipshape, dummyshape=None, column=None, crs=None, name=None, all_touched=False):
    """
    Clips an xarray DataArray to a geopandas GeoDataFrame with Polygons. If dummyshape is a GeoDataFrame with Polygons,
    raster cells inside the polygons are set to 1 while raster cells outside are set to 0.
    :param rasterarray: an xarray DataArray
    :param clipshape: a GeoDataFrame with Polygon-geometries to which the rasterarray is clipped
    :param dummyshape: a GeoDataFrame with Polygon geometries. Raster cells inside these Polygons are set to 1
    :param crs: a coordinate reference system
    :param all_touched: option from rioxarray clip()-function. If all_touched is True, all raster cells touched by
    Polygon are affected. Otherwise, only raster cells where center point is inside GeoDataFrame-polygons
    :return: an xarray DataArray
    """
    # set everything to the same projection
    if crs is None:
        crs = clipshape.crs
    elif clipshape.crs != crs:
        clipshape = clipshape.to_crs(crs)
    if rasterarray.rio.crs != crs:
        rasterarray = rasterarray.rio.reproject(crs)
    if dummyshape.crs != crs:
        dummyshape = dummyshape.to_crs(crs)

    # set all values initially to 0
    rasterarray.data[~np.isnan(rasterarray.data)] = 0

    if (dummyshape is not None) & (column is None):
        # for dummy, set all raster cells in shapefile to 1
        rasterarray = rasterarray.where(rasterarray.rio.clip(dummyshape.geometry.values, crs, drop=False,
                                                             all_touched=all_touched), 1)
    elif (dummyshape is not None) & (column is not None):
        # for each geometry, set value
        for ix, row in dummyshape.iterrows():
            rasterarray = xr.where(rasterarray.rio.set_crs(crs).rio.clip([row.geometry], crs, drop=False,
                                                                         all_touched=False) == rasterarray,
                                   row[column], rasterarray, keep_attrs=True)
    # clip raster to clipshape
    if rasterarray.rio.crs == crs:
        clippedraster = rasterarray.rio.clip(clipshape.geometry.values, crs, drop=True, all_touched=all_touched)
    else:
        clippedraster = rasterarray.rio.set_crs(crs).rio.clip(clipshape.geometry.values, crs, drop=True,
                                                              all_touched=all_touched)
    clippedraster = clippedraster.squeeze()
    if name is not None:
        clippedraster.name = name
    return clippedraster


def dataset_to_pandas(dataset, digits=4, x='x', y='y', drop_labels=None):

    raw_df_list = []
    for var_name, values in dataset.items():
        if len(values.shape) == 2:
            tmp = values.to_dataframe().dropna()
            ix0 = np.round(tmp.index.get_level_values(y), digits)
            ix1 = np.round(tmp.index.get_level_values(x), digits)
            tmp.index = pd.MultiIndex.from_arrays([ix0, ix1])
            raw_df_list.append(tmp)

        elif len(values.shape) == 3:
            iterdim = [d for d in list(values.dims) if d not in ['x', 'y']]
            tmp = []
            for itm in values[iterdim[0]]:
                subdim = values[values[iterdim[0]] == itm]
                subdim = subdim.squeeze()
                subdim.name = f'{var_name}_{itm.values.item()}'
                tmp = subdim.to_dataframe().dropna()
                ix0 = np.round(tmp.index.get_level_values(y), digits)
                ix1 = np.round(tmp.index.get_level_values(x), digits)
                tmp.index = pd.MultiIndex.from_arrays([ix0, ix1])
                raw_df_list.append(tmp)

    if drop_labels is not None:
        df_list = []
        for j in raw_df_list:
            for lbl in drop_labels:
                if lbl in j.columns:
                    j = j.drop(lbl, axis=1)
            df_list.append(j)
    else:
        df_list = raw_df_list

    df = pd.concat(df_list, axis=1)
    return df

def concat_to_pandas(dataarray_list, digits=4, drop_labels=None):
    """
    converts list of named xr.DataArrays to a pd.DataFrame
    :param dataarray_list: list of named xr.DataArrays. If not named, set name with dataarray.name = 'name'
    :return: pandas DataFrame
    """
    raw_df_list = []
    for i in dataarray_list:
        tmp = i.to_dataframe().dropna()
        ix0 = np.round(tmp.index.get_level_values(0), digits)
        ix1 = np.round(tmp.index.get_level_values(1), digits)
        tmp.index = pd.MultiIndex.from_arrays([ix0, ix1])
        raw_df_list.append(tmp)

    if drop_labels is not None:
        df_list = []
        for j in raw_df_list:
            for lbl in drop_labels:
                if lbl in j.columns:
                    j = j.drop(lbl, axis=1)
            df_list.append(j)
    else:
        df_list = raw_df_list

    df = pd.concat(df_list, axis=1)
    return df


def outside(points, polygons, criterion, splitter):
    """
    returns points outside of polygons.
    Uses column 'criterion' to split up dataset according to df[criterion] == splitter
    :param points: GeoDataFrame
    :param polygons: GeoDataFrame
    :param criterion: string:: Column name
    :param splitter: string:: element of df[criterion]
    :return:
    """
    b = points.loc[points[criterion] == splitter, :]
    s = polygons.loc[polygons[criterion] == splitter, :]
    b = b.drop(['index_right'], axis=1)
    s = s.drop(['index_right'], axis=1)
    remo_builds = gpd.sjoin(b, s, predicate='within', how='left')
    remo_builds = remo_builds.loc[remo_builds['index_right'].isna(), :]
    return remo_builds


def splitlines(lines, num_splits):
    sublines = []
    for i in range(len(lines)):
        for n in range(num_splits):
            sublines.append(substring(lines.iloc[i].geometry,
                                      (n/num_splits)*lines.iloc[i].geometry.length,
                                      ((n+1)/num_splits)*lines.iloc[i].geometry.length))
    linegdf = gpd.GeoDataFrame(geometry=sublines)
    return linegdf


def segments(curve):
    lstlst = []
    curve = curve.reset_index(drop=True)
    for k in range(len(curve)):
        lst = list(map(LineString, zip(curve[k].coords[:-1], curve[k].coords[1:])))
        lstlst.extend(lst)
    return lstlst


def kdnearest(gdfA, gdfB, gdfB_cols=['Place']):
    # resetting the index of gdfA and gdfB here.
    gdfA = gdfA.reset_index(drop=True)
    gdfB = gdfB.reset_index(drop=True)
    # original code snippet from https://gis.stackexchange.com/questions/222315/finding-nearest-point-in-other-geodataframe-using-geopandas/301935#301935
    A = np.concatenate(
        [np.array(geom.coords) for geom in gdfA.geometry.to_list()])
    B = [np.array(geom.coords) for geom in gdfB.geometry.to_list()]
    B_ix = tuple(itertools.chain.from_iterable(
        [itertools.repeat(i, x) for i, x in enumerate(list(map(len, B)))]))
    B = np.concatenate(B)
    kd_tree = KDTree(B)
    dist, idx = kd_tree.query(A, k=1)
    idx = itemgetter(*idx)(B_ix)
    gdf = pd.concat(
        [gdfA, gdfB.loc[idx, gdfB_cols].reset_index(drop=True),
         pd.Series(dist, name='dist')], axis=1)
    return gdf


def calculate_distance(data_array, geo_data_frame, cols=[], crs='epsg:3416', digits=4):
    """
    Calculates the nearest distance from each grid cell center in data_array to each geometry in geo_data_frame.
    :param digits: number of digits to round coordinates to
    :param data_array: a xarray DataArray with (x,y)-coordinate index
    :param geo_data_frame: a GeoDataFrame with the geometries to calculate distances
    :param cols: see kdnearest()-function
    :param crs: a coordinate reference system in which distances are calculate. Should be in meters.
    :return: a xarray DataArray with distances
    """
    data_array = data_array.stack(z=('x', 'y'))
    centers = gpd.GeoDataFrame(geometry=gpd.points_from_xy(data_array[data_array.notnull()].indexes['z'].get_level_values(0),
                                                           data_array[data_array.notnull()].indexes['z'].get_level_values(1),
                                                           crs=data_array.rio.crs))
    centers = centers.to_crs(crs)
    geo_data_frame = geo_data_frame.to_crs(crs)
    if any(geo_data_frame.geom_type.str.contains('Polygon')):
        geo_data_frame = geo_data_frame.boundary.explode()
        geo_data_frame = geo_data_frame.reset_index(drop=True)
        geo_data_frame = gpd.GeoDataFrame(geometry=segments(geo_data_frame))
    distances = kdnearest(centers, geo_data_frame, gdfB_cols=cols)

    distances.index = pd.MultiIndex.from_arrays([distances.geometry.x.round(digits), distances.geometry.y.round(digits)])
    # get x,y-multiindex from data array
    dadf = pd.MultiIndex.from_arrays([data_array.coords.indexes['z'].get_level_values(0).values.round(digits),
                                      data_array.coords.indexes['z'].get_level_values(1).values.round(digits)])
    dadf = pd.DataFrame(np.nan, dadf, ['dist'])
    # broadcast values from distance-dataframe to data array-indexed dataframe
    dadf.loc[distances.index, 'dist'] = distances['dist']
    # assign values from data array-data frame to actual, stacked data array
    data_array.values = dadf['dist'].values
    # unstack data array
    data_array = data_array.unstack()
    data_array = data_array.T
    data_array.name = 'distance'
    return data_array


def count_points_in_raster(points, raster):
    # pts is a 1-D data array of x and y values
    pts = pd.DataFrame(data={'x': points.loc[:, 'geometry'].x, 'y': points.loc[:, 'geometry'].y})
    pts = pts.to_xarray()

    # prepare 2-D data array of x- and y-values
    rstr = raster.copy()
    rstr.name = 'point_count'
    rstr = rstr.stack(z=('x', 'y'))
    rstr.data = rstr.z.data
    rstr = rstr.unstack()

    # count
    tally = rstr.sel(x=pts.x, y=pts.y, method='nearest')
    tally = tally.groupby(tally).count()

    # re-use 2-D template and assign count to raster cells
    rstr.data = xr.zeros_like(rstr)
    rstr = rstr.stack(z=('x', 'y'))
    rstr.loc[dict(z=tally.point_count)] = tally.data
    rstr = rstr.unstack()
    rstr.data = rstr.data.astype(float)
    rstr = xr.where(raster.isnull(), raster, rstr)
    rstr.name = 'point_count'
    rstr = rstr.rio.write_crs(raster.rio.crs)
    rstr = rstr.interp_like(raster)
    rstr = xr.where(rstr < 0, 0, rstr)
    return rstr


# %% optimal wind turbine locations
def sliced_location_optimization(gams_dict, gams_transfer_container, lcoe_array, num_slices, num_turbines, space_px,
                                 min_turbines=1000, max_turbines=10000, lcoe_thresh=85, gdx_out_string='base', read_only=False, axis=0):

    locations = pd.DataFrame()
    lcoe_df = pd.DataFrame(data=lcoe_array.data)
    num_pixels = np.ceil(np.round(lcoe_array.count().data, 0) / (num_slices))
    os.chdir(gams_dict['gdx_output'])
    if axis == 0:
        dim_max = lcoe_array.sizes['x']
    elif axis == 1:
        dim_max = lcoe_array.sizes['y']
    else:
        raise ValueError('axis must be 0 or 1')
    j = 0
    i = 0
    for n in range(0, num_slices):
        size = 0
        while size <= num_pixels and i < dim_max:
            if axis == 0:
                lcoe_slice = lcoe_df.iloc[:, i]
            elif axis == 1:
                lcoe_slice = lcoe_df.iloc[i, :]
            else:
                raise ValueError('Only 2-dimensional arrays allowed')
            lcoe_slice = lcoe_slice.dropna()
            i += 1
            size = size + len(lcoe_slice)

        if axis == 0:
            lcoe_map = lcoe_df.iloc[:, j:i]
        else:
            lcoe_map = lcoe_df.iloc[j:i, :]

        if num_turbines == 'auto':
            nt = np.ceil((lcoe_map < lcoe_thresh).sum().sum() / 1000) * 1000
            ntmax = np.floor(len(lcoe_map.index) * len(lcoe_map.columns) / space_px**3 * 0.5 / 1000) * 1000
            nturb = np.min([nt, ntmax, max_turbines])
            nturb = np.max([nturb, min_turbines])
        else:
            nturb = num_turbines
        print(f'Number of turbines: {nturb}')

        gdx_out = gams_dict['gdx_output'] / f'locations_{gdx_out_string}_{n}.gdx'
        if not read_only:
            for sym in ['l', 'b', 'i', 'j', 'lcoe', 'num_turbines']:
                if gams_transfer_container.hasSymbols(sym):
                    gams_transfer_container.removeSymbols(sym)
            laenge = gams_transfer_container.addSet('l', records=list(lcoe_map.index), description='laenge')
            breite = gams_transfer_container.addSet('b', records=list(lcoe_map.columns), description='breite')
            abstd_l = gams_transfer_container.addSet('i', records=list(range(0, space_px)), description='abstand in laenge')
            abstd_b = gams_transfer_container.addSet('j', records=list(range(0, space_px)), description='abstand in  breite')
            gams_transfer_container.addParameter('lcoe', domain=[laenge, breite],
                                                                  records=lcoe_map.stack().reset_index())
            gams_transfer_container.addParameter('num_turbines', domain=[], records=nturb)
            gams_transfer_container.write(str(gams_dict['gdx_input']))
            # run optimization
            gms_exe_dir = gams_dict['gams_exe']
            gms_model = gams_dict['gams_model']
            subprocess.run(f'{gms_exe_dir}\\gams {gms_model} gdx={gdx_out} lo=3 o=nul')

        results = gt.Container()
        results.read(str(gdx_out), 'build')
        locs = results.data['build'].records[['l', 'b', 'level']]
        locs = locs.loc[locs['level'] > 0]
        locations = pd.concat([locations, locs])  #locations.append(locs)
        # print(f'j: {j}, i: {i}. Using rows {j} to {i} with pixel {num_pixels} of size of {size}')
        j = i
    locations['l'] = locations['l'].astype('int')
    locations['b'] = locations['b'].astype('int')
    return locations


def locations_to_gdf(criterion_array, locations, energy_array=None, power_array=None):

    lcoe_gdf = criterion_array[locations['l'], locations['b']]
    lcoe_gdf = gpd.GeoDataFrame(data=lcoe_gdf.data.diagonal(), columns=['LCOE'],
                                geometry=gpd.points_from_xy(lcoe_gdf.x, lcoe_gdf.y), crs=lcoe_gdf.rio.crs)
    if energy_array is not None:
        energy = energy_array[locations['l'], locations['b']].data.diagonal()
        lcoe_gdf['Energy'] = energy
    if power_array is not None:
        power = power_array[locations['l'], locations['b']].data.diagonal()
        lcoe_gdf['Power'] = power
    lcoe_gdf = lcoe_gdf.sort_values(by='LCOE')
    lcoe_gdf.reset_index(inplace=True, drop=True)
    return lcoe_gdf


# %% LCOE functions
def turbine_overnight_cost(power, hub_height, rotor_diameter, year):
    """
    calculates wind turbine investment cost in EUR per MW based on Rinne et al. (201x)
    :param power: im MW
    :param hub_height: in m
    :param rotor_diameter: in m
    :return: overnight investment cost in EUR per kW
    """ 
    rotor_area = np.pi * (rotor_diameter / 2) ** 2
    spec_power = power * 10**6 / rotor_area
    cost = ((620 * np.log(hub_height)) - (1.68 * spec_power) + (182 * (2016 - year) ** 0.5) - 1005) * 1000
    return cost.astype('float')


def grid_connect_cost(power):
    """
    Calculates grid connection cost according to §54 (3,4) ElWOG https://www.ris.bka.gv.at/GeltendeFassung.wxe?Abfrage=Bundesnormen&Gesetzesnummer=20007045
    :param power: power in kW
    :return:
    """
    cost = 50 * power
    return cost


def discount_factor(discount_rate, period):
    dcf_numerator = 1 - (1 + discount_rate) ** (-period)
    dcf_denominator = 1 - (1 + discount_rate) ** (-1)
    dcf = dcf_numerator / dcf_denominator
    return dcf


def levelized_cost(capacity_factor, overnight_cost, grid_cost, fix_om, var_om, discount_rate, lifetime):
    """
    Calculates wind turbines' levelized cost of electricity in EUR per MWh
    :param capacity_factor: xarray DataArray
    :param overnight_cost: in EUR/MW
    :param grid_cost: xarray DataArray
    :param fix_om: EUR/MW
    :param var_om: EUR/MWh
    :param discount_rate: percent
    :param lifetime: years
    :return:
    """
    npv_energy = capacity_factor * 8760 * discount_factor(discount_rate, lifetime)

    npv_cost = capacity_factor.copy()
    npv_cost = npv_cost.where(npv_cost.isnull(), (var_om * capacity_factor * 8760 + fix_om) * discount_factor(discount_rate, lifetime))
    npv_cost = npv_cost + overnight_cost + grid_cost
    lcoe = npv_cost / npv_energy
    return lcoe


def compute_lcoe(cap_factors, power_curves, turbines, turbine_cost_share, om_fixed, om_variable, discount_rate, lifetime):
    LCOE = []
    onc = []

    for key, value in turbines.items():
        if key in power_curves.columns:
            power = value[0]/1000
            hub_height = value[1]
            rotor_diameter = value[2]
            overnight_cost = np.round(turbine_overnight_cost(power, hub_height, rotor_diameter, value[3]), 0) * turbine_cost_share
            onc.append(overnight_cost)
            cap_fac = cap_factors.sel(turbine_models=key)
            lcoe = levelized_cost(cap_fac, overnight_cost, 0, om_fixed, om_variable, discount_rate, lifetime)
            # grid_connect_cost(power*1000) - set grid connection cost to zero for quasi-lcoe
            # overnight_cost, grid_cost, fix_om, var_om, discount_rate, lifetime)
            LCOE.append(lcoe)
        else:
            pass

    lcoe = xr.concat(LCOE, dim='turbine_models')
    lc_min = lcoe.min(dim='turbine_models', keep_attrs=True)
    lc_min = lc_min.assign_coords({'turbine_models': 'min_lcoe'})
    lcoe = xr.concat([lcoe, lc_min], dim='turbine_models')
    return lcoe
