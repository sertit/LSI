# This file is part of LSI.
# LSI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# LSI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with LSI. If not, see <https://www.gnu.org/licenses/>.
"""Utils"""

import os
import warnings
from typing import Union

import dask.array as da
import geopandas as gpd

# import jenkspy
import numpy as np
import pandas as pd
import rasterio

# import rasterio as rio
import xarray as xr
from rasterio.enums import Resampling
from rasterstats import zonal_stats

# from rasterio.merge import merge
from sertit import AnyPath, rasters

PATH_ARR_DS = Union[str, tuple, rasterio.DatasetReader]

SIEVE_THRESH = 30  # pixels
MMU = 200  # pixels

BREAKS_GLOBAL = [0, 0.12, 0.15, 0.175, 0.2, 0.35]
BREAKS_EUROPE = [0, 0.08, 0.12, 0.15, 0.175, 0.30]


def xr_to_gdf(
    x_array,
    crs,
    column_name: str = None,
    column_rename: str = None,
):
    """
    This function allows to transform from Xarray to GeoDataFrame
    Inputs:
        x_array: Xarray
        crs: CRS
        column_name: str. Current name of band
        column_rename: str. String to Rename the band
    Output:
        GeoDataFrame
    """
    df = x_array.to_dataframe().reset_index()
    if column_name is not None and column_rename is not None:
        df.rename(columns={column_name: column_rename}, inplace=True)
    df_geometry = gpd.points_from_xy(df.x, df.y)
    return gpd.GeoDataFrame(df, crs=crs, geometry=df_geometry)


def np_to_xr(raster, reference_raster, band_name="Value"):
    """
    This function allows to transform from numpy to xarray
    Inputs:
        raster: Numpy Array
        reference_raster: in Xarray format with available coordinates to
        create the new xarray
        crs: CRS
        band_name: Str, Name of the band to be settled in the Xarray
    """
    if len(raster.shape) == 2:
        raster = np.broadcast_to(raster, (1, raster.shape[0], raster.shape[1]))
    raster = xr.DataArray(raster, coords=reference_raster.coords, name=band_name)
    raster = raster.rio.set_spatial_dims(x_dim="x", y_dim="y")
    return raster


# Extracted from hillshade function from sertit package
def aspect(ds: PATH_ARR_DS):
    """
    This function was extracted from the hillshade function available at the sertit package
    available at: https://sertit-utils.readthedocs.io/en/
    It allows to calculate only the aspect

     Args:
         ds (AnyRasterType): Path to the raster, its dataset, its :code:`xarray` or a tuple containing its array and metadata
         proj_crs (string): Coordinate Reference System

     Returns:
         (np.ma.masked_array, dict): Hillshade and its metadata
    """
    array = ds
    # Squeeze if needed
    if len(array.shape) == 3 and array.shape[0] == 1:
        array = np.squeeze(array)
    # Compute slope and aspect
    dx, dy = np.gradient(array, *array.rio.resolution())
    aspect = np.arctan2(dx, dy)
    # from numpy to xarray
    aspect = np_to_xr(aspect, ds)
    # collocate
    aspect = rasters.collocate(ds, aspect, Resampling.bilinear)
    return aspect


def produce_a_reclass_arr(a_xarr, location):  # downsample_factor=200
    """
    Produce reclassified a
    Args:
        a_xarr : a xarray

    Returns:
        xarray of the reclassified a raster
    """

    breaks = BREAKS_GLOBAL if location == "Global" else BREAKS_EUROPE

    # -- List conditions and choices
    a_arr = a_xarr.data
    conditions = [
        (a_arr <= breaks[1]),
        (a_arr > breaks[1]) & (a_arr <= breaks[2]),
        (a_arr > breaks[2]) & (a_arr <= breaks[3]),
        (a_arr > breaks[3]) & (a_arr <= breaks[4]),
        (a_arr > breaks[4]) & (a_arr <= breaks[5]),
        (a_arr > breaks[5]),
    ]
    choices = [
        1,
        2,
        3,
        4,
        5,
        5,
    ]  # The last one should be 6, but we set 5 to consider that anything above the last break
    # is also Severe

    a_reclass_arr = da.zeros_like(a_arr, dtype=int)
    for condition, choice in zip(conditions, choices):
        a_reclass_arr = da.where(condition, choice, a_reclass_arr)

    return a_xarr.copy(data=a_reclass_arr).astype(np.int32)


def mosaicing(raster_list, output_path, name):
    """
    This function allows to mosaic the rasters by zone
    Args:
        raster_list: list of xarrays
        output_path: Path to be written the new raster
        name: str, name for the raster
    Returns:
        Path for the Mosaic raster in xarray format
    """
    output_path = os.path.join(output_path, AnyPath(str(name) + "_mosaic.tif"))
    mosaic_raster = rasters.merge_gtiff(raster_list, output_path)
    return output_path, mosaic_raster


def raster_postprocess(x_raster: xr.DataArray, resolution) -> gpd.GeoDataFrame:
    """
    This function allows a postprocessing of sieving and MMU in the LSI raster
    Args:
        x_raster: xarrays
        proj_crs: CRS
        output_path: Path to be written the new raster
        name: str, name for the raster
    Returns:
        Path for the Mosaic raster in xarray format
    """
    # Sieve
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        raster_sieved = rasters.sieve(
            x_raster, sieve_thresh=SIEVE_THRESH, connectivity=8
        )

    # Vectorise
    raster_vectorized = rasters.vectorize(raster_sieved)

    return raster_sieved, raster_vectorized


def compute_statistics(gadm_layer, susceptibility_path, location):
    """
    This function allows the zonal statistics and formatting of the geodataframe
    for the LSI statistics based on Deparments level 0, 1 and 2 for the GADM layer
    Args:
        gadm_layer: geodataframe already cropped to the AOI
        raster_path: Path for the LSI raster
    Returns:
        Geodataframe with the statistics data for each of the Levels 0,1 and 2 availables
        in the AOI.
    """
    # Prepare the geodataframe structure
    gadm_df = gadm_layer[["NAME_0", "NAME_1", "NAME_2", "geometry"]]

    # Prepare the three (0, 1, 2) levels of deparments:
    # Level0
    gadm_0 = gadm_df.dissolve(by="NAME_0").reset_index()
    gadm_0["LEVL_CODE"] = 0
    gadm_0 = gadm_0[["NAME_0", "LEVL_CODE", "geometry"]].rename(
        columns={"NAME_0": "NUTS_NAME"}
    )
    # Level1
    gadm_1 = gadm_df.dissolve(by="NAME_1").reset_index()
    gadm_1["LEVL_CODE"] = 1
    gadm_1 = gadm_1[["NAME_1", "LEVL_CODE", "geometry"]].rename(
        columns={"NAME_1": "NUTS_NAME"}
    )
    # Level2
    gadm_2 = gadm_df.dissolve(by="NAME_2").reset_index()
    gadm_2["LEVL_CODE"] = 2
    gadm_2 = gadm_2[["NAME_2", "LEVL_CODE", "geometry"]].rename(
        columns={"NAME_2": "NUTS_NAME"}
    )

    # GADM layer for our AOI
    lsi_stats = pd.concat([gadm_0, gadm_1, gadm_2]).reset_index()
    lsi_stats["FER_LR_ave"] = 0.0
    lsi_stats = lsi_stats[["LEVL_CODE", "NUTS_NAME", "FER_LR_ave", "geometry"]]

    # Compute zonal statistics
    stats = zonal_stats(lsi_stats, susceptibility_path, stats=["mean"])

    breaks = BREAKS_GLOBAL if location == "Global" else BREAKS_EUROPE

    # Add reclassification of Code (1 to 5) and Class (Very low to Severe)
    def reclassify_code(value):
        try:
            if value > breaks[0] and value <= breaks[1]:
                return 1.0
            elif value > breaks[1] and value <= breaks[2]:
                return 2.0
            elif value > breaks[2] and value <= breaks[3]:
                return 3.0
            elif value > breaks[3] and value <= breaks[4]:
                return 4.0
            elif value > breaks[4] and value <= breaks[5]:
                return 5.0
            else:
                return None
        except TypeError:
            return None

    def reclassify_class(value):
        classes = {1: "Very low", 2: "Low", 3: "Moderate", 4: "High", 5: "Severe"}
        try:
            return classes[value]
        except:  # noqa
            return "No data"

    lsi_code = [{"lsi_code": reclassify_code(stat["mean"])} for stat in stats]
    lsi_class = [{"lsi_class": reclassify_class(lsi["lsi_code"])} for lsi in lsi_code]

    # Write average, code and class to GeoDataFrame
    lsi_stats["FER_LR_ave"] = pd.DataFrame(stats)
    lsi_stats["LR_code"] = pd.DataFrame(lsi_code)
    lsi_stats["LR_class"] = pd.DataFrame(lsi_class)

    return lsi_stats
