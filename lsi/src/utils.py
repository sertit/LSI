# -*- coding: utf-8 -*-
# (LSI) © 2024 by Sertit is licensed under Attribution-NonCommercial-NoDerivatives 4.0 International.
# To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-nd/4.0/
""" Utils """

from typing import Optional, Union # , Any, Callable

import geopandas as gpd
import numpy as np
import rasterio
import xarray as xr
#from sertit import AnyPath, geometry, rasters, rasters_rio, unistra, vectors
from sertit import rasters
from sertit.types import AnyPathType # ,AnyPathStrType
from whitebox import WhiteboxTools
from enum import Enum #, unique
from rasterio.enums import Resampling

def classify_raster(raster, raster_steps, raster_classes):
    """ """
    # Conditions
    conds = (
        [raster <= raster_steps[1]]
        + [
            (raster > raster_steps[i]) & (raster <= raster_steps[i + 1])
            for i in range(1, len(raster_steps) - 1)
        ]
        + [raster > raster_steps[-1]]
    )

    # Create classified array
    class_arr = np.select(
        conds, raster_classes.keys(), default=1  # Minimum class by default
    )
    return class_arr


def xr_to_gdf(
    x_array,
    crs,
    column_name: str = None,
    column_rename: str = None,
):
    df = x_array.to_dataframe().reset_index()
    if column_name is not None and column_rename is not None:
        df.rename(columns={column_name: column_rename}, inplace=True)
    df_geometry = gpd.points_from_xy(df.x, df.y)
    return gpd.GeoDataFrame(df, crs=crs, geometry=df_geometry)

def np_to_xr(raster, reference_raster, crs, band_name ="Value"):
    if len(raster.shape) == 2:
        raster = np.broadcast_to(raster, (1, raster.shape[0], raster.shape[1]))
        #print(raster.shape)
    raster = xr.DataArray(raster, coords = reference_raster.coords, name=band_name)
    raster = raster.rio.set_spatial_dims(x_dim = "x", y_dim = "y")
#    raster = raster.rio.write_crs(crs, inplace=True)
    return raster


def initialize_whitebox_tools(
    whitebox_tools_path: Optional[AnyPathType] = None,
    is_verbose: bool = False,
    compress_rasters: bool = True,
) -> WhiteboxTools:
    wbt = WhiteboxTools()
    wbt.set_verbose_mode(is_verbose)
    wbt.set_compress_rasters(compress_rasters)
    wbt.set_default_callback(None)
    wbt.start_minimized = True

    if whitebox_tools_path is not None:
        wbt.set_whitebox_dir = str(whitebox_tools_path)

    return wbt

class RoutingAlgorithm(Enum):
    D8 = "d8"
    DINF = "dinf"

    def __str__(self) -> str:
        return self.value


# Extracted from hillshade function from sertit package
PATH_ARR_DS = Union[str, tuple, rasterio.DatasetReader]

def aspect(ds: PATH_ARR_DS, proj_crs):
    DEG_2_RAD = np.pi / 180
    array = ds
    # Squeeze if needed
    expand = False
    if len(array.shape) == 3 and array.shape[0] == 1:
        array = np.squeeze(array)
        expand = True
    # Compute slope and aspect
    dx, dy = np.gradient(array, *array.rio.resolution())
    x2_y2 = dx**2 + dy**2
    aspect = np.arctan2(dx, dy)
    # from numpy to xarray
    aspect = np_to_xr(aspect, ds, proj_crs)
    # collocate
    aspect = rasters.collocate(ds, aspect, Resampling.bilinear)
    return aspect

def compute_flow_direction(
    input_dtm_path: AnyPathType,
    output_path: AnyPathType,
    wbt_instance: WhiteboxTools,
    routing: RoutingAlgorithm = RoutingAlgorithm.D8,
) -> None:
    # Select routing algorithm
    if routing == RoutingAlgorithm.D8:
        method = wbt_instance.d8_pointer
    elif routing == RoutingAlgorithm.DINF:
        method = wbt_instance.d_inf_pointer
    else:
        raise ValueError(
            f"Unknown routing algorithm: {routing}. "
            "It should be either 'd8' or 'dinf'."
        )

    # Compute flow accumulation
    result = method(
        str(input_dtm_path),
        str(output_path),
    )


def compute_flow_accumulation(
    input_dtm_path: AnyPathType,
    output_path: AnyPathType,
    wbt_instance: WhiteboxTools,
    routing: RoutingAlgorithm = RoutingAlgorithm.D8,
    pointer: bool = False,
) -> None:
    # Select routing algorithm
    if routing == RoutingAlgorithm.D8:
        method = wbt_instance.d8_flow_accumulation
    elif routing == RoutingAlgorithm.DINF:
        method = wbt_instance.d_inf_flow_accumulation
    else:
        raise ValueError(
            f"Unknown routing algorithm: {routing}. "
            "It should be either 'd8' or 'dinf'."
        )

    # Compute flow accumulation
    result = method(
        str(input_dtm_path),
        str(output_path),
        out_type="cells",
        pntr=pointer,
    )
