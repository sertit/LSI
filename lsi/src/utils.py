"""

This file is part of LSI.

LSI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

LSI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with LSI. If not, see <https://www.gnu.org/licenses/>.

"""

""" Utils """


from enum import Enum  # , unique
from typing import Optional, Union  # , Any, Callable

import dask.array as da
import geopandas as gpd
import jenkspy
import numpy as np
import rasterio
import xarray as xr
from rasterio.enums import Resampling

# from sertit import AnyPath, geometry, rasters, rasters_rio, unistra, vectors
from sertit import rasters
from sertit.types import AnyPathType  # ,AnyPathStrType
from whitebox import WhiteboxTools


def classify_raster(raster, raster_steps, raster_classes):
    """
    This function allows to calculate a classification of raster based on STEPS and CLASSES
    Inputs:
        raster: Raster in xarray
        raster_steps: List of STEPS
        raster_classes: Dictionary of CLASSES based on the STEPS
    Outputs:
        classified raster in numpy
    """
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


def np_to_xr(raster, reference_raster, crs, band_name="Value"):
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
        # print(raster.shape)
    raster = xr.DataArray(raster, coords=reference_raster.coords, name=band_name)
    raster = raster.rio.set_spatial_dims(x_dim="x", y_dim="y")
    #    raster = raster.rio.write_crs(crs, inplace=True)
    return raster


def initialize_whitebox_tools(
    whitebox_tools_path: Optional[AnyPathType] = None,
    is_verbose: bool = False,
    compress_rasters: bool = True,
) -> WhiteboxTools:
    """
    Function extracted from the Compute_Hand tool from SERTIT.
    """
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
    """
    This function was extracted from the hillshade function availble at the sertit package
    available at: https://sertit-utils.readthedocs.io/en/
    It allows to calculate only the aspect
    """
    array = ds
    # Squeeze if needed
    if len(array.shape) == 3 and array.shape[0] == 1:
        array = np.squeeze(array)
    # Compute slope and aspect
    dx, dy = np.gradient(array, *array.rio.resolution())
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
    """
    Function extracted from the Compute_Hand tool from SERTIT.
    """
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
    method(
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
    """
    Function extracted from the Compute_Hand tool from SERTIT.
    Compute flow accumulation from a DTM.

    Parameters
    ----------
    input_dtm_path : AnyPathType
        Path to input DTM raster file.
    output_path : AnyPathType
        Path to output flow accumulation file.
    wbt_instance : WhiteboxTools
        Instance for Whitebox tools.
    routing : RoutingAlgorithm, optional
        Routing algorithm, by default RoutingAlgorithm.D8.
    apply_ln : bool, optional
        Whether to apply a natural logarithm transform to the flow accumulation
        raster, by default False.

    Raises
    ------
    ValueError
        Unknown routing algorithm. It should be either "D8" or "DINF".
    RuntimeError
        Failed to compute flow accumulation.

    """
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
    method(
        str(input_dtm_path),
        str(output_path),
        out_type="cells",
        pntr=pointer,
    )


def produce_a_reclass_arr(a_xarr, downsample_factor=200):
    """
    Produce reclassified a
    Args:
        a_xarr : a xarray

    Returns:
        xarray of the reclassified a raster
    """

    # jenks breaks computation
    a_xarr_downsampled = a_xarr[:, ::downsample_factor, ::downsample_factor]
    a_xarr_flatten = a_xarr_downsampled.stack(stacked=[...]).values
    a_xarr_finite = a_xarr_flatten[np.isfinite(a_xarr_flatten)]
    nb_class = 5
    breaks = jenkspy.jenks_breaks(a_xarr_finite, nb_class)

    # get max value from the a_xarr
    a_xarr_max = a_xarr.stack(stacked=[...]).values
    a_xarr_max = a_xarr_max[np.isfinite(a_xarr_max)]
    breaks[5] = a_xarr_max.max()

    # -- List conditions and choices
    a_arr = a_xarr.data
    conditions = [
        (a_arr < breaks[1]),
        (a_arr >= breaks[1]) & (a_arr < breaks[2]),
        (a_arr >= breaks[2]) & (a_arr < breaks[3]),
        (a_arr >= breaks[3]) & (a_arr < breaks[4]),
        (a_arr >= breaks[4]) & (a_arr < breaks[5]),
        (a_arr >= breaks[5]),
    ]
    choices = [1, 2, 3, 4, 5, 6]

    a_reclass_arr = da.zeros_like(a_arr, dtype=int)
    for condition, choice in zip(conditions, choices):
        a_reclass_arr = da.where(condition, choice, a_reclass_arr)

    return a_xarr.copy(data=a_reclass_arr)
