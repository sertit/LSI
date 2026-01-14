# This file is part of LSI.
# LSI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# LSI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with LSI. If not, see <https://www.gnu.org/licenses/>.
"""Reclass"""

from enum import unique

import numpy as np
import xarray as xr
from sertit.misc import ListEnum

from lsi.src.utils import load_lulc_reclass_table


@unique
class LandcoverType(ListEnum):
    """
    List of the Landcover type
    """

    CLC = "Corine Land Cover - 2018 (100m)"
    ESAWC = "ESA WorldCover - 2021 (10m)"
    GLC = "Global Land Cover - Copernicus 2019 (100m)"
    ELC = "ESRI Annual Land Cover 2021 (10m)"


def reclass_custom_lulc(landcover, reclass_lulc_path: str, location: str):
    """
    Reclassify a custom LULC raster using an external Excel table.

    The Excel table must contain columns 'LULC', 'EUROPE', 'GLOBAL'. The
    'location' argument (\"Europe\" or \"Global\") selects which target column
    to use.
    """
    df = load_lulc_reclass_table(reclass_lulc_path)

    if location == "Europe":
        target_col = "EUROPE"
    elif location == "Global":
        target_col = "GLOBAL"
    else:
        raise ValueError(
            f"Unsupported location '{location}' for custom LULC reclassification. "
            "Expected 'Europe' or 'Global'."
        )

    mapping = dict(zip(df["LULC"], df[target_col]))

    landcover_reclass = xr.apply_ufunc(
        lambda x: mapping.get(x, x),
        landcover,
        vectorize=True,
        dask="parallelized",
        output_dtypes=[landcover.dtype],
    )

    return landcover_reclass


def classify_raster(raster, raster_steps, raster_classes):
    """
    Reclassification (used for Slope, Aspect and Distance to River)

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
        conds,
        raster_classes.keys(),
        default=1,  # Minimum class by default
    )
    return class_arr


def reclass_landcover(landcover, landcover_name):
    """
    Reclassify landcover data from ESA WorldCover / Corine Land Cover / Global Copernicus / ESRI Land Cover
    to the Final Weights Standard for LSI calculation in the GLOBAL standard.

    Args:
        landcover: xarray.DataArray containing landcover data
        landcover_name: Name of the landcover type (from LandcoverType Enum)
    Returns:
        xarray.DataArray with reclassified landcover values
    """

    # Mapping definitions for reclassification
    reclassification_rules = {
        LandcoverType.ESAWC.value: {
            10: 3,
            20: 3,
            30: 2,
            40: 1,
            50: 5,
            60: 4,
            70: 997,
            80: 6,
            90: 997,
            95: 997,
            100: 2,
        },
        LandcoverType.CLC.value: {
            **{k: 5 for k in [111, 112, 121, 122, 123, 124, 131, 132, 133, 141, 142]},
            **{k: 1 for k in [211, 212, 213, 221, 222, 223, 241, 242, 243]},
            231: 2,
            244: 3,
            **{k: 3 for k in [311, 312, 313, 324]},
            **{k: 2 for k in [321, 322, 323]},
            **{k: 4 for k in [331, 332, 333, 334]},
            **{k: 997 for k in [335, 411, 412, 421, 422, 423, 999]},
            **{k: 6 for k in [511, 512, 521, 522, 523]},
        },
        LandcoverType.GLC.value: {
            0: 997,
            20: 2,
            30: 2,
            40: 1,
            50: 5,
            60: 4,
            70: 997,
            80: 6,
            90: 6,
            100: 2,
            **{
                k: 3
                for k in [111, 112, 113, 114, 115, 116, 121, 122, 123, 124, 125, 126]
            },
        },
        LandcoverType.ELC.value: {
            1: 6,
            2: 3,
            4: 997,
            5: 1,
            7: 5,
            8: 4,
            9: 997,
            10: 997,
            11: 2,
        },
    }

    # Get the reclassification mapping for the given landcover_name
    rules = reclassification_rules.get(landcover_name)

    # Apply reclassification using the rules
    landcover_reclass = xr.apply_ufunc(
        lambda x: rules.get(x, x),
        landcover,
        vectorize=True,
        dask="parallelized",
        output_dtypes=[landcover.dtype],
    )

    return landcover_reclass


def reclass_landcover_elsus(landcover, proj_crs, landcover_name):
    """
    Reclassify landcover data for ELSUS analysis based on the specified landcover type.

    Args:
        landcover: xarray.DataArray containing landcover data
        proj_crs: CRS for the dataset
        landcover_name: Landcover type (ESAWC, CLC, GLC, ELC)

    Returns:
        Reclassified landcover xarray.DataArray
    """
    # Define reclassification rules as dictionaries
    reclassification_rules = {
        LandcoverType.ESAWC.value: {
            10: 3,
            20: 4,
            30: 5,
            40: 1,
            50: 7,
            60: 6,
            70: 997,
            80: 997,
            90: 997,
            95: 997,
            100: 5,
        },
        LandcoverType.CLC.value: {
            **{k: 5 for k in [111, 112, 121, 122, 123, 124, 131, 132, 133, 141, 142]},
            **{k: 1 for k in [211, 212, 213, 221, 222, 223, 241, 242, 243]},
            231: 2,
            244: 3,
            **{k: 3 for k in [311, 312, 313, 324]},
            **{k: 2 for k in [321, 322, 323]},
            **{k: 4 for k in [331, 332, 333, 334]},
            **{k: 997 for k in [335, 411, 412, 421, 422, 423]},
            **{k: 6 for k in [511, 512, 521, 522, 523]},
            999: 997,
        },
        LandcoverType.GLC.value: {
            0: 997,
            20: 4,
            30: 5,
            40: 1,
            50: 7,
            60: 6,
            70: 997,
            80: 997,
            90: 997,
            100: 5,
            **{k: 3 for k in range(111, 117)},
            **{k: 2 for k in range(121, 127)},
        },
        LandcoverType.ELC.value: {
            1: 997,
            2: 3,
            4: 997,
            5: 1,
            7: 7,
            8: 6,
            9: 997,
            10: 997,
            11: 5,
        },
    }

    # Select the appropriate rules for the landcover type
    rules = reclassification_rules.get(landcover_name)

    # Apply reclassification using the rules
    landcover_reclass = xr.apply_ufunc(
        lambda x: rules.get(x, x),
        landcover,
        vectorize=True,
        dask="parallelized",
        output_dtypes=[landcover.dtype],
    )

    # Write CRS to the reclassified data
    landcover_reclass = landcover_reclass.rio.write_crs(proj_crs, inplace=True)
    return landcover_reclass
