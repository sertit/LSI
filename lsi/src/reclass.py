# This file is part of LSI.
# LSI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# LSI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with LSI. If not, see <https://www.gnu.org/licenses/>.
""" Reclass """

from enum import unique

import numpy as np
import xarray as xr
from sertit.misc import ListEnum


@unique
class LandcoverType(ListEnum):
    """
    List of the Landcover type
    """

    CLC = "Corine Land Cover - 2018 (100m)"
    ESAWC = "ESA WorldCover - 2021 (10m)"
    GLC = "Global Land Cover - Copernicus 2019 (100m)"
    ELC = "ESRI Annual Land Cover 2021 (10m)"


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
        conds, raster_classes.keys(), default=1  # Minimum class by default
    )
    return class_arr


def reclass_landcover(landcover, landcover_name):
    """
    This function returns the landcover reclassified from ESA WorldCover / Corine Land Cover to the Final Weights Standard
    for the LSI calculation in the GLOBAL standard.
    Args:
        landcover: landcover xarray
        landcover_name: Name of the Landcover to be reclassified
    Returns:
        landcover xarray reclassified
    """
    if landcover_name == LandcoverType.ESAWC.value:
        # Reclassification of LULC for ESA WorldCover
        landcover_reclass = xr.where(
            landcover == 10, 3, landcover
        )  # Tree cover -> Forest
        landcover_reclass = xr.where(
            landcover_reclass == 20, 3, landcover_reclass
        )  # Shrubland -> Forest
        landcover_reclass = xr.where(
            landcover_reclass == 30, 2, landcover_reclass
        )  # Grassland -> Grassland
        landcover_reclass = xr.where(
            landcover_reclass == 40, 1, landcover_reclass
        )  # Cropland -> Arable land
        landcover_reclass = xr.where(
            landcover_reclass == 50, 5, landcover_reclass
        )  # Built-up -> Urban areas
        landcover_reclass = xr.where(
            landcover_reclass == 60, 4, landcover_reclass
        )  # Bare/Sparse vegetation -> Bare
        landcover_reclass = xr.where(
            landcover_reclass == 70, 997, landcover_reclass
        )  # Snow and ice -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 80, 6, landcover_reclass
        )  # Permanent water bodies -> Water areas
        landcover_reclass = xr.where(
            landcover_reclass == 90, 997, landcover_reclass
        )  # Herbaceous wetland -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 95, 997, landcover_reclass
        )  # Mangroves -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 100, 2, landcover_reclass
        )  # Moss and lichen -> Grassland

    elif landcover_name == LandcoverType.CLC.value:
        # Reclassification of LULC for Corine Land Cover
        landcover_reclass = xr.where(
            landcover == 111, 5, landcover
        )  # Artificial -> Urban areas
        landcover_reclass = xr.where(
            landcover_reclass == 112, 5, landcover_reclass
        )  # Artificial -> Urban areas
        landcover_reclass = xr.where(
            landcover_reclass == 121, 5, landcover_reclass
        )  # Artificial -> Urban areas
        landcover_reclass = xr.where(
            landcover_reclass == 122, 5, landcover_reclass
        )  # Artificial -> Urban areas
        landcover_reclass = xr.where(
            landcover_reclass == 123, 5, landcover_reclass
        )  # Artificial -> Urban areas
        landcover_reclass = xr.where(
            landcover_reclass == 124, 5, landcover_reclass
        )  # Artificial -> Urban areas
        landcover_reclass = xr.where(
            landcover_reclass == 131, 5, landcover_reclass
        )  # Artificial -> Urban areas
        landcover_reclass = xr.where(
            landcover_reclass == 132, 5, landcover_reclass
        )  # Artificial -> Urban areas
        landcover_reclass = xr.where(
            landcover_reclass == 133, 5, landcover_reclass
        )  # Artificial -> Urban areas
        landcover_reclass = xr.where(
            landcover_reclass == 141, 5, landcover_reclass
        )  # Artificial -> Urban areas
        landcover_reclass = xr.where(
            landcover_reclass == 142, 5, landcover_reclass
        )  # Artificial -> Urban areas
        landcover_reclass = xr.where(
            landcover_reclass == 211, 1, landcover_reclass
        )  # Cropland -> Arable Land
        landcover_reclass = xr.where(
            landcover_reclass == 212, 1, landcover_reclass
        )  # Cropland -> Arable Land
        landcover_reclass = xr.where(
            landcover_reclass == 213, 1, landcover_reclass
        )  # Cropland -> Arable Land
        landcover_reclass = xr.where(
            landcover_reclass == 221, 1, landcover_reclass
        )  # Cropland -> Arable Land
        landcover_reclass = xr.where(
            landcover_reclass == 222, 1, landcover_reclass
        )  # Cropland -> Arable Land
        landcover_reclass = xr.where(
            landcover_reclass == 223, 1, landcover_reclass
        )  # Cropland -> Arable Land
        landcover_reclass = xr.where(
            landcover_reclass == 231, 2, landcover_reclass
        )  # Pasture/Meadow -> Grassland
        landcover_reclass = xr.where(
            landcover_reclass == 241, 1, landcover_reclass
        )  # Cropland -> Arable Land
        landcover_reclass = xr.where(
            landcover_reclass == 242, 1, landcover_reclass
        )  # Cropland -> Arable Land
        landcover_reclass = xr.where(
            landcover_reclass == 243, 1, landcover_reclass
        )  # Cropland -> Arable Land
        landcover_reclass = xr.where(
            landcover_reclass == 244, 3, landcover_reclass
        )  # Open Forest -> Forest
        landcover_reclass = xr.where(
            landcover_reclass == 311, 3, landcover_reclass
        )  # Closed Forest -> Forest
        landcover_reclass = xr.where(
            landcover_reclass == 312, 3, landcover_reclass
        )  # Closed Forest -> Forest
        landcover_reclass = xr.where(
            landcover_reclass == 313, 3, landcover_reclass
        )  # Closed Forest -> Forest
        landcover_reclass = xr.where(
            landcover_reclass == 321, 2, landcover_reclass
        )  # Pasture/Meadow -> Grassland
        landcover_reclass = xr.where(
            landcover_reclass == 322, 2, landcover_reclass
        )  # Shrub -> Grassland
        landcover_reclass = xr.where(
            landcover_reclass == 323, 2, landcover_reclass
        )  # Shrub -> Grassland
        landcover_reclass = xr.where(
            landcover_reclass == 324, 3, landcover_reclass
        )  # Open Forest -> Forest
        landcover_reclass = xr.where(
            landcover_reclass == 331, 4, landcover_reclass
        )  # Bare -> Bare
        landcover_reclass = xr.where(
            landcover_reclass == 332, 4, landcover_reclass
        )  # Bare -> Bare
        landcover_reclass = xr.where(
            landcover_reclass == 333, 4, landcover_reclass
        )  # Bare -> Bare
        landcover_reclass = xr.where(
            landcover_reclass == 334, 4, landcover_reclass
        )  # Bare -> Bare
        landcover_reclass = xr.where(
            landcover_reclass == 335, 997, landcover_reclass
        )  # Not applicable -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 411, 997, landcover_reclass
        )  # Not applicable -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 412, 997, landcover_reclass
        )  # Not applicable -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 421, 997, landcover_reclass
        )  # Not applicable -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 422, 997, landcover_reclass
        )  # Not applicable -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 423, 997, landcover_reclass
        )  # Not applicable -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 511, 6, landcover_reclass
        )  # Not applicable -> Water areas
        landcover_reclass = xr.where(
            landcover_reclass == 512, 6, landcover_reclass
        )  # Not applicable -> Water areas
        landcover_reclass = xr.where(
            landcover_reclass == 521, 6, landcover_reclass
        )  # Not applicable -> Water areas
        landcover_reclass = xr.where(
            landcover_reclass == 522, 6, landcover_reclass
        )  # Not applicable -> Water areas
        landcover_reclass = xr.where(
            landcover_reclass == 523, 6, landcover_reclass
        )  # Not applicable -> Water areas
        landcover_reclass = xr.where(
            landcover_reclass == 999, 997, landcover_reclass
        )  # Not applicable -> Not applicable

    elif landcover_name == LandcoverType.GLC.value:
        # Transform to landcover_dbf scale (GLOBAL COPERNICUS LAND COVER)
        landcover_reclass = xr.where(
            landcover == 0, 997, landcover
        )  # No input data available -> NA
        landcover_reclass = xr.where(
            landcover_reclass == 20, 2, landcover_reclass
        )  # Shrubs -> Grassland
        landcover_reclass = xr.where(
            landcover_reclass == 30, 2, landcover_reclass
        )  # Herbaceous vegetation -> Grassland
        landcover_reclass = xr.where(
            landcover_reclass == 40, 1, landcover_reclass
        )  # Cultivated and managed vegetation/agriculture (cropland) -> Arable Land
        landcover_reclass = xr.where(
            landcover_reclass == 50, 5, landcover_reclass
        )  # Urban / built-up -> Urban areas
        landcover_reclass = xr.where(
            landcover_reclass == 60, 4, landcover_reclass
        )  # Bare / sparse vegetation -> Bare
        landcover_reclass = xr.where(
            landcover_reclass == 70, 997, landcover_reclass
        )  # Snow and ice -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 80, 6, landcover_reclass
        )  # Permanent water bodies -> Water areas
        landcover_reclass = xr.where(
            landcover_reclass == 90, 6, landcover_reclass
        )  # Herbaceous wetland -> Water areas
        landcover_reclass = xr.where(
            landcover_reclass == 100, 2, landcover_reclass
        )  # Moss and lichen -> Grassland
        landcover_reclass = xr.where(
            landcover_reclass == 111, 3, landcover_reclass
        )  # Closed forest, evergreen needle leaf -> Forest
        landcover_reclass = xr.where(
            landcover_reclass == 112, 3, landcover_reclass
        )  # Closed forest, evergreen broad leaf -> Forest
        landcover_reclass = xr.where(
            landcover_reclass == 113, 3, landcover_reclass
        )  # Closed forest, deciduous needle leaf -> Forest
        landcover_reclass = xr.where(
            landcover_reclass == 114, 3, landcover_reclass
        )  # Closed forest, deciduous broad leaf -> Forest
        landcover_reclass = xr.where(
            landcover_reclass == 115, 3, landcover_reclass
        )  # Closed forest, mixed -> Forest
        landcover_reclass = xr.where(
            landcover_reclass == 116, 3, landcover_reclass
        )  # Closed forest, unknown -> Forest
        landcover_reclass = xr.where(
            landcover_reclass == 121, 3, landcover_reclass
        )  # Open forest, evergreen needle leaf -> Forest
        landcover_reclass = xr.where(
            landcover_reclass == 122, 3, landcover_reclass
        )  # Open forest, evergreen broad leaf -> Forest
        landcover_reclass = xr.where(
            landcover_reclass == 123, 3, landcover_reclass
        )  # Open forest, deciduous needle leaf -> Forest
        landcover_reclass = xr.where(
            landcover_reclass == 124, 3, landcover_reclass
        )  # Open forest, deciduous broad leaf -> Forest
        landcover_reclass = xr.where(
            landcover_reclass == 125, 3, landcover_reclass
        )  # Open forest, mixed -> Forest
        landcover_reclass = xr.where(
            landcover_reclass == 126, 3, landcover_reclass
        )  # Open forest, unknown -> Forest

    elif landcover_name == LandcoverType.ELC.value:
        # Transform to landcover_dbf scale (ESRI LAND COVER)
        landcover_reclass = xr.where(
            landcover == 1, 6, landcover
        )  # Water -> Water areas
        landcover_reclass = xr.where(
            landcover_reclass == 2, 3, landcover_reclass
        )  # Trees -> Forest
        landcover_reclass = xr.where(
            landcover_reclass == 4, 997, landcover_reclass
        )  # Flooded vegetation -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 5, 1, landcover_reclass
        )  # Crops -> Arable Land
        landcover_reclass = xr.where(
            landcover_reclass == 7, 5, landcover_reclass
        )  # Built Area -> Urban areas
        landcover_reclass = xr.where(
            landcover_reclass == 8, 4, landcover_reclass
        )  # Bare ground -> Bare
        landcover_reclass = xr.where(
            landcover_reclass == 9, 997, landcover_reclass
        )  # Snow/Ice -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 10, 997, landcover_reclass
        )  # Clouds -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 11, 2, landcover_reclass
        )  # Rangeland -> Grassland

    return landcover_reclass


def reclass_landcover_elsus(landcover, proj_crs, landcover_name):
    """
    This function returns the landcover reclassified from Corine Land Cover to the landcover
    classes for ELSUS.
    Args:
        landcover: landcover xarray
        proj_crs: CRS
        landcover_name: Landcover name (ESA WC or CLC)
    Returns:
        landcover xarray reclassified
    """
    if landcover_name == LandcoverType.ESAWC.value:
        # Reclassification of LULC for ESA WorldCover
        landcover_reclass = xr.where(
            landcover == 10, 3, landcover
        )  # Tree cover -> Closed Forest
        landcover_reclass = xr.where(
            landcover_reclass == 20, 4, landcover_reclass
        )  # Shrubland -> Shrub
        landcover_reclass = xr.where(
            landcover_reclass == 30, 5, landcover_reclass
        )  # Grassland -> Pasture/Meadow
        landcover_reclass = xr.where(
            landcover_reclass == 40, 1, landcover_reclass
        )  # Arable land -> Cropland
        landcover_reclass = xr.where(
            landcover_reclass == 50, 7, landcover_reclass
        )  # Built-up -> Artificial
        landcover_reclass = xr.where(
            landcover_reclass == 60, 6, landcover_reclass
        )  # Bare/Sparse vegetation -> Bare
        landcover_reclass = xr.where(
            landcover_reclass == 70, 997, landcover_reclass
        )  # Snow and ice -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 80, 997, landcover_reclass
        )  # Permanent water bodies -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 90, 997, landcover_reclass
        )  # Herbaceous wetland -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 95, 997, landcover_reclass
        )  # Mangroves -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 100, 5, landcover_reclass
        )  # Moss and lichen -> Pasture/Meadow

    elif landcover_name == LandcoverType.CLC.value:
        # Transform to landcover_dbf scale (CORINE LAND COVER)
        landcover_reclass = xr.where(
            landcover == 111, 5, landcover
        )  # Artificial -> Urban areas
        landcover_reclass = xr.where(
            landcover_reclass == 112, 5, landcover_reclass
        )  # Artificial -> Urban areas
        landcover_reclass = xr.where(
            landcover_reclass == 121, 5, landcover_reclass
        )  # Artificial -> Urban areas
        landcover_reclass = xr.where(
            landcover_reclass == 122, 5, landcover_reclass
        )  # Artificial -> Urban areas
        landcover_reclass = xr.where(
            landcover_reclass == 123, 5, landcover_reclass
        )  # Artificial -> Urban areas
        landcover_reclass = xr.where(
            landcover_reclass == 124, 5, landcover_reclass
        )  # Artificial -> Urban areas
        landcover_reclass = xr.where(
            landcover_reclass == 131, 5, landcover_reclass
        )  # Artificial -> Urban areas
        landcover_reclass = xr.where(
            landcover_reclass == 132, 5, landcover_reclass
        )  # Artificial -> Urban areas
        landcover_reclass = xr.where(
            landcover_reclass == 133, 5, landcover_reclass
        )  # Artificial -> Urban areas
        landcover_reclass = xr.where(
            landcover_reclass == 141, 5, landcover_reclass
        )  # Artificial -> Urban areas
        landcover_reclass = xr.where(
            landcover_reclass == 142, 5, landcover_reclass
        )  # Artificial -> Urban areas

        landcover_reclass = xr.where(
            landcover_reclass == 211, 1, landcover_reclass
        )  # Cropland -> Arable Land
        landcover_reclass = xr.where(
            landcover_reclass == 212, 1, landcover_reclass
        )  # Cropland -> Arable Land
        landcover_reclass = xr.where(
            landcover_reclass == 213, 1, landcover_reclass
        )  # Cropland -> Arable Land
        landcover_reclass = xr.where(
            landcover_reclass == 221, 1, landcover_reclass
        )  # Cropland -> Arable Land
        landcover_reclass = xr.where(
            landcover_reclass == 222, 1, landcover_reclass
        )  # Cropland -> Arable Land
        landcover_reclass = xr.where(
            landcover_reclass == 223, 1, landcover_reclass
        )  # Cropland -> Arable Land

        landcover_reclass = xr.where(
            landcover_reclass == 231, 2, landcover_reclass
        )  # Pasture/Meadow -> Grassland
        landcover_reclass = xr.where(
            landcover_reclass == 241, 1, landcover_reclass
        )  # Cropland -> Arable Land
        landcover_reclass = xr.where(
            landcover_reclass == 242, 1, landcover_reclass
        )  # Cropland -> Arable Land
        landcover_reclass = xr.where(
            landcover_reclass == 243, 1, landcover_reclass
        )  # Cropland -> Arable Land
        landcover_reclass = xr.where(
            landcover_reclass == 244, 3, landcover_reclass
        )  # Open Forest -> Forest

        landcover_reclass = xr.where(
            landcover_reclass == 311, 3, landcover_reclass
        )  # Closed Forest -> Forest
        landcover_reclass = xr.where(
            landcover_reclass == 312, 3, landcover_reclass
        )  # Closed Forest -> Forest
        landcover_reclass = xr.where(
            landcover_reclass == 313, 3, landcover_reclass
        )  # Closed Forest -> Forest
        landcover_reclass = xr.where(
            landcover_reclass == 321, 2, landcover_reclass
        )  # Pasture/Meadow -> Grassland
        landcover_reclass = xr.where(
            landcover_reclass == 322, 2, landcover_reclass
        )  # Shrub -> Grassland
        landcover_reclass = xr.where(
            landcover_reclass == 323, 2, landcover_reclass
        )  # Shrub -> Grassland
        landcover_reclass = xr.where(
            landcover_reclass == 324, 3, landcover_reclass
        )  # Open Forest -> Forest

        landcover_reclass = xr.where(
            landcover_reclass == 331, 4, landcover_reclass
        )  # Bare -> Bare
        landcover_reclass = xr.where(
            landcover_reclass == 332, 4, landcover_reclass
        )  # Bare -> Bare
        landcover_reclass = xr.where(
            landcover_reclass == 333, 4, landcover_reclass
        )  # Bare -> Bare
        landcover_reclass = xr.where(
            landcover_reclass == 334, 4, landcover_reclass
        )  # Bare -> Bare

        landcover_reclass = xr.where(
            landcover_reclass == 335, 997, landcover_reclass
        )  # Not applicable -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 411, 997, landcover_reclass
        )  # Not applicable -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 412, 997, landcover_reclass
        )  # Not applicable -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 421, 997, landcover_reclass
        )  # Not applicable -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 422, 997, landcover_reclass
        )  # Not applicable -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 423, 997, landcover_reclass
        )  # Not applicable -> Not applicable

        landcover_reclass = xr.where(
            landcover_reclass == 511, 6, landcover_reclass
        )  # Not applicable -> Water areas
        landcover_reclass = xr.where(
            landcover_reclass == 512, 6, landcover_reclass
        )  # Not applicable -> Water areas
        landcover_reclass = xr.where(
            landcover_reclass == 521, 6, landcover_reclass
        )  # Not applicable -> Water areas
        landcover_reclass = xr.where(
            landcover_reclass == 522, 6, landcover_reclass
        )  # Not applicable -> Water areas
        landcover_reclass = xr.where(
            landcover_reclass == 523, 6, landcover_reclass
        )  # Not applicable -> Water areas
        landcover_reclass = xr.where(
            landcover_reclass == 999, 997, landcover_reclass
        )  # Not applicable -> Not applicable

    elif landcover_name == LandcoverType.GLC.value:
        # Transform to landcover_dbf scale (GLOBAL COPERNICUS LAND COVER)
        landcover_reclass = xr.where(
            landcover == 0, 997, landcover
        )  # No input data available -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 20, 4, landcover_reclass
        )  # Shrubs -> Shrub
        landcover_reclass = xr.where(
            landcover_reclass == 30, 5, landcover_reclass
        )  # Herbaceous vegetation -> Pasture/Meadow
        landcover_reclass = xr.where(
            landcover_reclass == 40, 1, landcover_reclass
        )  # Cultivated and managed vegetation/agriculture (cropland) -> Cropland
        landcover_reclass = xr.where(
            landcover_reclass == 50, 7, landcover_reclass
        )  # Urban / built-up -> Artificial
        landcover_reclass = xr.where(
            landcover_reclass == 60, 6, landcover_reclass
        )  # Bare / sparse vegetation -> Bare
        landcover_reclass = xr.where(
            landcover_reclass == 70, 997, landcover_reclass
        )  # Snow and ice -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 80, 997, landcover_reclass
        )  # Permanent water bodies -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 90, 997, landcover_reclass
        )  # Herbaceous wetland -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 100, 5, landcover_reclass
        )  # Moss and lichen -> Pasture/Meadow
        landcover_reclass = xr.where(
            landcover_reclass == 111, 3, landcover_reclass
        )  # Closed forest, evergreen needle leaf -> Closed Forest
        landcover_reclass = xr.where(
            landcover_reclass == 112, 3, landcover_reclass
        )  # Closed forest, evergreen broad leaf -> Closed Forest
        landcover_reclass = xr.where(
            landcover_reclass == 113, 3, landcover_reclass
        )  # Closed forest, deciduous needle leaf -> Closed Forest
        landcover_reclass = xr.where(
            landcover_reclass == 114, 3, landcover_reclass
        )  # Closed forest, deciduous broad leaf -> Closed Forest
        landcover_reclass = xr.where(
            landcover_reclass == 115, 3, landcover_reclass
        )  # Closed forest, mixed -> Closed Forest
        landcover_reclass = xr.where(
            landcover_reclass == 116, 3, landcover_reclass
        )  # Closed forest, unknown -> Closed Forest
        landcover_reclass = xr.where(
            landcover_reclass == 121, 2, landcover_reclass
        )  # Open forest, evergreen needle leaf -> Open Forest
        landcover_reclass = xr.where(
            landcover_reclass == 122, 2, landcover_reclass
        )  # Open forest, evergreen broad leaf -> Open Forest
        landcover_reclass = xr.where(
            landcover_reclass == 123, 2, landcover_reclass
        )  # Open forest, deciduous needle leaf -> Open Forest
        landcover_reclass = xr.where(
            landcover_reclass == 124, 2, landcover_reclass
        )  # Open forest, deciduous broad leaf -> Open Forest
        landcover_reclass = xr.where(
            landcover_reclass == 125, 2, landcover_reclass
        )  # Open forest, mixed -> Open Forest
        landcover_reclass = xr.where(
            landcover_reclass == 126, 2, landcover_reclass
        )  # Open forest, unknown -> Open Forest

    elif landcover_name == LandcoverType.ELC.value:
        # Transform to landcover_dbf scale (ESRI LAND COVER)
        landcover_reclass = xr.where(
            landcover == 1, 997, landcover
        )  # Water -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 2, 3, landcover_reclass
        )  # Trees -> Closed Forest
        landcover_reclass = xr.where(
            landcover_reclass == 4, 997, landcover_reclass
        )  # Flooded vegetation -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 5, 1, landcover_reclass
        )  # Crops -> Cropland
        landcover_reclass = xr.where(
            landcover_reclass == 7, 7, landcover_reclass
        )  # Built Area -> Artificial
        landcover_reclass = xr.where(
            landcover_reclass == 8, 6, landcover_reclass
        )  # Bare ground -> Bare
        landcover_reclass = xr.where(
            landcover_reclass == 9, 997, landcover_reclass
        )  # Snow/Ice -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 10, 997, landcover_reclass
        )  # Clouds -> Not applicable
        landcover_reclass = xr.where(
            landcover_reclass == 11, 5, landcover_reclass
        )  # Rangeland -> Pasture/Meadow

    landcover_reclass = landcover_reclass.rio.write_crs(proj_crs, inplace=True)
    return landcover_reclass
