""" lsi """

import logging
import os
import warnings
from enum import unique

import geopandas as gpd
import numpy as np
import rasterio as rio
import xarray as xr
from rasterio.enums import Resampling
from rasterio.merge import merge
from scipy.ndimage import distance_transform_edt
from sertit import AnyPath, geometry, rasters, vectors
from sertit.misc import ListEnum
from sertit.rasters import FLOAT_NODATA
from sertit.unistra import get_geodatastore, s3_env

from lsi.src.utils import (
    RoutingAlgorithm,
    aspect,
    classify_raster,
    compute_flow_accumulation,
    initialize_whitebox_tools,
    np_to_xr,
    xr_to_gdf,
)

# from lsi import LOGGER_NAME # on the meantime to solve the acces to lsi.py

DEBUG = False
LOGGING_FORMAT = "%(asctime)s - [%(levelname)s] - %(message)s"
LOGGER = logging.getLogger("LSI")  # on the meantime to solve the acces to lsi.py

buffer_n = 0.1


def geodatastore(ftep=False):
    """
    This function returns the root path to the geo data store (DEM, Soil database...).
    Args:
        ftep: If True, the path to the s3 bucket for the ftep platform is returned. Else, get_geodatastore from sertit utils module is called.

    Returns:
        The path to the geo data store.

    """
    if ftep:
        return AnyPath("s3://eo4sdg-data")  # TODO
    else:
        return get_geodatastore()


class DataPath:
    GLOBAL_DIR = None

    @classmethod
    def load_paths(cls, ftep=False):
        cls.GLOBAL_DIR = geodatastore(ftep) / "GLOBAL"
        cls.ESAWC_PATH = (
            cls.GLOBAL_DIR / "ESA_WorldCover" / "2021" / "ESA_WorldCover_10m_2021.vrt"
        )
        cls.CLC_PATH = (
            cls.GLOBAL_DIR
            / "Corine_Land_Cover"
            / "CLC_2018"
            / "clc2018_clc2018_v2018_20_raster100m"
            / "CLC2018_CLC2018_V2018_20.tif"
        )

        #        cls.EUDEM_PATH = cls.GLOBAL_DIR / "EUDEM_v2" / "eudem_dem_3035_europe.tif"
        cls.SRTM30_PATH = cls.GLOBAL_DIR / "SRTM_30m_v4" / "index.vrt"
        cls.COPDEM30_PATH = cls.GLOBAL_DIR / "COPDEM_30m" / "COPDEM_30m.vrt"
        cls.FABDEM_PATH = cls.GLOBAL_DIR / "FABDEM" / "FABDEM.vrt"
        cls.ELSUS_PATH = (
            cls.GLOBAL_DIR / "ELSUSv2" / "ELSUSv2_six_datasets" / "elsus_v2_2.tif"
        )
        cls.WEIGHTS_GLOBAL_PATH = (
            cls.GLOBAL_DIR / "ELSUSv2" / "Weights_Outside"
        )  # This is used as root folder for Weights Outside of Europe
        cls.WEIGHTS_EUROPE_PATH = (
            cls.GLOBAL_DIR / "ELSUSv2" / "DBF_weights"
        )  # This is used as root folder for Weights in Europe
        cls.LITHOLOGY_PATH = cls.GLOBAL_DIR / "LiMW_GIS 2015" / "LiMW_GIS_2015.gdb"
        cls.ELSUS_ZONES_PATH = (
            cls.GLOBAL_DIR
            / "ELSUSv2"
            / "ELSUSv2_six_datasets"
            / "climate_phys_regions.shp"
        )


@unique
class InputParameters(ListEnum):
    """
    List of the input parameters
    """

    AOI_PATH = "aoi_path"
    LOCATION = "location"
    DEM_NAME = "dem_name"
    OTHER_DEM_PATH = "other_dem_path"
    #   LITHOLOGY_PATH = "lithology_path" #not anymore
    LANDCOVER_NAME = "landcover_name"
    EUROPE_METHOD = "europe_method"
    #   WEIGHTS_PATH = "weights_path" #not anymore
    OUTPUT_RESOLUTION = "output_resolution"
    REF_EPSG = "ref_epsg"
    OUTPUT_DIR = "output_path"
    TMP_DIR = "temp_dir"


@unique
class LandcoverType(ListEnum):
    """
    List of the Landcover type
    """

    CLC = "Corine Land Cover - 2018 (100m)"
    ESAWC = "ESA WorldCover - 2021 (10m)"


@unique
class LocationType(ListEnum):
    """
    List of the location
    """

    EUROPE = "Europe"
    GLOBAL = "Global"


@unique
class EUMethodType(ListEnum):
    """
    List of the Methods (similar to location)
    """

    REFINED = "Europe"
    FAST = "Global"


@unique
class DemType(ListEnum):
    """
    List of the DEM
    """

    #    EUDEM = "EUDEM 25m"
    SRTM = "SRTM 30m"
    #    MERIT = "MERIT 5 deg"
    COPDEM_30 = "COPDEM 30m"
    FABDEM = "FABDEM"
    OTHER = "Other"


class LandcoverStructure(ListEnum):
    """
    List of the Landcover type (Coding)
    """

    CLC = "Corine Land Cover - 2018 (100m)"

    ESAWC = "ESA WorldCover - 2021 (10m)"


def reclass_landcover(landcover, landcover_name):
    """
    This function returns the landcover reclassified from ESA WorldCover / Corine Land Cover to the Final Weights Standard
    for the LSI calculation in the GLOBAL standard.
    Args:
        landcover: landcover xarray
        location: [EUROPE, GLOBAL] which indicates whether the reclassification is based on ESAWorldCover or CorineLandCover
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
        )  # Cropland -> Cropland
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

    landcover_reclass = landcover_reclass.rio.write_crs(proj_crs, inplace=True)
    return landcover_reclass


def check_parameters(input_dict: dir) -> None:
    """
     Check if parameters values are valid
    Args:
        input_dict (dict) : dict with parameters values

    Returns:

    """
    # --- Extract parameters ---

    # aoi_path = input_dict.get(InputParameters.AOI_PATH.value)
    location = input_dict.get(InputParameters.LOCATION)
    dem_name = input_dict.get(InputParameters.DEM_NAME.value)
    other_dem_path = input_dict.get(InputParameters.OTHER_DEM_PATH.value)
    # lithology_path = input_dict.get(InputParameters.LITHOLOGY_PATH.value)
    landcover_name = input_dict.get(InputParameters.LANDCOVER_NAME.value)
    # weights_path = input_dict.get(InputParameters.WEIGHTS_PATH.value)
    # epsg_code = input_dict.get(InputParameters.REF_EPSG.value)
    # output_path = input_dict.get(InputParameters.OUTPUT_DIR.value)

    # Check if other_dem_path is needed
    if (dem_name == DemType.OTHER.value) and (other_dem_path is None):
        raise ValueError(f'{"Dem path is needed !"}')
    if (location == LocationType.GLOBAL.value) and (
        landcover_name is LandcoverType.CLC.value
    ):
        raise ValueError(
            f'{"Corine Land Cover can not be used for GLOBAL calculations, only in Europe !"}'
        )
    # if landcover_name != LandcoverType.CLC.value and landcover_name != LandcoverType.CLC.value:
    #     raise ValueError(f"The Landcover was not found with the following name ({landcover_name}) was not identified, please check whether you have written the correct Landcover name!")

    return


# --- GLOBAL LSI method functions


def geology_raster(geology_dbf, litho_shp, dem, aoi, proj_crs, output_path):
    """ """
    LOGGER.info("-- Produce the Geology/Lithology raster for the LSI model")
    if not os.path.exists(os.path.join(output_path, "geology_weight.tif")):
        # if location == LocationType.GLOBAL.value:
        litho_raster = rasters.rasterize(
            path_or_ds=dem, vector=litho_shp, value_field="Rating"
        )

        litho_raster = litho_raster.fillna(997)
        litho_raster = rasters.crop(litho_raster, aoi)

        litho_gdf = xr_to_gdf(
            litho_raster, proj_crs, column_name=litho_raster.name, column_rename="Value"
        )

        # -- JOIN with Geology_dbf
        geology_tif = litho_gdf.merge(geology_dbf, on="Value")
        geology_tif = geology_tif.set_index(["y", "x"]).Weights.to_xarray()
        geology_tif = geology_tif.rio.write_crs(litho_raster.rio.crs)
        geology_tif = rasters.crop(geology_tif, aoi)

        rasters.write(geology_tif, os.path.join(output_path, "geology_weight.tif"))

        return geology_tif
    else:
        return rasters.read(os.path.join(output_path, "geology_weight.tif"))


def slope_raster(slope_dbf, dem_b, aoi, proj_crs, output_path):
    """ """
    LOGGER.info("-- Produce the Slope raster for the LSI model")
    if not os.path.exists(os.path.join(output_path, "slope_weight.tif")):
        slope = rasters.slope(dem_b, in_rad=False)

        # -- Classify
        SLOPE_STEPS = [0, 2, 5, 15, 35, 90]
        SLOPE_CLASSES = {
            1: f"{SLOPE_STEPS[0]} - {SLOPE_STEPS[1]}",
            2: f"{SLOPE_STEPS[1]} - {SLOPE_STEPS[2]}",
            3: f"{SLOPE_STEPS[2]} - {SLOPE_STEPS[3]}",
            4: f"{SLOPE_STEPS[3]} - {SLOPE_STEPS[4]}",
            5: f"{SLOPE_STEPS[4]} - {SLOPE_STEPS[5]}",
            6: f"{SLOPE_STEPS[5]}",
        }
        slope_name = "Value"
        slope_arr = classify_raster(slope, SLOPE_STEPS, SLOPE_CLASSES)
        slope_d = slope.copy(data=slope_arr).astype(np.float32).rename(slope_name)
        slope_d.attrs["long_name"] = slope_name

        slope_gdf = xr_to_gdf(slope_d, proj_crs)
        slope_tif = slope_gdf.merge(slope_dbf, on="Value")
        slope_tif = slope_tif.set_index(["y", "x"]).Weights.to_xarray()
        slope_tif = slope_tif.rio.write_crs(slope_d.rio.crs)
        slope_tif = rasters.crop(slope_tif, aoi)

        rasters.write(slope_tif, os.path.join(output_path, "slope_weight.tif"))

        return slope_tif
    else:
        return rasters.read(os.path.join(output_path, "slope_weight.tif"))


def landcover_raster(
    landuse_dbf,
    lulc_path,
    landcover_name,
    aoi,
    proj_crs,
    output_resolution,
    output_path,
):
    """ """
    LOGGER.info("-- Produce the Land use raster for the LSI model")

    if not os.path.exists(os.path.join(output_path, "landcover_weight.tif")):
        aoi_b = geometry.buffer(aoi, 200)  # buffer (for LandUse + Hydro rasters)
        lulc = rasters.crop(lulc_path, aoi_b)
        # lulc = rasters.collocate(dem, lulc, Resampling.nearest)

        # Reclassification of LULC for LSI calculation
        landcover_reclass = reclass_landcover(lulc, landcover_name)
        landcover_reclass = landcover_reclass.rio.write_crs(lulc.rio.crs, inplace=True)

        landcover_gdf = xr_to_gdf(
            landcover_reclass, lulc.rio.crs, landcover_reclass.name, "Value"
        )
        lulc_tif = landcover_gdf.merge(landuse_dbf, on="Value")
        lulc_tif = lulc_tif.set_index(["y", "x"]).Weights.to_xarray()
        lulc_tif = lulc_tif.rio.write_crs(lulc.rio.crs)
        # Reproject to proj_crs and resolution
        lulc_tif = lulc_tif.rio.reproject(
            proj_crs, resolution=output_resolution, resampling=Resampling.bilinear
        )
        lulc_tif = rasters.crop(lulc_tif, aoi)

        rasters.write(lulc_tif, os.path.join(output_path, "landcover_weight.tif"))
        return lulc_tif
    else:
        return rasters.read(os.path.join(output_path, "landcover_weight.tif"))


def elevation_raster(elevation_dbf, dem_b, aoi, proj_crs, output_path):
    """ """
    LOGGER.info("-- Produce the Elevation raster for the LSI model")
    if not os.path.exists(os.path.join(output_path, "elevation_weight.tif")):

        # Reclassify
        ELEVATION_STEPS = [0, 500, 600, 700, 800, 9000]
        ELEVATION_CLASSES = {
            1: f"{ELEVATION_STEPS[0]} - {ELEVATION_STEPS[1]}",
            2: f"{ELEVATION_STEPS[1]} - {ELEVATION_STEPS[2]}",
            3: f"{ELEVATION_STEPS[2]} - {ELEVATION_STEPS[3]}",
            4: f"{ELEVATION_STEPS[3]} - {ELEVATION_STEPS[4]}",
            5: f"{ELEVATION_STEPS[4]} - {ELEVATION_STEPS[5]}",
            6: f"{ELEVATION_STEPS[5]}",
        }
        elevation_name = "Value"
        elevation_arr = classify_raster(dem_b, ELEVATION_STEPS, ELEVATION_CLASSES)
        elevation_d = (
            dem_b.copy(data=elevation_arr).astype(np.float32).rename(elevation_name)
        )
        elevation_d.attrs["long_name"] = elevation_name

        # JOIN with elevation dbf
        elevation_gdf = xr_to_gdf(elevation_d, proj_crs)
        elevation_tif = elevation_gdf.merge(elevation_dbf, on="Value")
        elevation_tif = elevation_tif.set_index(["y", "x"]).Weights.to_xarray()
        elevation_tif = elevation_tif.rio.write_crs(elevation_d.rio.crs)
        elevation_tif = rasters.crop(elevation_tif, aoi)

        rasters.write(elevation_tif, os.path.join(output_path, "elevation_weight.tif"))

        return elevation_tif
    else:
        return rasters.read(os.path.join(output_path, "elevation_weight.tif"))


def hydro_raster(
    hydro_dbf,
    dem_buff,
    aoi,
    proj_crs,
    ref_raster,
    output_resolution,
    tmp_dir,
    output_path,
):
    """
    Make raster of hydro_weights
    """

    LOGGER.info("-- Produce the Distance to Hydro raster for the LSI model")
    if not os.path.exists(os.path.join(output_path, "hydro_weight.tif")):
        LOGGER.info("-- -- Preprocessing the DEM")
        # -- Pre-process DEM

        # Prepare DEM
        # utm_zone = aoi.estimate_utm_crs()

        # reproject
        dem_b = dem_buff.rio.reproject(
            dst_crs=proj_crs,
            nodata=rasters.FLOAT_NODATA,
            resampling=Resampling.bilinear,
        )
        # no data
        dem_b = xr.where(dem_buff <= -700, FLOAT_NODATA, dem_buff)

        dem_max = dem_b.max()
        dem_min = dem_b.min()

        # LOGGER.info(dem_min, dem_max)

        # dem_max = dem_buff.max()
        # dem_min = dem_buff.min()

        dem_b = dem_b.rio.write_crs(proj_crs)
        # dem_b = dem_b.rio.write_nodata(0, encoded=True)
        # dem_b = dem_b.rio.write_nodata(FLOAT_NODATA, encoded=True)

        # reproject
        dem_b = dem_b.rio.reproject(
            dst_crs=proj_crs,
            # nodata=rasters.FLOAT_NODATA,
            resampling=Resampling.bilinear,
        )
        dem_b = rasters.crop(dem_b, aoi)
        # resolution
        # x_res, y_res = dem_b.rio.resolution()
        # dem_resolution = (abs(x_res) + abs(y_res)) / 2

        # Write DEM in memory
        rasters.write(
            dem_b,
            os.path.join(tmp_dir, "dem_d.tif"),
            compress="deflate",
            predictor=1,
            nodata=FLOAT_NODATA,
            dtype=np.float32,
        )
        # -- Hydro processing

        wbt = initialize_whitebox_tools()

        LOGGER.info("-- -- Preprocessing the DEM: Filling Pits")
        # -- Fill pits
        wbt.fill_single_cell_pits(
            os.path.join(tmp_dir, "dem_d.tif"),
            os.path.join(tmp_dir, "filled_pits.tif"),
        )
        LOGGER.info("-- -- Preprocessing the DEM: Filling Depressions")
        # -- Fill depressions
        wbt.fill_depressions(
            os.path.join(tmp_dir, "filled_pits.tif"),
            os.path.join(tmp_dir, "filled_depressions.tif"),
        )
        # -- Flow accumulation
        LOGGER.info("-- -- Compute Flow Accumulation")
        compute_flow_accumulation(
            os.path.join(tmp_dir, "filled_depressions.tif"),
            os.path.join(tmp_dir, "flow_acc.tif"),
            wbt,
            RoutingAlgorithm.D8,
        )

        flow_acc = rasters.read(os.path.join(tmp_dir, "flow_acc.tif"))
        # Thresholding the flow accumulation
        elevation_threshold = (dem_max - abs(dem_min)).values
        flow_acc_thresh = xr.where(flow_acc > elevation_threshold, 1, 0)
        rasters.write(
            flow_acc_thresh,
            os.path.join(tmp_dir, "flow_acc_thresh.tif"),
            compress="deflate",
            predictor=1,
            dtype=np.float32,
        )
        # Flow_acc raster to polyline
        wbt.raster_to_vector_lines(
            os.path.join(tmp_dir, "flow_acc_thresh.tif"),
            os.path.join(tmp_dir, "flowacc_thresh_polyline.shp"),
        )
        flowacc_thresh_lines = vectors.read(
            os.path.join(tmp_dir, "flowacc_thresh_polyline.shp")
        )
        flowacc_thresh_lines = rasters.rasterize(
            path_or_ds=dem_b, vector=flowacc_thresh_lines, value_field="VALUE"
        )
        rasters.write(
            flowacc_thresh_lines,
            os.path.join(tmp_dir, "flowacc_thresh_lines.tif"),
            compress="deflate",
            predictor=1,
        )
        LOGGER.info("-- -- Compute Distance to rivers")
        # -- Euclidean Distance to River
        with rio.open(os.path.join(tmp_dir, "flowacc_thresh_lines.tif")) as src:
            river_streams = src.read(1)
            nodata = src.nodata

        # Transform raster values
        # Invert the values so that river cells become background
        river_streams_inverted = np.where(river_streams == nodata, 1, 0)

        # Euclidean distance
        euclidean_distance = distance_transform_edt(
            river_streams_inverted
            # , sampling=profile['transform'][0]
        )

        euclidean_distance_xr = np_to_xr(euclidean_distance, dem_b, proj_crs)
        # transform from pixel to meters

        LOGGER.info("-- -- Compute Flow Accumulation")
        euclidean_distance_xr = euclidean_distance_xr * int(output_resolution)
        # -- Reclassify
        ED_STEPS = [0, 100, 200, 300, 400, 20000]  # 5000
        ED_CLASSES = {
            1: f"{ED_STEPS[0]} - {ED_STEPS[1]}",  #
            2: f"{ED_STEPS[1]} - {ED_STEPS[2]}",  #
            3: f"{ED_STEPS[2]} - {ED_STEPS[3]}",  #
            4: f"{ED_STEPS[3]} - {ED_STEPS[4]}",  #
            5: f"{ED_STEPS[4]} - {ED_STEPS[5]}",  #
            6: f"{ED_STEPS[5]}",
        }
        ed_class = classify_raster(euclidean_distance_xr, ED_STEPS, ED_CLASSES)
        ed_class = np_to_xr(ed_class, euclidean_distance_xr, proj_crs)
        ed_reclass = rasters.crop(ed_class, aoi)  # el proj_crs jode todo

        # -- JOIN with Hydro.dbf
        hydro_gdf = xr_to_gdf(ed_reclass, ed_reclass.rio.crs)
        hydro_tif = hydro_gdf.merge(hydro_dbf, on="Value")
        hydro_tif = hydro_tif.set_index(["y", "x"]).Weights.to_xarray()
        hydro_tif = hydro_tif.rio.write_crs(ed_reclass.rio.crs)

        rasters.write(hydro_tif, os.path.join(output_path, "hydro_weight_utm.tif"))

        # LOGGER.info(
        #     "-- -- Transform raster from UTM to LatLon"
        # )
        # # From UTM to LatLon

        # # dst_transform = ref_raster.rio.transform
        # # dst_width = ref_raster.rio.width
        # # dst_height = ref_raster.rio.height
        # # dst_crs = ref_raster.rio.crs

        # # hydro_tif = hydro_tif.rio.reproject(dst_crs = dst_crs,
        # #                                     dst_transfom = dst_transform,
        # #                                     resampling=Resampling.bilinear)
        # # hydro_tif = rasters.collocate(ref_raster, hydro_tif, Resampling.bilinear)

        # with rio.open(os.path.join(output_path, "hydro_weight_utm.tif")) as src:
        #     src_crs = src.crs
        #     src_transform = src.transform
        #     src_width = src.width
        #     src_height = src.height

        #     dst_transform, dst_width, dst_height = calculate_default_transform(
        #         src_crs, proj_crs, src_width, src_height, *src.bounds)
        #     # Define the metadata for the destination file
        #     dst_meta = src.meta.copy()
        #     dst_meta.update({
        #         'crs': proj_crs,
        #         'transform': dst_transform,
        #         'width': dst_width,
        #         'height': dst_height
        #     })

        #     # Create a new file to save the reprojected raster
        #     with rio.open(os.path.join(output_path, "hydro_weight_latlon.tif"), 'w', **dst_meta) as dst:
        #         for i in range(1, src.count + 1):
        #             reproject(
        #                 source=rio.band(src, i),
        #                 destination=rio.band(dst, i),
        #                 src_transform=src_transform,
        #                 src_crs=src_crs,
        #                 dst_transform=dst_transform,
        #                 dst_crs=proj_crs,
        #                 resampling=Resampling.nearest)
        # hydro_tif = rasters.read(os.path.join(output_path, "hydro_weight_latlon.tif"))
        # hydro_tif = rasters.collocate(ref_raster, hydro_tif, Resampling.bilinear)
        # # Write in memory

        hydro_tif = hydro_tif.rio.reproject(
            proj_crs, resolution=output_resolution, resampling=Resampling.bilinear
        )
        rasters.write(hydro_tif, os.path.join(output_path, "hydro_weight.tif"))
        return hydro_tif
    else:
        return rasters.read(os.path.join(output_path, "hydro_weight.tif"))


def aspect_raster(aspect_dbf, dem_b, aoi, proj_crs, output_path):
    """ """
    LOGGER.info("-- Produce the Aspect raster for the LSI model")
    if not os.path.exists(os.path.join(output_path, "aspect_weight.tif")):

        aspect_tif = aspect(dem_b, proj_crs)

        aspect_tif_deg = np.round(aspect_tif * (180 / np.pi)) + 180

        # Taking the maximum Azimuth as the Flat (360 -> -1) for classification purposes
        aspect_tif_deg = xr.where(
            aspect_tif_deg == aspect_tif_deg.max(), -1, aspect_tif_deg
        )

        # -- Classify from degrees [0, 359] Flat/North/Northeast/etc...
        ASPECT_STEPS = [-1, 0, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5]
        ASPECT_CLASSES = {
            1: f"{ASPECT_STEPS[0]} - {ASPECT_STEPS[1]}",  # Flat
            2: f"{ASPECT_STEPS[1]} - {ASPECT_STEPS[2]}",  # North
            3: f"{ASPECT_STEPS[2]} - {ASPECT_STEPS[3]}",  # Northeast
            4: f"{ASPECT_STEPS[3]} - {ASPECT_STEPS[4]}",  # East
            5: f"{ASPECT_STEPS[4]} - {ASPECT_STEPS[5]}",  # Southeast
            6: f"{ASPECT_STEPS[5]} - {ASPECT_STEPS[6]}",  # South
            7: f"{ASPECT_STEPS[6]} - {ASPECT_STEPS[7]}",  # Southwest
            8: f"{ASPECT_STEPS[7]} - {ASPECT_STEPS[8]}",  # West
            9: f"{ASPECT_STEPS[8]} - {ASPECT_STEPS[9]}",  # Northwest
            10: f"{ASPECT_STEPS[9]}",  # North
        }
        aspect_class = classify_raster(aspect_tif_deg, ASPECT_STEPS, ASPECT_CLASSES)

        # -- Classify following the aspect classes from dbf
        # the following codes are not considered = {1: Flat, 2: North, 3: Northeast, 4: East}
        # as they already are within the class defined in Aspect_dbf

        # Transform to Aspect_dbf scale
        aspect_reclass = xr.where(aspect_class == 5, 4, aspect_class)  # Southeast
        aspect_reclass = xr.where(aspect_reclass == 6, 5, aspect_reclass)  # South
        aspect_reclass = xr.where(aspect_reclass == 7, 5, aspect_reclass)  # Southwest
        aspect_reclass = xr.where(aspect_reclass == 8, 3, aspect_reclass)  # West
        aspect_reclass = xr.where(aspect_reclass == 9, 3, aspect_reclass)  # Northeast
        aspect_reclass = xr.where(aspect_reclass == 10, 2, aspect_reclass)  # North

        aspect_reclass_xr = np_to_xr(aspect_reclass, dem_b, proj_crs)
        aspect_gdf = xr_to_gdf(aspect_reclass_xr, proj_crs)

        # JOIN with aspect_dbf
        aspect_tif = aspect_gdf.merge(aspect_dbf, on="Value")
        aspect_tif = aspect_tif.set_index(["y", "x"]).Weights.to_xarray()
        aspect_tif = aspect_tif.rio.write_crs(aspect_reclass_xr.rio.crs)
        aspect_tif = rasters.crop(aspect_tif, aoi)

        rasters.write(aspect_tif, os.path.join(output_path, "aspect_weight.tif"))

        return aspect_tif
    else:
        return rasters.read(os.path.join(output_path, "aspect_weight.tif"))


# --- ELSUS LSI method


def landcover_raster_eu(
    landcover_weight_path,
    landcover_path,
    reference_raster,
    aoi_zone,
    proj_crs,
    zone_id,
    counter,
    final_weight_dbf,
    tmp_dir,
    output_resolution,
    landcover_name,
):
    """
    This function process the landcover to produce tha raster of weights for landcover according to ELSUS.
    Part of the idea of this function is to iterate over the different Zones [0, 6] in the ELSUS method,
    to choose the right weights for each function. As result, the Area of Interest in the project is segmented
    in the different aoi_zone per Zone available and each weights are calculated sepparetly and written as a new
    raster.
    Args:
        landcover_weight_path: Path to ELSUS weights for Landcover
        landcover_path: Path to the Landcover
        aoi_zone: Polygon, the AOI to clip the Landcover (for the ELSUS method it should be the polygon correspondent to
                  the climatic zone)
        proj_crs: CRS
        zone_id: String or Int, the ID of the Zone. i.e: 1
        counter: Int; a counter to number the rasters to be written in memory per Zone
        final_weight_dbf: the GeoDataFrame of the Final Weights for ELSUS.

    Returns:
        Path for the Landcover Final Weights
    """
    LOGGER.info(" --- LSI : Landcover")
    landcover_weight_dbf = gpd.read_file(landcover_weight_path)
    landcover_weight_dbf.loc[len(landcover_weight_dbf)] = [
        997,
        "Not Applicable",
        0.0,
        None,
    ]
    with warnings.catch_warnings():
        warnings.simplefilter(
            "ignore"
        )  # the buffer is being applied in degrees, not meters
        buffer_lc = 200
        aoi_b = geometry.buffer(aoi_zone, buffer_lc)
    try:
        landcover = rasters.crop(landcover_path, aoi_b)
    except ValueError:
        raise ValueError("Your AOI doesn't cover your Landcover layer.")

    # Reclassification based on ELSUS
    landcover = rasters.collocate(reference_raster, landcover, Resampling.nearest)
    landcover_reclass = reclass_landcover_elsus(landcover, proj_crs, landcover_name)
    landcover_reclass = rasters.crop(landcover_reclass, aoi_zone)
    landcover_gdf = xr_to_gdf(
        landcover_reclass, proj_crs, landcover_reclass.name, "Value"
    )

    # JOIN Landcover with Weights gdf
    landcover_tif = landcover_gdf.merge(landcover_weight_dbf, on="Value")
    landcover_tif = landcover_tif.set_index(["y", "x"]).Weight.to_xarray()
    landcover_tif = landcover_tif.rio.write_crs(proj_crs)
    landcover_tif = rasters.crop(landcover_tif, aoi_zone)

    # Calculating final Weights
    zone_class = "Z" + str(zone_id)  # Class for the zone
    final_weight_factor = final_weight_dbf[final_weight_dbf.Factor == "Landcover"][
        zone_class
    ].iloc[0]
    landcover_tif = landcover_tif * final_weight_factor

    # Write in memory
    temp_dir = os.path.join(
        tmp_dir, AnyPath("landcover_" + str(zone_id) + "_" + str(counter) + "_utm.tif")
    )
    output_dir = os.path.join(
        tmp_dir, AnyPath("landcover_" + str(zone_id) + "_" + str(counter) + ".tif")
    )

    rasters.write(landcover_tif, temp_dir, compress="deflate", predictor=1)
    rasters.write(
        landcover_tif.rio.reproject(
            proj_crs, resolution=output_resolution, resampling=Resampling.bilinear
        ),
        output_dir,
        predictor=1,
    )
    # From UTM to LatLon
    # with rio.open(temp_dir) as src:
    #        src_crs = src.crs
    #        src_transform = src.transform
    #        src_height = src.height
    #        src_width = src.width

    #        dst_transform, dst_width, dst_height = calculate_default_transform(
    #            src_crs, proj_crs, src_width, src_height, *src.bounds)
    #        # Define the metadata for the destination file
    #         dst_meta = src.meta.copy()
    #         dst_meta.update({
    #             'crs': proj_crs,
    #             'transform': dst_transform,
    #             'width': dst_width,
    #             'height': dst_height
    #         })

    #         # Create a new file to save the reprojected raster
    #         with rio.open(output_dir, 'w', **dst_meta) as dst:
    #             for i in range(1, src.count + 1):
    #                 reproject(
    #                     source=rio.band(src, i),
    #                     destination=rio.band(dst, i),
    #                     src_transform=src_transform,
    #                     src_crs=src_crs,
    #                     dst_transform=dst_transform,
    #                     dst_crs=proj_crs,
    #                     resampling=Resampling.nearest)

    return output_dir


def slope_raster_eu(
    slope_weight_path,
    dem,
    aoi_zone,
    proj_crs,
    zone_id,
    counter,
    final_weight_dbf,
    tmp_dir,
    output_resolution,
):

    LOGGER.info(" --- LSI : Slope")
    # -- Slope raster computation

    slope_dbf = gpd.read_file(slope_weight_path)
    slope_degrees = rasters.slope(dem, in_rad=False)
    with warnings.catch_warnings():
        warnings.simplefilter(
            "ignore"
        )  # the buffer is being applied in degrees, not meters
        buffer_slope = 200
        aoi_b = geometry.buffer(aoi_zone, buffer_slope)
    # slope_degrees = rasters.crop(slope_degrees, aoi_b)

    # -- Slope Reclassification:
    # Define steps of classification depending on the zone
    if int(zone_id) == 0:
        SLOPE_STEPS = [0, 1, 9, 13, 21, 27, 35, 42, 90]
    elif int(zone_id) in range(1, 4):
        SLOPE_STEPS = [0, 1, 5, 9, 13, 17, 21, 31, 90]
    elif int(zone_id) in range(5, 6):
        SLOPE_STEPS = [0, 1, 4, 8, 13, 18, 26, 38, 90]
    SLOPE_CLASSES = {
        1: f"{SLOPE_STEPS[0]} - {SLOPE_STEPS[1]}",
        2: f"{SLOPE_STEPS[1]} - {SLOPE_STEPS[2]}",
        3: f"{SLOPE_STEPS[2]} - {SLOPE_STEPS[3]}",
        4: f"{SLOPE_STEPS[3]} - {SLOPE_STEPS[4]}",
        5: f"{SLOPE_STEPS[4]} - {SLOPE_STEPS[5]}",
        6: f"{SLOPE_STEPS[4]} - {SLOPE_STEPS[5]}",
        7: f"{SLOPE_STEPS[4]} - {SLOPE_STEPS[5]}",
        8: f"{SLOPE_STEPS[4]} - {SLOPE_STEPS[5]}",
        9: f"{SLOPE_STEPS[5]}",
    }
    # Apply classification
    slope_name = "Value"
    slope_arr = classify_raster(slope_degrees, SLOPE_STEPS, SLOPE_CLASSES)
    slope_d = slope_degrees.copy(data=slope_arr).astype(np.float32).rename(slope_name)
    slope_d.attrs["long_name"] = slope_name
    slope_d = rasters.crop(slope_d, aoi_b)
    slope_gdf = xr_to_gdf(slope_d, proj_crs)
    # -- JOIN Slope with Weights
    slope_tif = slope_gdf.merge(slope_dbf, on="Value")
    slope_tif = slope_tif.set_index(["y", "x"]).Weight.to_xarray()
    slope_tif = slope_tif.rio.write_crs(slope_d.rio.crs)
    slope_tif = rasters.crop(slope_tif, aoi_zone)

    # -- Apply Final Weights
    zone_class = "Z" + str(zone_id)  # Class for the zone
    final_weight_factor = final_weight_dbf[final_weight_dbf.Factor == "Slope"][
        zone_class
    ].iloc[0]
    # Apply factor
    slope_tif = slope_tif * final_weight_factor

    # Write in memory
    temp_dir = os.path.join(
        tmp_dir, AnyPath("slope_" + str(zone_id) + "_" + str(counter) + "utm.tif")
    )
    output_dir = os.path.join(
        tmp_dir, AnyPath("slope_" + str(zone_id) + "_" + str(counter) + ".tif")
    )

    rasters.write(slope_tif, temp_dir, compress="deflate", predictor=1)
    rasters.write(
        slope_tif.rio.reproject(
            proj_crs, resolution=output_resolution, resampling=Resampling.bilinear
        ),
        output_dir,
        predictor=1,
    )

    # # From UTM to LatLon
    # with rio.open(temp_dir) as src:
    #         src_crs = src.crs
    #         src_transform = src.transform
    #         src_height = src.height
    #         src_width = src.width

    #         dst_transform, dst_width, dst_height = calculate_default_transform(
    #             src_crs, proj_crs, src_width, src_height, *src.bounds)
    #         # Define the metadata for the destination file
    #         dst_meta = src.meta.copy()
    #         dst_meta.update({
    #             'crs': proj_crs,²
    #             'transform': dst_transform,
    #             'width': dst_width,
    #             'height': dst_height
    #         })

    #         # Create a new file to save the reprojected raster
    #         with rio.open(output_dir, 'w', **dst_meta) as dst:
    #             for i in range(1, src.count + 1):
    #                 reproject(
    #                     source=rio.band(src, i),
    #                     destination=rio.band(dst, i),
    #                     src_transform=src_transform,
    #                     src_crs=src_crs,
    #                     dst_transform=dst_transform,
    #                     dst_crs=proj_crs,
    #                     resampling=Resampling.nearest)

    return output_dir


def mosaicing(raster_list, proj_crs, output_path, name):
    """
    This function allows to mosaic the rasters by zone
    Args:
        raster_list: list of xarrays
        proj_crs: CRS
        output_path: Path to be written the new raster
        name: str, name for the raster
    Returns:
        Path for the Mosaic raster in xarray format
    """
    src_files_to_mosaic = []
    for raster_path in raster_list:
        src = rio.open(raster_path)
        src_files_to_mosaic.append(src)
    # Merge the rasters
    mosaic, out_trans = merge(src_files_to_mosaic)

    # Copy the metadata
    out_meta = src_files_to_mosaic[0].meta.copy()
    out_meta.update(
        {
            "driver": "GTiff",
            "height": mosaic.shape[1],
            "width": mosaic.shape[2],
            "transform": out_trans,
            "crs": proj_crs,
        }
    )

    output_path = os.path.join(output_path, AnyPath(str(name) + "_mosaic.tif"))

    with rio.open(output_path, "w", **out_meta) as dest:
        dest.write(mosaic)

    # Close all source files
    for src in src_files_to_mosaic:
        src.close()

    return output_path

# --- LSI computation
def lsi_core(input_dict: dict) -> None:
    """
    Make a dict with the rasters:
        1. Geology
        2. Slope (in degrees)
        3. Landuse
        4. Elevation
        5. Distance to hydro_L
        6. Aspect
    from the input dict
    Args:
        input_dict (dict) : Dict that store parameters values
    Returns:
        dict : Dict that store raster to be processed for LSI
        copmputation
    """
    # --- Check parameters ---

    check_parameters(input_dict)

    # --- Extract parameters ---

    aoi_path = input_dict.get(InputParameters.AOI_PATH.value)
    location = input_dict.get(InputParameters.LOCATION.value)
    dem_name = input_dict.get(InputParameters.DEM_NAME.value)
    other_dem_path = input_dict.get(InputParameters.OTHER_DEM_PATH.value)
    #    lithology_path = input_dict.get(InputParameters.LITHOLOGY_PATH.value) #not anymore
    landcover_name = input_dict.get(InputParameters.LANDCOVER_NAME.value)
    europe_method = input_dict.get(InputParameters.EUROPE_METHOD.value)
    #    weights_path = input_dict.get(InputParameters.WEIGHTS_PATH.value) #not anymore
    output_resolution = input_dict.get(InputParameters.OUTPUT_RESOLUTION.value)
    epsg_code = input_dict.get(InputParameters.REF_EPSG.value)
    output_path = input_dict.get(InputParameters.OUTPUT_DIR.value)

    # Define folder for temporal files
    tmp_dir = os.path.join(output_path, "temp_dir")
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # -- Read AOI
    aoi = vectors.read(aoi_path)

    # -- Define EPSG
    if epsg_code:
        proj_crs = str("EPSG:" + str(epsg_code))
    else:
        proj_crs = aoi.estimate_utm_crs()

    # Reproject aoi to UTM crs
    aoi = aoi.to_crs(proj_crs)
    # maybe you'll have to re-write a utm crs projection

    # --- Define standard paths for both locations EUROPE and GLOBAL
    lithology_path = str(DataPath.LITHOLOGY_PATH)

    #  Dict that store dem name and dem path
    dem_path_dict = {
        DemType.COPDEM_30.value: DataPath.COPDEM30_PATH,
        DemType.FABDEM.value: DataPath.FABDEM_PATH,
        DemType.SRTM.value: DataPath.SRTM30_PATH,
        DemType.OTHER.value: other_dem_path,
    }
    # Store DEM path in a variable
    dem_path = dem_path_dict[dem_name]

    LOGGER.info(dem_path)

    # Reading and Checking errors in DEM
    try:
        buffer = 300  # A buffer that will work for all processes. (Can be recropped depending on the need)
        aoi_b = geometry.buffer(aoi, buffer)
        # dem_b = rasters.crop(dem_path, aoi_b)
        dem_b = rasters.read(dem_path, window=aoi_b)
    except ValueError:
        raise ValueError(
            "Your AOI doesn't cover your DTM. Make sure your input data is valid."
        )
    try:
        if np.unique(dem_b) == 0 or np.unique(dem_b) == np.nan:
            raise ValueError(
                "Your DEM has no elevations. Make sure your input DEM is valid"
            )
    except:
        if np.unique(dem_b).all == 0 or np.unique(dem_b).all == np.nan:
            raise ValueError(
                "Your DEM has no elevations. Make sure your input DEM is valid"
            )
    LOGGER.info(dem_path)
    # -- Define Resolution
    if output_resolution:
        output_resolution = int(output_resolution)
    else:
        x_res, y_res = dem_b.rio.resolution()
        output_resolution = int(np.round(abs((x_res) + abs(y_res) / 2)))

    dem_b = rasters.crop(dem_path, aoi)

    #LOGGER.info(dem_b.min(), dem_b.max())

    rasters.write(dem_b, AnyPath(output_path) / "dem_original.tif")
    # -- Reprojecting DEMS
    # DEM will be used as input for SLOPE, ELEVATION, ASPECT and HYDRO layers. Also as reference for Rasterization of Geology Layer.
    dem_b = dem_b.rio.reproject(
        proj_crs, resolution=output_resolution, resampling=Resampling.bilinear
    )
    dem_b = rasters.crop(dem_b, aoi_b)
    dem_b = dem_b.rio.write_nodata(FLOAT_NODATA)
    rasters.write(dem_b, AnyPath(output_path) / "dem_b.tif")

    # Define inputs for ELSUS (inside of Europe) method
    physio_zones_path = str(DataPath.ELSUS_ZONES_PATH)

    # -- Define Weights dbfs paths:
    if location == LocationType.EUROPE.value:
        elsus_weights_path = (
            DataPath.WEIGHTS_EUROPE_PATH
        )  # The path for each Zone weights will later
        # created based on this path
        elsus_final_weights_dbf_path = str(
            DataPath.WEIGHTS_EUROPE_PATH / "FinalWeigths.dbf"
        )

        geology_dbf_path = str(DataPath.WEIGHTS_GLOBAL_PATH / "Geology.dbf")
        global_final_weights_dbf_path = str(
            DataPath.WEIGHTS_GLOBAL_PATH / "Final_weights.dbf"
        )

    elif location == LocationType.GLOBAL.value:
        geology_dbf_path = str(DataPath.WEIGHTS_GLOBAL_PATH / "Geology.dbf")
        slope_dbf_path = str(DataPath.WEIGHTS_GLOBAL_PATH / "Slope.dbf")
        elevation_dbf_path = str(DataPath.WEIGHTS_GLOBAL_PATH / "Elevation.dbf")
        aspect_dbf_path = str(DataPath.WEIGHTS_GLOBAL_PATH / "Aspect.dbf")
        landuse_dbf_path = str(DataPath.WEIGHTS_GLOBAL_PATH / "Land use.dbf")
        hydro_dbf_path = str(DataPath.WEIGHTS_GLOBAL_PATH / "Hydro.dbf")
        global_final_weights_dbf_path = str(
            DataPath.WEIGHTS_GLOBAL_PATH / "Final_weights.dbf"
        )

    # Store landcover path in a variable depending on LOCATION
    # -- Check location (EUROPE or OUTSIDE of Europe)
    if location == LocationType.EUROPE.value:
        LOGGER.info("-- LSI - LOCATION : EUROPE")
        LOGGER.info(
            "-- LSI - CRS: "
            + str(proj_crs)
            + " | resolution: "
            + str(output_resolution)
        )

        if europe_method == EUMethodType.FAST.value:
            LOGGER.info("-- LSI -- Method: Europe fast method based on ELSUS layer")

            # IF EUROPE, we take the ELSUS layer at the AOI
            # -- Define path for LULC
            elsus_path = DataPath.ELSUS_PATH

            # -- Crop classic ELSUS layer by AOI
            lsi_tif = rasters.crop(elsus_path, aoi)
            lsi_tif = lsi_tif.rio.write_crs(proj_crs)
            lsi_tif = lsi_tif.rio.reproject(
                proj_crs, resolution=output_resolution, resampling=Resampling.bilinear
            )
            lsi_tif.rasters.crop(lsi_tif, aoi)
            lsi_tif.name = "Value"

            LOGGER.info("-- Writing lsi.tif in memory")

            # Write in memory
            rasters.write(lsi_tif, os.path.join(output_path, "lsi.tif"))
            return
        elif landcover_name == LandcoverType.ESAWC.value:
            LOGGER.info(
                "-- LSI -- Method: Europe Refined layer based on ESA WorldCover"
            )
            lulc_path = DataPath.ESAWC_PATH

        elif landcover_name == LandcoverType.CLC.value:
            LOGGER.info(
                "-- LSI -- Method: Europe Refined layer based on Corine Land Cover"
            )
            # -- Define path for LULC
            lulc_path = DataPath.CLC_PATH

        # -- 0. DEM
        LOGGER.info("-- Reading DEM")
        # with warnings.catch_warnings():
        #    warnings.simplefilter("ignore") # the buffer is being applied in degrees, not meters
        #    aoi_b = geometry.buffer(aoi, buffer_n) # buffer (for LandUse + Hydro rasters)
        # dem_b = rasters.crop(dem_path, aoi_b)

        LOGGER.info("-- Crop DEM")
        dem = rasters.crop(dem_b, aoi)

        # -- 1. Geology
        # Reading geology database, clip to aoi and reproject to proj_crs
        litho_db = gpd.read_file(lithology_path, driver="FileGDB", mask=aoi)
        LOGGER.info("-- AOI to litho_db crs")
        aoi_db = aoi.to_crs(litho_db.crs)
        litho_shp = gpd.clip(litho_db, aoi_db)
        LOGGER.info("-- Litho shp to proj_crs")
        litho_shp = litho_shp.to_crs(proj_crs)

        geology_dbf = gpd.read_file(geology_dbf_path)
        # momentaneous line to add Not Applicable class
        geology_dbf.loc[len(geology_dbf)] = ["Not Applicable", 997, 0.0, 0.0, None]
        geology_tif = geology_raster(
            geology_dbf, litho_shp, dem, aoi, proj_crs, output_path
        )

        # -- 2. Landcover + Slope (calculations per zone in the AOI according to ELSUS)

        LOGGER.info("-- Defining physio zones for Europe Refined method")

        # -- Read Climate-physiographic zones
        physio_zones = gpd.read_file(physio_zones_path)
        physio_zones = physio_zones.to_crs(proj_crs)
        physio_zones_aoi = gpd.clip(physio_zones, aoi)
        physio_zones_aoi = physio_zones_aoi.explode(index_parts=False)

        # Define a mapping of Zones to the Zone_Weights database file
        zone_to_dbf = {
            0: "Zone0",
            1: "Zone1",
            2: "Zone2",
            3: "Zone3",
            4: "Zone4",
            5: "Zone5",
            6: "Zone6",
        }

        LOGGER.info(
            str(
                "-- Produce the Land use and Slope raster for the LSI model | Amount of sub-zones to be processed: "
                + str(physio_zones_aoi["Zone"].count())
            )
        )

        i = 0
        landcover_list = []
        slope_list = []
        for index, row in physio_zones_aoi.iterrows():
            i += 1
            # The Climatic Zone
            zone = row["Zone"]
            db_file = zone_to_dbf[zone]
            # Extracting the shapefile for the Climatic Zone
            zone_geom = gpd.GeoDataFrame(row).T.set_geometry("geometry")
            zone_geom = zone_geom.set_crs(proj_crs)

            # -- Final weights
            fw_dbf = gpd.read_file(str(elsus_final_weights_dbf_path))

            # ---- LANDCOVER CASE ----
            if zone != 0:  # For Zone0 the Landcover is not used
                # -- Define weights path
                # landcover_w_path = os.path.join(elsus_weights_path,  str(db_file) + "Landcover_Z" + str(zone) + ".dbf")
                landcover_w_path = str(
                    elsus_weights_path
                    / str(str(db_file) + "/Landcover_Z" + str(zone) + ".dbf")
                )
            else:
                # landcover_w_path = os.path.join(elsus_weights_path, "Zone1/Landcover_Z1.dbf") # A random file is chosen. At the end
                # the weights are set to 0 for Zone0
                landcover_w_path = str(
                    elsus_weights_path / "Zone1/Landcover_Z1.dbf"
                )  # A random file is chosen. At the end
                # the weights are set to 0 for Zone0
            # -- Compute the Rasters
            landcover_dir = landcover_raster_eu(
                landcover_w_path,
                lulc_path,
                dem,
                zone_geom,
                proj_crs,
                zone,
                i,
                fw_dbf,
                tmp_dir,
                output_resolution,
                landcover_name,
            )
            landcover_list.append(landcover_dir)

            # ---- SLOPE and LITHOLOGY ----
            # -- Define weights path
            # slope_w_path = os.path.join(elsus_weights_path,  str(db_file) + "/Slope_Z" + str(zone) + ".dbf")
            slope_w_path = str(
                elsus_weights_path / str(str(db_file) + "/Slope_Z" + str(zone) + ".dbf")
            )
            # -- Compute the Rasters
            slope_dir = slope_raster_eu(
                slope_w_path,
                dem_b,
                zone_geom,
                proj_crs,
                zone,
                i,
                fw_dbf,
                tmp_dir,
                output_resolution,
            )

            slope_list.append(slope_dir)

        # Mosaicing the rasters by zone into the complete Raster
        slope_mosaic_path = mosaicing(
            slope_list, proj_crs, output_path, "slope_weight.tif"
        )
        landcover_mosaic_path = mosaicing(
            landcover_list, proj_crs, output_path, "landcover_weight.tif"
        )

        slope_tif = rasters.read(slope_mosaic_path)
        landcover_tif = rasters.read(landcover_mosaic_path)

        # Interpolate landcover and geology layers
        landcover_tif = landcover_tif.interp_like(slope_tif, method="nearest")
        geology_tif = geology_tif.interp_like(slope_tif, method="nearest")
        geology_tif = rasters.collocate(
            slope_tif, geology_tif, Resampling.bilinear
        )  # Collocate Geology for LSI computation

        # LSI computation:
        LOGGER.info("-- Computing LSI")
        fw_dbf = gpd.read_file(global_final_weights_dbf_path)
        geology_weights = fw_dbf[fw_dbf.Factors == "Geology"].Weights.iloc[0]

        lsi_tif = slope_tif + geology_tif * float(geology_weights) + landcover_tif

    elif location == LocationType.GLOBAL.value:
        LOGGER.info("-- LSI - LOCATION : GLOBAL")
        # Print general information about CRS and Resolution
        LOGGER.info(
            "-- LSI - CRS: "
            + str(proj_crs)
            + " | resolution: "
            + str(output_resolution)
        )
        # -- Dict that store landcover name and landcover path
        landcover_path_dict = {
            LandcoverType.ESAWC.value: DataPath.ESAWC_PATH,
            LandcoverType.CLC.value: DataPath.CLC_PATH,
        }

        lulc_path = landcover_path_dict[landcover_name]

        # 0. DEM
        LOGGER.info("-- Reading DEM")
        # with warnings.catch_warnings():
        #    warnings.simplefilter("ignore") # the buffer is being applied in degrees, not meters
        #    aoi_b = geometry.buffer(aoi, buffer_n) # buffer (for LandUse + Hydro rasters)
        # try:
        #    dem_buff = rasters.crop(dem_path, aoi_b)
        # except ValueError:
        #    raise ValueError(
        #       "Your AOI doesn't cover your DEM. Make sure your input data is valid."
        #    )

        # dem_b = rasters.crop(dem_b, aoi_b) # second crop

        # with warnings.catch_warnings():
        #    warnings.simplefilter("ignore") # the buffer is being applied in degrees, not meters
        dem = rasters.crop(dem_b, aoi)
        dem = dem.rio.reproject(
            proj_crs, resolution=output_resolution, resampling=Resampling.bilinear
        )

        # with warnings.catch_warnings():
        #     warnings.simplefilter("ignore") # the buffer is being applied in degrees, not meters
        #     dem = rasters.crop(dem, aoi)

        # -- 1. Geology
        # Reading geology database, clip to aoi and reproject to proj_crs
        litho_db = gpd.read_file(lithology_path, driver="FileGDB", mask=aoi)
        aoi_db = aoi.to_crs(litho_db.crs)
        litho_shp = gpd.clip(litho_db, aoi_db)
        litho_shp = litho_shp.to_crs(proj_crs)

        geology_dbf = gpd.read_file(geology_dbf_path)
        # momentaneous line to add Not Applicable class
        geology_dbf.loc[len(geology_dbf)] = ["Not Applicable", 997, 0.0, 0.0, None]
        geology_layer = geology_raster(
            geology_dbf, litho_shp, dem, aoi, proj_crs, output_path
        )

        # -- 2. Slope
        slope_dbf = gpd.read_file(slope_dbf_path)
        slope_layer = slope_raster(slope_dbf, dem_b, aoi, proj_crs, output_path)

        # -- 3. Landcover
        landuse_dbf = gpd.read_file(landuse_dbf_path)
        landuse_dbf.loc[len(landuse_dbf)] = ["Not Applicable", 997, 0.0, 0.0, None]
        # Landcover
        landuse_layer = landcover_raster(
            landuse_dbf,
            lulc_path,
            landcover_name,
            aoi,
            proj_crs,
            output_resolution,
            output_path,
        )

        # -- 4. Elevation
        elevation_dbf = gpd.read_file(elevation_dbf_path)
        elevation_layer = elevation_raster(
            elevation_dbf, dem_b, aoi, proj_crs, output_path
        )

        # -- 5. Hydro
        hydro_dbf = gpd.read_file(hydro_dbf_path)
        hydro_layer = hydro_raster(
            hydro_dbf,
            dem_b,
            aoi,
            proj_crs,
            elevation_layer,
            output_resolution,
            tmp_dir,
            output_path,
        )

        # -- 6. Aspect
        aspect_dbf = gpd.read_file(aspect_dbf_path)
        aspect_layer = aspect_raster(aspect_dbf, dem_b, aoi, proj_crs, output_path)

        # -- Final weights computing
        fw_dbf = gpd.read_file(global_final_weights_dbf_path)

        # Extracting final weights
        slope_weights = fw_dbf[fw_dbf.Factors == "Slope"].Weights.iloc[0]
        geology_weights = fw_dbf[fw_dbf.Factors == "Geology"].Weights.iloc[0]
        aspect_weights = fw_dbf[fw_dbf.Factors == "Slope aspect"].Weights.iloc[0]
        elevation_weights = fw_dbf[fw_dbf.Factors == "Elevation"].Weights.iloc[0]
        hydro_weights = fw_dbf[fw_dbf.Factors == "Distance from river"].Weights.iloc[0]
        landuse_weights = fw_dbf[fw_dbf.Factors == "Land use"].Weights.iloc[0]

        # Final weight
        LOGGER.info("-- Computing LSI")

        # Read layers
        slope_layer = rasters.read(os.path.join(output_path, "slope_weight.tif"))
        geology_layer = rasters.read(os.path.join(output_path, "geology_weight.tif"))
        elevation_layer = rasters.read(
            os.path.join(output_path, "elevation_weight.tif")
        )
        aspect_layer = rasters.read(os.path.join(output_path, "aspect_weight.tif"))
        landuse_layer = rasters.read(os.path.join(output_path, "landcover_weight.tif"))
        hydro_layer = rasters.read(os.path.join(output_path, "hydro_weight.tif"))

        # Collocate layers

        geology_layer = rasters.collocate(
            slope_layer, geology_layer, Resampling.bilinear
        )
        aspect_layer = rasters.collocate(slope_layer, aspect_layer, Resampling.bilinear)
        aspect_layer = rasters.crop(aspect_layer, aoi)
        landuse_layer = rasters.collocate(
            slope_layer, landuse_layer, Resampling.bilinear
        )
        hydro_layer = rasters.collocate(slope_layer, hydro_layer, Resampling.bilinear)

        lsi_tif = (
            slope_layer * float(slope_weights)
            + geology_layer * float(geology_weights)
            + elevation_layer * float(elevation_weights)
            + aspect_layer * float(aspect_weights)
            + landuse_layer * float(landuse_weights)
            + hydro_layer * float(hydro_weights)
        )

    # Clean errors in LSI due to reprojection [?]

    # Fix coordinates
    lsi_tif = lsi_tif.rio.write_crs(proj_crs, inplace=True)
    # lsi_tif = rasters.collocate(  rasters.read(os.path.join(output_path, "aspect_weight.tif"))
    #                             , lsi_tif
    #                             , Resampling.bilinear)
    # lsi_tif = lsi_tif.rio.write_nodata(np.nan, encoded=True)
    lsi_tif = rasters.crop(lsi_tif, aoi)
    lsi_tif = xr.where(
        lsi_tif > 10, np.nan, lsi_tif
    )  # There should not be values greater than 10 (border effect)
    lsi_tif = xr.where(
        lsi_tif < 0, np.nan, lsi_tif
    )  # There should not be negative values (border effect)
    lsi_tif = lsi_tif.rio.write_crs(proj_crs, inplace=True)

    LOGGER.info("-- Writing lsi.tif in memory")

    # Write in memory
    rasters.write(lsi_tif, os.path.join(output_path, "lsi.tif"))

    return
    # raise NotImplementedError
