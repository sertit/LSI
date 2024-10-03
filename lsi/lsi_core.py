# -*- coding: utf-8 -*-
# This file is part of LSI.
# LSI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# LSI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with LSI. If not, see <https://www.gnu.org/licenses/>.
""" lsi_core """

import logging
import os
import shutil
import warnings
from enum import unique

import geopandas as gpd
import numpy as np

# import rasterio as rio
import xarray as xr
from rasterio.enums import Resampling
from sertit import AnyPath, rasters, vectors
from sertit.misc import ListEnum
from sertit.rasters import FLOAT_NODATA
from sertit.unistra import get_geodatastore

from lsi.src.lsi_calculator import (  # lithology_raster_eu,; hydro_raster,
    aspect_raster,
    elevation_raster,
    geology_raster,
    hydro_raster_wbw,
    landcover_raster,
    landcover_raster_eu,
    slope_raster,
    slope_raster_eu,
)
from lsi.src.reclass import LandcoverType
from lsi.src.utils import (
    compute_statistics,
    mosaicing,
    produce_a_reclass_arr,
    raster_postprocess,
)

DEBUG = False
LOGGING_FORMAT = "%(asctime)s - [%(levelname)s] - %(message)s"
LOGGER = logging.getLogger("LSI")

REGULAR_BUFFER = 1000
GADM_BUFFER = 15000  # Big buffer to avoid missing departments


def geodatastore(ftep=False):
    """
    This function returns the root path to the geo data store (DEM, Soil database...).
    Args:
        ftep: If True, the path to the s3 bucket for the ftep platform is returned. Else, get_geodatastore from sertit utils module is called.

    Returns:
        The path to the geo data store.

    """
    if ftep:
        return AnyPath("s3://eo4sdg-data")
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
            / "CLC2018_CLC2018_V2018_20_epsg4326.tif"  # projected to EPSG:4326
        )
        cls.GLC_PATH = (
            cls.GLOBAL_DIR
            / "Global_Land_Cover"
            / "2019"
            / "PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif"
        )
        cls.ELC_PATH = (
            cls.GLOBAL_DIR
            / "2021"
            / "lulc2021"
            / "lc2021"
            / "38T_20230101-20240101.tif"
            / ".vrt"  # Where to find it??
        )
        #        cls.EUDEM_PATH = cls.GLOBAL_DIR / "EUDEM_v2" / "eudem_dem_3035_europe.tif"
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
        cls.LITHOLOGY_PATH = cls.GLOBAL_DIR / "LiMW_GIS_2015" / "LiMW_GIS_2015.gdb"
        cls.ELSUS_ZONES_PATH = (
            cls.GLOBAL_DIR
            / "ELSUSv2"
            / "ELSUSv2_six_datasets"
            / "climate_phys_regions.shp"
        )
        cls.LITHOLOGY_ELSUS_PATH = (
            cls.GLOBAL_DIR
            / "ELSUSv2"
            / "ELSUSv2_six_datasets"
            / "lithology_epsg4326.tif"  # projected to epsg4326
        )
        cls.GADM_PATH = cls.GLOBAL_DIR / "GADM" / "gadm_410.gdb"


@unique
class InputParameters(ListEnum):
    """
    List of the input parameters
    """

    AOI_PATH = "aoi_path"
    LOCATION = "location"
    DEM_NAME = "dem_name"
    OTHER_DEM_PATH = "other_dem_path"
    LANDCOVER_NAME = "landcover_name"
    EUROPE_METHOD = "europe_method"
    OUTPUT_RESOLUTION = "output_resolution"
    REF_EPSG = "ref_epsg"
    OUTPUT_DIR = "output_path"
    TMP_DIR = "temp_dir"
    TEMP = "temp"
    JENKS = "jenks_break"


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

    COPDEM_30 = "COPDEM 30m"
    FABDEM = "FABDEM"
    OTHER = "Other"


class LandcoverStructure(ListEnum):
    """
    List of the Landcover type (Coding)
    """

    CLC = "Corine Land Cover - 2018 (100m)"

    ESAWC = "ESA WorldCover - 2021 (10m)"

    GLC = "Global Land Cover - Copernicus 2019 (100m)"

    ELC = "ESRI Annual Land Cover 2021 (10m)"


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
    landcover_name = input_dict.get(InputParameters.LANDCOVER_NAME.value)

    # Check if other_dem_path is needed
    if (dem_name == DemType.OTHER.value) and (other_dem_path is None):
        raise ValueError(f'{"Dem path is needed !"}')
    if (location == LocationType.GLOBAL.value) and (
        landcover_name is LandcoverType.CLC.value
    ):
        raise ValueError(
            f'{"Corine Land Cover can not be used for GLOBAL calculations, only in Europe !"}'
        )
    return


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
    landcover_name = input_dict.get(InputParameters.LANDCOVER_NAME.value)
    europe_method = input_dict.get(InputParameters.EUROPE_METHOD.value)
    output_resolution = input_dict.get(InputParameters.OUTPUT_RESOLUTION.value)
    epsg_code = input_dict.get(InputParameters.REF_EPSG.value)
    output_path = input_dict.get(InputParameters.OUTPUT_DIR.value)
    temp = input_dict.get(InputParameters.TEMP.value)
    jenks_break = input_dict.get(InputParameters.JENKS.value)

    # Read GADM path for statistics
    gadm_path = str(DataPath.GADM_PATH)

    # Define folder for temporal files
    tmp_dir = os.path.join(output_path, "temp_dir")
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # -- Read AOI

    # Write wkt string input to shapefile
    if aoi_path.startswith("POLYGON"):
        aoi_gpd_wkt = gpd.GeoSeries.from_wkt([aoi_path])
        aoi_raw_path_wkt = os.path.join(tmp_dir, "aoi_from_wkt.shp")

        aoi_gpd_wkt_4326 = aoi_gpd_wkt.set_crs(epsg=4326)
        aoi_gpd_wkt_4326.to_file(aoi_raw_path_wkt)
        aoi_raw_path = aoi_raw_path_wkt

        aoi = gpd.read_file(aoi_raw_path)
    else:
        aoi = vectors.read(aoi_path)

    # -- Define EPSG
    if epsg_code:
        proj_crs = str("EPSG:" + str(epsg_code))
    else:
        proj_crs = aoi.estimate_utm_crs()

    # Reproject aoi to UTM crs
    try:  # In some cases the AOI has a corrupted crs
        aoi = aoi.to_crs(proj_crs)
    except:  # noqa
        aoi = aoi.set_crs(proj_crs, allow_override=True)
        # aoi = aoi.to_crs(proj_crs)

    #  Dict that store dem name and dem path
    dem_path_dict = {
        DemType.COPDEM_30.value: DataPath.COPDEM30_PATH,
        DemType.FABDEM.value: DataPath.FABDEM_PATH,
        DemType.OTHER.value: other_dem_path,
    }
    # Store DEM path in a variable
    dem_path = dem_path_dict[dem_name]

    # Reading and Checking errors in DEM
    try:
        aoi_b = aoi
        aoi_b = aoi_b.geometry.buffer(REGULAR_BUFFER)
        # dem_b = rasters.read(dem_path, window=aoi_b)
        dem_b = rasters.crop(dem_path, aoi_b)
    except ValueError:
        raise ValueError(
            "Your AOI doesn't cover your DTM. Make sure your input data is valid."
        )
    try:
        if np.unique(dem_b) == 0 or np.unique(dem_b) == np.nan:
            raise ValueError(
                "Your DEM has no elevations. Make sure your input DEM is valid"
            )
    except ValueError:
        if np.unique(dem_b).all == 0 or np.unique(dem_b).all == np.nan:
            raise ValueError(
                "Your DEM has no elevations. Make sure your input DEM is valid"
            )
        dem_max = dem_b.max()
        dem_min = dem_b.min()

    # -- Define Resolution
    if output_resolution:
        output_resolution = int(output_resolution)
    else:
        x_res, y_res = dem_b.rio.resolution()
        output_resolution = int(np.round(abs((x_res) + abs(y_res) / 2)))

    dem_b = rasters.crop(dem_path, aoi_b)

    # -- Reprojecting DEM
    # DEM will be used as input for SLOPE, ELEVATION, ASPECT and HYDRO layers. Also as reference for Rasterization of Geology Layer.
    dem_b = dem_b.rio.reproject(
        proj_crs, resolution=output_resolution, resampling=Resampling.bilinear
    )
    dem_b = dem_b.rio.write_nodata(FLOAT_NODATA)
    rasters.write(dem_b, AnyPath(tmp_dir) / "dem_b.tif")

    # -- DEFINING WEIGHT DBF PATHs and LITHOLOGY PATH:
    if location == LocationType.EUROPE.value:
        # Define lithology
        lithology_path = str(DataPath.LITHOLOGY_PATH)
        # lithology_path = str(DataPath.LITHOLOGY_ELSUS_PATH)

        # Define inputs for ELSUS (inside of Europe) method
        physio_zones_path = str(DataPath.ELSUS_ZONES_PATH)

        elsus_weights_path = (
            DataPath.WEIGHTS_EUROPE_PATH
        )  # The path for each Zone weights will later
        # created based on this path
        elsus_final_weights_dbf_path = str(
            DataPath.WEIGHTS_EUROPE_PATH / "FinalWeigths.dbf"
        )

        # geology_dbf_path = str(DataPath.WEIGHTS_GLOBAL_PATH / "Geology.dbf")
        global_final_weights_dbf_path = str(
            DataPath.WEIGHTS_GLOBAL_PATH / "Final_weights.dbf"
        )

    elif location == LocationType.GLOBAL.value:
        # Define lithology
        lithology_path = str(DataPath.LITHOLOGY_PATH)

        # Define weights path
        # geology_dbf_path = str(DataPath.WEIGHTS_GLOBAL_PATH / "Geology.dbf")
        # slope_dbf_path = str(DataPath.WEIGHTS_GLOBAL_PATH / "Slope.dbf")
        # elevation_dbf_path = str(DataPath.WEIGHTS_GLOBAL_PATH / "Elevation.dbf")
        # aspect_dbf_path = str(DataPath.WEIGHTS_GLOBAL_PATH / "Aspect.dbf")
        # landuse_dbf_path = str(DataPath.WEIGHTS_GLOBAL_PATH / "Land use.dbf")
        # hydro_dbf_path = str(DataPath.WEIGHTS_GLOBAL_PATH / "Hydro.dbf")
        global_final_weights_dbf_path = str(
            DataPath.WEIGHTS_GLOBAL_PATH / "Final_weights.dbf"
        )

    # --  LSI COMPUTATION ACCORDING TO METHOD
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
        LOGGER.info("-- Crop DEM")
        dem = rasters.crop(dem_b, aoi)

        # -- 1. Geology (to be erased for new Lithology dataset if ELSUS lithology is used)

        # Reading geology database, clip to aoi and reproject to proj_crs
        litho_db = gpd.read_file(lithology_path, driver="FileGDB", mask=aoi)
        aoi_db = aoi.to_crs(litho_db.crs)
        litho_shp = gpd.clip(litho_db, aoi_db)
        litho_shp = litho_shp.to_crs(proj_crs)

        # geology_dbf = gpd.read_file(geology_dbf_path)
        # geology_dbf.loc[len(geology_dbf)] = ["Not Applicable", 997, 0.0, 0.0, None]
        geology_tif = geology_raster(litho_shp, dem, aoi, proj_crs, tmp_dir)

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
                "-- Produce the Landuse, Slope and Lithology raster for the LSI model | Amount of sub-zones to be processed: "
                + str(physio_zones_aoi["Zone"].count())
            )
        )

        # A counter i is set to separate polygons from the same Zone
        i = 0

        # The following loop crops the rasters by each physiographical zone in the AOI and stores the path of each
        # raster in the lists below
        landcover_list = []
        slope_list = []
        # lithology_list = []
        for _, row in physio_zones_aoi.iterrows():
            i += 1
            # The Climatic Zone
            zone = row["Zone"]
            db_file = zone_to_dbf[zone]

            # Extracting the shapefile for the Climatic Zone
            zone_geom = gpd.GeoDataFrame(row).T.set_geometry("geometry")
            zone_geom = zone_geom.set_crs(proj_crs)

            # -- Read Final weights
            fw_dbf = gpd.read_file(str(elsus_final_weights_dbf_path))

            # ---- LANDCOVER CASE ----
            if (
                zone != 0
            ):  # For Zone0 the Landcover is not used as it representes water bodies
                # -- Define weights path
                landcover_w_path = str(
                    elsus_weights_path
                    / str(str(db_file) + "/Landcover_Z" + str(zone) + ".dbf")
                )
            else:  # the weights are set to 0 for Zone0
                # A random file is chosen (in this case Zone1). which file i is, is not important because at the end
                # the weights are set to 0 for Zone0
                landcover_w_path = str(elsus_weights_path / "Zone1/Landcover_Z1.dbf")

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

            # ---- SLOPE CASE ----
            # -- Define weights path
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

            # ---- LITHOLOGY CASE ----
            # -- Define weights path
            # lithology_w_path = str(
            #     elsus_weights_path / str(str(db_file) + "/Lithology_Z" + str(zone) + ".dbf")
            # )
            # -- Compute the Rasters
            # lithology_dir = lithology_raster_eu(
            #     lithology_w_path,
            #     lithology_path,
            #     zone_geom,
            #     proj_crs,
            #     zone,
            #     i,
            #     fw_dbf,
            #     tmp_dir,
            #     output_resolution,
            # )
            # lithology_list.append(lithology_dir)

        # Mosaicing the rasters by zone available in the lists previously mentioned
        slope_mosaic_path, _ = mosaicing(slope_list, tmp_dir, "slope_weight")
        landcover_mosaic_path, _ = mosaicing(
            landcover_list, tmp_dir, "landcover_weight"
        )
        # lithology_mosaic_path, _ = mosaicing(
        #     lithology_list, tmp_dir, "lithology_weight"
        # )

        # Read the mosaiced rasters
        slope_tif = rasters.crop(slope_mosaic_path, aoi)
        landcover_tif = rasters.crop(landcover_mosaic_path, aoi)
        # lithology_tif = rasters.crop(lithology_mosaic_path, aoi)

        # Interpolate landcover and lithology layers to be equivalent to Slope layer
        landcover_tif = landcover_tif.interp_like(slope_tif, method="nearest")
        geology_tif = geology_tif.interp_like(slope_tif, method="nearest")
        # lithology_tif = lithology_tif.interp_like(slope_tif, method="nearest")

        # Collocate Geology for LSI computation
        geology_tif = rasters.collocate(slope_tif, geology_tif, Resampling.bilinear)

        # LSI computation:
        LOGGER.info("-- Computing LSI")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fw_dbf = gpd.read_file(global_final_weights_dbf_path)
            geology_weights = fw_dbf[fw_dbf.Factors == "Geology"].Weights.iloc[0]
            lsi_tif = slope_tif + geology_tif * float(geology_weights) + landcover_tif
            # lsi_tif = slope_tif + lithology_tif + landcover_tif

    elif location == LocationType.GLOBAL.value:
        LOGGER.info("-- LSI - LOCATION : GLOBAL")
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
            LandcoverType.GLC.value: DataPath.GLC_PATH,
            LandcoverType.ELC.value: DataPath.ELC_PATH,
        }
        lulc_path = landcover_path_dict[landcover_name]

        # 0. DEM
        LOGGER.info("-- Reading DEM")

        dem = rasters.crop(dem_path, aoi)
        dem = dem.rio.reproject(
            proj_crs, resolution=output_resolution, resampling=Resampling.bilinear
        )
        dem = dem.rio.write_nodata(FLOAT_NODATA)

        # -- 1. Geology
        # Reading geology database, clip to aoi and reproject to proj_crs
        litho_db = gpd.read_file(lithology_path, driver="FileGDB", mask=aoi)
        aoi_db = aoi.to_crs(litho_db.crs)
        litho_shp = gpd.clip(litho_db, aoi_db)
        litho_shp = litho_shp.to_crs(proj_crs)

        # Compute Geology layer
        geology_layer = geology_raster(litho_shp, dem, aoi, proj_crs, tmp_dir)

        # -- 2. Slope
        slope_layer = slope_raster(dem_b, aoi, proj_crs, tmp_dir)

        # -- 3. Landcover
        lulc = rasters.crop(lulc_path, aoi_b)
        lulc = rasters.collocate(dem_b, lulc, Resampling.nearest)

        # Compute landcover layer
        landuse_layer = landcover_raster(
            lulc,
            landcover_name,
            aoi,
            proj_crs,
            output_resolution,
            tmp_dir,
        )
        landuse_layer = rasters.crop(landuse_layer, aoi)
        landuse_layer = landuse_layer.rio.write_nodata(FLOAT_NODATA)

        # -- 4. Elevation
        elevation_layer = elevation_raster(dem_b, aoi, proj_crs, tmp_dir)

        # -- 5. Hydro
        hydro_layer = hydro_raster_wbw(
            dem_b,
            aoi,
            proj_crs,
            dem_max,
            dem_min,
            output_resolution,
            tmp_dir,
        )

        # -- 6. Aspect
        aspect_layer = aspect_raster(dem_b, aoi, proj_crs, tmp_dir)

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

        # Read all layers from temporal directory
        slope_layer = rasters.read(os.path.join(tmp_dir, "slope_weight.tif"))
        geology_layer = rasters.read(os.path.join(tmp_dir, "geology_weight.tif"))
        elevation_layer = rasters.read(os.path.join(tmp_dir, "elevation_weight.tif"))
        aspect_layer = rasters.read(os.path.join(tmp_dir, "aspect_weight.tif"))
        landuse_layer = rasters.read(os.path.join(tmp_dir, "landcover_weight.tif"))
        hydro_layer = rasters.read(os.path.join(tmp_dir, "hydro_weight.tif"))

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

        # Calculate lsi
        lsi_tif = (
            slope_layer * float(slope_weights)
            + geology_layer * float(geology_weights)
            + elevation_layer * float(elevation_weights)
            + aspect_layer * float(aspect_weights)
            + landuse_layer * float(landuse_weights)
            + hydro_layer * float(hydro_weights)
        )

    # -- LSI postprocessing: Fix coordinates
    lsi_tif = lsi_tif.rio.write_crs(proj_crs, inplace=True)
    lsi_tif = rasters.crop(lsi_tif, aoi)
    lsi_tif = xr.where(
        lsi_tif > 10, np.nan, lsi_tif
    )  # There should not be values greater than 10 (border effect)
    lsi_tif = xr.where(
        lsi_tif < 0, np.nan, lsi_tif
    )  # There should not be negative values (border effect)
    lsi_tif = lsi_tif.rio.write_crs(proj_crs, inplace=True)

    if (
        jenks_break
    ):  # Apply jenks only if user requires it (this option takes a longer time)
        # -- LSI reclassification: 1 to 5

        # Write in memory LSI with unclassified values
        rasters.write(lsi_tif, os.path.join(output_path, "lsi_unclassified.tif"))

        LOGGER.info("-- Reclassification of LSI in 5 vulnerability classes")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            lsi_reclass = produce_a_reclass_arr(lsi_tif, downsample_factor=450)

            # Vectorizing
            lsi_reclass_vector = rasters.vectorize(lsi_reclass)

            # Delete class 6 for values above the max value of LSI
            lsi_reclass_vector_c = lsi_reclass_vector[
                lsi_reclass_vector["raster_val"] != 6
            ]

            # Going back to raster
            lsi_reclass_tif = rasters.rasterize(
                lsi_tif, lsi_reclass_vector_c, value_field="raster_val"
            )

            lsi_reclass_tif = rasters.crop(lsi_reclass_tif, aoi)
            lsi_reclass_tif = xr.where(
                lsi_reclass_tif == 0, np.nan, lsi_reclass_tif
            )  # There should not be 0 values for classes) (border effect)
            lsi_tif = lsi_reclass_tif.rio.write_crs(proj_crs, inplace=True)

        LOGGER.info("-- Writing LSI in memory")
        # Post-processing
        # Currently there is an error with the sieving
        lsi_tif_sieved, lsi_vector = raster_postprocess(lsi_tif)
        vectors.write(lsi_vector, os.path.join(output_path, "lsi.shp"))

        # Write in memory
        rasters.write(lsi_tif, os.path.join(output_path, "lsi.tif"))
    else:
        LOGGER.info("-- Writing LSI in memory")
        # Write in memory LSI with unclassified values
        rasters.write(lsi_tif, os.path.join(output_path, "lsi_unclassified.tif"))

    LOGGER.info("-- Computing LSI statistics for unclassified LSI")

    # Read GADM layer and overlay with AOI
    aoi_gadm = aoi.geometry.buffer(GADM_BUFFER)
    gadm = vectors.read(gadm_path, window=aoi_gadm)
    gadm = gadm.to_crs(proj_crs)
    gadm = gpd.overlay(gadm, aoi)

    lsi_stats = compute_statistics(
        gadm, os.path.join(output_path, "lsi_unclassified.tif")
    )

    LOGGER.info("-- Writing LSI statistics in memory")
    # Write statistics in memory
    vectors.write(lsi_stats, os.path.join(output_path, "lsi_unclassified_stats.shp"))

    if not temp:
        LOGGER.info("-- Deleting temporary files")
        shutil.rmtree(tmp_dir)

    return
    # raise NotImplementedError
