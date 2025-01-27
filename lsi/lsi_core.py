"""
This file is part of LSI.
LSI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
LSI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with LSI. If not, see <https://www.gnu.org/licenses/>.
"""

"""lsi_core"""

import logging
import os
import shutil
import warnings
from enum import unique

import geopandas as gpd
import numpy as np
import xarray as xr
from rasterio.enums import Resampling
from sertit import AnyPath, rasters, vectors
from sertit.misc import ListEnum
from sertit.rasters import FLOAT_NODATA
from sertit.unistra import get_geodatastore

from lsi import __version__
from lsi.src.lsi_calculator import (
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
        os.environ["AWS_VIRTUAL_HOSTING"] = "False"
        return AnyPath("s3://eo4sdg-data-sertit")
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
        cls.GADM_PATH = cls.GLOBAL_DIR / "GADM" / "gadm_410.gdb"
        cls.COASTLINE_PATH = (
            cls.GLOBAL_DIR / "GAUL_2006_coast" / "gaul2006_coast_fixed.shp"
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


def check_parameters(input_dict: dir) -> None:
    """
     Check if parameters values are valid
    Args:
        input_dict (dict) : dict with parameters values

    Returns:

    """
    # --- Extract parameters ---
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

    LOGGER.info(f"*** Tool version: {__version__} ***")
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

    LOGGER.info("-- Reading inputs")
    # Read GADM path for statistics
    gadm_path = str(DataPath.GADM_PATH)
    # Read Coastline path for cropping AOIs that include water (GLOBAL case)
    coastline_path = str(DataPath.COASTLINE_PATH)

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

    # Removing the open water (or most of it at least) from the AOI by clip with coastline
    # (!) Only for GLOBAL method as the European method takes into consideration Coastlines
    # that could be erased for some AOIs where there is a coarser coastline delineation
    # from the layer used (GAUL 2006)
    if location == LocationType.GLOBAL.value:
        coastline = vectors.read(coastline_path, window=aoi)  # read with AOI window
        aoi = gpd.clip(aoi, coastline)

    # -- Define EPSG
    if epsg_code:
        proj_crs = str("EPSG:" + str(epsg_code))
    else:
        try:
            proj_crs = aoi.estimate_utm_crs()
        except ValueError:
            raise ValueError(
                "The AOI seems to be corrupted, please check it is not empty and try again."
            ) from ValueError

    # Reproject aoi to UTM crs
    try:  # In some cases the AOI has a corrupted crs
        aoi = aoi.to_crs(proj_crs)
    except:  # noqa
        aoi = aoi.set_crs(proj_crs, allow_override=True)

    # - Write original AOI to file
    aoi_o_path = os.path.join(tmp_dir, "aoi.shp")
    aoi.to_file(aoi_o_path)

    #  Dict that store dem name and dem path
    dem_path_dict = {
        DemType.COPDEM_30.value: DataPath.COPDEM30_PATH,
        DemType.FABDEM.value: DataPath.FABDEM_PATH,
        DemType.OTHER.value: other_dem_path,
    }
    # Store DEM path in a variable
    dem_path = dem_path_dict[dem_name]

    LOGGER.info("-- Reading DEM")
    # Reading and Checking errors in DEM + Buff AOI
    try:
        aoi_b = aoi.copy()
        aoi_b.geometry = aoi_b.geometry.buffer(REGULAR_BUFFER)

        # Export the new aoi
        aoi_b_path = os.path.join(tmp_dir, "aoi_buffREGULAR.shp")
        aoi_b.to_file(aoi_b_path)

        dem_b = rasters.read(dem_path, window=aoi_b)
        dem_b = rasters.crop(dem_b.copy(), aoi_b)
    except ValueError:
        raise ValueError(
            "Your AOI doesn't cover your DTM. Make sure your input data is valid."
        ) from ValueError
    try:
        if np.unique(dem_b) == 0 or np.unique(dem_b) == np.nan:
            raise ValueError(
                "Your DEM has no elevations. Make sure your input DEM is valid"
            ) from ValueError
    except ValueError:
        if np.unique(dem_b).all == 0 or np.unique(dem_b).all == np.nan:
            raise ValueError(
                "Your DEM has no elevations. Make sure your input DEM is valid"
            ) from ValueError
        dem_max = dem_b.max()
        dem_min = dem_b.min()

    aoi = vectors.read(aoi_o_path)
    aoi_b = vectors.read(aoi_b_path)

    # -- Define Resolution
    if output_resolution:
        output_resolution = int(output_resolution)
    else:
        x_res, y_res = dem_b.rio.resolution()
        output_resolution = int(np.round(abs((x_res) + abs(y_res) / 2)))

    # -- Reprojecting DEM
    # DEM will be used as input for SLOPE, ELEVATION, ASPECT and HYDRO layers. Also as reference for Rasterization of Geology Layer.
    dem_b = dem_b.astype(np.float32).rio.reproject(
        proj_crs, resolution=output_resolution, resampling=Resampling.bilinear
    )
    dem_b = dem_b.rio.write_nodata(FLOAT_NODATA)
    rasters.write(dem_b, AnyPath(tmp_dir) / "dem_b.tif")

    # --  LSI COMPUTATION ACCORDING TO METHOD
    # -- Check location (EUROPE or GLOBAL)
    if location == LocationType.EUROPE.value:
        LOGGER.info("-- LSI - LOCATION : EUROPE")

        # Define lithology
        lithology_path = str(DataPath.LITHOLOGY_PATH)

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
        elif landcover_name == LandcoverType.GLC.value:
            LOGGER.info(
                "-- LSI -- Method: Europe Refined layer based on Global Land Cover"
            )
            # -- Define path for LULC
            lulc_path = DataPath.GLC_PATH

        # -- 0. DEM
        LOGGER.info("-- Crop DEM")
        dem = rasters.crop(dem_b, aoi)

        # -- 1. Geology (to be erased for new Lithology dataset if ELSUS lithology is used)
        geology_tif = geology_raster(lithology_path, dem, aoi, proj_crs, tmp_dir)

        # -- 2. Landcover + Slope (calculations per zone in the AOI according to ELSUS)
        LOGGER.info("-- Defining physio zones for Europe Refined method")

        # -- Read Climate-physiographic zones
        physio_zones = vectors.read(physio_zones_path, window=aoi)
        physio_zones = physio_zones.to_crs(proj_crs)
        physio_zones_aoi = gpd.clip(physio_zones, aoi)

        vectors.write(physio_zones_aoi, os.path.join(tmp_dir, "physio_zone_AOI.shp"))
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

        total_zones = physio_zones_aoi["Zone"].count()
        LOGGER.info(
            str(
                "-- Produce the Landuse, Slope and Lithology raster for the LSI model | Amount of sub-zones to be processed: "
                + str(total_zones)
            )
        )

        # A counter i is set to separate polygons from the same Zone
        i = 0

        # The following loop crops the rasters by each physiographical zone in the AOI and stores the path of each
        # raster in the lists below
        landcover_list = []
        slope_list = []
        # lithology_list = []
        counter = enumerate(physio_zones_aoi.iterrows(), 1)
        for _, row in physio_zones_aoi.iterrows():
            # i += 1
            i = next(counter)[0]
            # The Climatic Zone
            zone = row["Zone"]
            db_file = zone_to_dbf[zone]

            LOGGER.info(
                str(
                    "-- -- Computing sub-zone type "
                    + str(zone)
                    + " | "
                    + str(i)
                    + "/"
                    + str(total_zones)
                )
            )

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
            )
            slope_list.append(slope_dir)

        # Mosaicing the rasters by zone available in the lists previously mentioned
        slope_mosaic_path, _ = mosaicing(slope_list, tmp_dir, "slope_weight")
        landcover_mosaic_path, _ = mosaicing(
            landcover_list, tmp_dir, "landcover_weight"
        )

        # Read the mosaiced rasters
        slope_tif = rasters.crop(slope_mosaic_path, aoi)
        landcover_tif = rasters.crop(landcover_mosaic_path, aoi)

        # Interpolate landcover and lithology layers to be equivalent to Slope layer
        landcover_tif = landcover_tif.interp_like(slope_tif, method="nearest")
        geology_tif = geology_tif.interp_like(slope_tif, method="nearest")

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
        # Define lithology
        lithology_path = str(DataPath.LITHOLOGY_PATH)

        # Define weights path
        global_final_weights_dbf_path = str(
            DataPath.WEIGHTS_GLOBAL_PATH / "Final_weights.dbf"
        )

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
        }
        lulc_path = landcover_path_dict[landcover_name]

        # 0. DEM
        LOGGER.info("-- DEM processing")

        dem = rasters.crop(dem_b.copy(), aoi)
        dem = dem.rio.reproject(
            proj_crs, resolution=output_resolution, resampling=Resampling.bilinear
        )
        dem = dem.rio.write_nodata(FLOAT_NODATA)

        # -- 1. Geology

        # Compute Geology layer
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            geology_layer = geology_raster(lithology_path, dem, aoi, proj_crs, tmp_dir)

        # -- 2. Slope
        slope_layer = slope_raster(dem_b, aoi, tmp_dir)

        # -- 3. Landcover
        lulc = rasters.read(lulc_path, window=aoi_b, as_type=np.float32)
        lulc = rasters.crop(lulc.copy(), aoi_b)
        lulc = rasters.collocate(dem_b, lulc.astype(np.float32), Resampling.nearest)
        lulc = rasters.crop(lulc.copy(), aoi_b)

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
        elevation_layer = elevation_raster(dem_b, aoi, tmp_dir)

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
        aspect_layer = aspect_raster(dem_b, aoi, tmp_dir)

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

    LOGGER.info("-- Writing LandslideSusceptibility.tif in memory")
    # Write in memory LSI with unclassified values
    rasters.write(lsi_tif, os.path.join(output_path, "LandslideSusceptibility.tif"))

    if (
        jenks_break
    ):  # Apply jenks only if user requires it (this option takes a longer time)
        # -- LSI reclassification: 1 to 5

        LOGGER.info("-- Reclassification of LSI in 5 vulnerability classes")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            lsi_reclass_tif = produce_a_reclass_arr(lsi_tif, location)

        LOGGER.info("-- Writing LandslideRisk.tif/shp in memory")
        # Post-processing
        lsi_tif_sieved, lsi_vector = raster_postprocess(
            lsi_reclass_tif, resolution=output_resolution
        )
        # Set 0 values to NaN
        lsi_tif_sieved = xr.where(lsi_tif_sieved <= 0, np.nan, lsi_tif_sieved)
        # There should not be 0 values for classes) (border effect)
        lsi_tif_sieved = lsi_tif_sieved.rio.write_crs(proj_crs, inplace=True)

        # Erase 255 no value
        lsi_vector = lsi_vector.drop(
            lsi_vector[lsi_vector["raster_val"].isin([0, 255])].index
        ).reset_index()
        vectors.write(lsi_vector, os.path.join(output_path, "LandslideRisk.shp"))

        # Write in memory
        rasters.write(lsi_tif_sieved, os.path.join(output_path, "LandslideRisk.tif"))

    LOGGER.info("-- Computing LSI statistics (FER_LR_av)")

    # Read GADM layer and overlay with AOI
    aoi_gadm = aoi_b
    aoi_gadm.geometry = aoi_gadm.geometry.buffer(GADM_BUFFER)
    with warnings.catch_warnings():  # For cases of polygons with more than 100 parts
        warnings.simplefilter("ignore")
        gadm = vectors.read(gadm_path, window=aoi_gadm)
    gadm = gadm.to_crs(proj_crs)
    gadm = gpd.clip(gadm, aoi)

    lsi_stats = compute_statistics(
        gadm, os.path.join(output_path, "LandslideSusceptibility.tif"), location
    )

    LOGGER.info("-- Writing LSI statistics in memory")

    # Write statistics in memory
    vectors.write(lsi_stats, os.path.join(output_path, "FER_LR_ave.shp"))

    if not temp:
        LOGGER.info("-- Deleting temporary files")
        shutil.rmtree(tmp_dir, ignore_errors=True)

    return
    # raise NotImplementedError
