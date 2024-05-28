""" lsi """
import logging
from enum import unique
import os

import numpy as np

from sertit.misc import ListEnum
from sertit import rasters, vectors, geometry, AnyPath
from sertit.rasters import FLOAT_NODATA
from sertit.unistra import get_geodatastore, s3_env

import rasterio as rio
from rasterio.enums import Resampling

import xarray as xr

import geopandas as gpd

from scipy.ndimage import distance_transform_edt

from src.utils import (initialize_whitebox_tools, compute_flow_accumulation
                       , RoutingAlgorithm, np_to_xr, xr_to_gdf, classify_raster)

# from lsi import LOGGER_NAME # on the meantime to solve the acces to lsi.py

DEBUG = False
LOGGING_FORMAT = "%(asctime)s - [%(levelname)s] - %(message)s"
LOGGER = logging.getLogger("LSI")  # on the meantime to solve the acces to lsi.py


def geodatastore(ftep=False):
    """
    This function returns the root path to the geo data store (DEM, Soil database...).
    Args:
        ftep: If True, the path to the s3 bucket for the ftep platform is returned. Else, get_geodatastore from sertit utils module is called.

    Returns:
        The path to the geo data store.

    """
    if ftep:
        return AnyPath("s3://eo4sdg-data") #TODO
    else:
        return get_geodatastore()


class DataPath:
    GLOBAL_DIR = None

    @classmethod
    def load_paths(cls, ftep=False):
        cls.GLOBAL_DIR = geodatastore(ftep) / "GLOBAL"
        cls.ESAWC_PATH = (
            cls.GLOBAL_DIR 
            / "ESA_WorldCover" 
            / "2021" 
            / "ESA_WorldCover_10m_2021.vrt")
        cls.CLC_PATH = (
            cls.GLOBAL_DIR
            / "Corine_Land_Cover"
            / "CLC_2018"
            / "clc2018_clc2018_v2018_20_raster100m"
            / "CLC2018_CLC2018_V2018_20.tif"
            )

#        cls.EUDEM_PATH = cls.GLOBAL_DIR / "EUDEM_v2" / "eudem_dem_3035_europe.tif"
        cls.SRTM30_PATH = cls.GLOBAL_DIR / "SRTM_30m_v4" / "index.vrt"
        cls.COPDEM30_PATH = (
            cls.GLOBAL_DIR / "COPDEM_30m" / "COPDEM_30m.vrt"
        )
        cls.FABDEM_PATH = cls.GLOBAL_DIR / "FABDEM" / "FABDEM.vrt"
        


@unique
class InputParameters(ListEnum):
    """
    List of the input parameters
    """

    AOI_PATH = "aoi_path"
    LOCATION = "location"
    DEM_NAME = "dem_name"
    OTHER_DEM_PATH = "other_dem_path"
    LITHOLOGY_PATH = "lithology_path"
    LANDCOVER_NAME = "landcover_name"
    WEIGHTS_PATH = "weights_path"
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

def check_parameters(input_dict: dir) -> None:
    """
     Check if parameters values are valid
    Args:
        input_dict (dict) : dict with parameters values

    Returns:

    """
    # --- Extract parameters ---

    aoi_path = input_dict.get(InputParameters.AOI_PATH.value)
    location = input_dict.get(InputParameters.LOCATION.value)
    dem_name = input_dict.get(InputParameters.DEM_NAME.value)
    other_dem_path = input_dict.get(InputParameters.OTHER_DEM_PATH.value)
    lithology_path = input_dict.get(InputParameters.LITHOLOGY_PATH.value)
    landcover_name = input_dict.get(InputParameters.LANDCOVER_NAME.value)
    weights_path = input_dict.get(InputParameters.WEIGHTS_PATH.value)
    epsg_code = input_dict.get(InputParameters.REF_EPSG.value)
    output_path = input_dict.get(InputParameters.OUTPUT_DIR.value)

    # Check if other_dem_path is needed
    if (dem_name == DemType.OTHER.value) and (other_dem_path is None):
        raise ValueError(f"Dem path is needed !")
    return  

def geology_raster(geology_dbf, litho_shp, dem, aoi, output_path):
    """
    
    """
    proj_crs = aoi.crs
    litho_raster = rasters.rasterize( path_or_ds = dem
                                    , vector = litho_shp
                                    , value_field = "Rating"
                                    )
    litho_shp_raster = litho_shp_raster.fillna(997)
    litho_shp_raster = rasters.crop(litho_shp_raster, aoi)

    litho_gdf = xr_to_gdf(litho_shp_raster
                      , proj_crs
                      , column_name = litho_shp_raster.name
                      , column_rename = "Value")
    
    # -- JOIN with Geology_dbf
    geology_tif = litho_gdf.merge(geology_dbf, on="Value")
    geology_tif = geology_tif.set_index(["y", "x"]).Weights.to_xarray()
    geology_tif = geology_tif.rio.write_crs(litho_shp_raster.rio.crs)
    geology_tif = rasters.crop(geology_tif, aoi)

    rasters.write(geology_tif, output_path / "geology_weight.tif")

    return geology_tif

def slope_raster(slope_dbf, dem, aoi, output_path):
    """
    """
    proj_crs = aoi.crs
    slope = rasters.slope(dem, in_rad = False)

    # -- Classify
    SLOPE_STEPS = [0, 2, 5, 15, 35, 90]
    SLOPE_CLASSES = {1: f"{SLOPE_STEPS[0]} - {SLOPE_STEPS[1]}", 
                    2: f"{SLOPE_STEPS[1]} - {SLOPE_STEPS[2]}", 
                    3: f"{SLOPE_STEPS[2]} - {SLOPE_STEPS[3]}", 
                    4: f"{SLOPE_STEPS[3]} - {SLOPE_STEPS[4]}", 
                    5: f"{SLOPE_STEPS[4]} - {SLOPE_STEPS[5]}",
                    6: f"{SLOPE_STEPS[5]}"
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

    rasters.write(slope_tif, output_path / "slope_weight.tif")

    return slope_tif

def hydro_raster(hydro_dbf, dem_buff, aoi, tmp_dir, output_path):
    """
    Make raster of hydro_weights
    """
    # -- Pre-process DEM
        ## Prepare DEM
    proj_crs = aoi.crs
    utm_zone = aoi.estimate_utm_crs()
        # no data 
    dem_b = xr.where(dem_buff <= 0, FLOAT_NODATA, dem_buff)
    dem_b = dem_b.rio.write_nodata(FLOAT_NODATA, encoded=True)
    dem_b = dem_b.rio.write_crs(proj_crs)
        # reproject
    dem_b = dem_b.rio.reproject(dst_crs=utm_zone,
                                nodata=rasters.FLOAT_NODATA,
                                resampling=Resampling.bilinear,
                                )
    # resolution
    x_res, y_res = dem_b.rio.resolution()
    dem_resolution = (abs(x_res) + abs(y_res)) / 2
    # Write DEM in memory
    rasters.write(
        dem_b, 
        tmp_dir / "dem_d.tif", 
        compress="deflate", 
        predictor=1, 
        nodata=FLOAT_NODATA,
        dtype=np.float32
    )
    # -- Hydro processing

    wbt = initialize_whitebox_tools()

    # -- Fill pits
    wbt.fill_single_cell_pits(
                            tmp_dir / "dem_d.tif",
                            tmp_dir / "filled_pits.tif",
                            )
    # -- Fill depressions
    wbt.fill_depressions(
                            tmp_dir / "filled_pits.tif", 
                            tmp_dir / "filled_depressions.tif",
                        )
    # -- Flow accumulation
    compute_flow_accumulation(
                            tmp_dir / "filled_depressions.tif",
                            tmp_dir / "flow_acc.tif",
                            wbt,
                            RoutingAlgorithm.D8,
                            )
    
    flow_acc = rasters.read(tmp_dir / "flow_acc.tif")
    # Thresholding the flow accumulation
    elevation_threshold = (dem_buff.max() - abs(dem_buff.min())).values
    flow_acc_thresh = xr.where(flow_acc > elevation_threshold, 1, 0)
    rasters.write(
                    flow_acc_thresh, 
                    tmp_dir / "flow_acc_thresh.tif", 
                    compress="deflate", 
                    predictor=1, 
                    dtype=np.float32
                )
    # Flow_acc raster to polyline
    wbt.raster_to_vector_lines(
                                tmp_dir / 'flow_acc_thresh.tif', 
                                tmp_dir / 'flowacc_thresh_polyline.shp')
    flowacc_thresh_lines = vectors.read(tmp_dir / "flowacc_thresh_polyline.shp")
    flowacc_thresh_lines = rasters.rasterize(  
                                            path_or_ds = dem_b
                                            , vector = flowacc_thresh_lines
                                            , value_field = "VALUE"
                                            )
    rasters.write(flowacc_thresh_lines
                    , tmp_dir / "flowacc_thresh_lines.tif"
                    , compress="deflate"
                    , predictor=1)
    
    # -- Euclidean Distance to River
    with rio.open(tmp_dir / "flowacc_thresh_lines.tif") as src:
        river_streams = src.read(1)
        profile = src.profile
        transform = src.transform
        nodata = src.nodata

    # Transform raster values
    # Invert the values so that river cells become background
    river_streams_inverted = np.where(river_streams == nodata, 1, 0)

    # Euclidean distance
    euclidean_distance = distance_transform_edt(river_streams_inverted
                                            # , sampling=profile['transform'][0]
                                            )

    euclidean_distance_xr = np_to_xr(euclidean_distance, dem_b, proj_crs)
        # transform from pixel to meters
    euclidean_distance_xr = euclidean_distance_xr * dem_resolution
    # -- Reclassify
    ED_STEPS = [0, 100, 200, 300, 400, 20000] #5000
    ED_CLASSES = {   1: f"{ED_STEPS[0]} - {ED_STEPS[1]}", #
                    2: f"{ED_STEPS[1]} - {ED_STEPS[2]}", #
                    3: f"{ED_STEPS[2]} - {ED_STEPS[3]}", #
                    4: f"{ED_STEPS[3]} - {ED_STEPS[4]}", #
                    5: f"{ED_STEPS[4]} - {ED_STEPS[5]}", #
                    6: f"{ED_STEPS[5]}"
                        }
    ed_class = classify_raster(euclidean_distance_xr, ED_STEPS, ED_CLASSES)
    ed_reclass = rasters.crop(np_to_xr(ed_class, euclidean_distance_xr, proj_crs), aoi) # el proj_crs jode todo

    # -- JOIN with Hydro.dbf
    hydro_gdf = xr_to_gdf(ed_reclass, ed_reclass.rio.crs)
    hydro_tif = hydro_gdf.merge(hydro_dbf, on="Value")
    hydro_tif = hydro_tif.set_index(["y", "x"]).Weights.to_xarray()
    hydro_tif = hydro_tif.rio.write_crs(ed_reclass.rio.crs)

    rasters.write(hydro_tif,  output_path / "hydro_weight.tif")
    return hydro_tif

def make_raster_list(input_dict):
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
    # --- Extract parameters ---

    aoi_path = input_dict.get(InputParameters.AOI_PATH.value)
    location = input_dict.get(InputParameters.LOCATION.value)
    dem_name = input_dict.get(InputParameters.DEM_NAME.value)
    other_dem_path = input_dict.get(InputParameters.OTHER_DEM_PATH.value)
    lithology_path = input_dict.get(InputParameters.LITHOLOGY_PATH.value)
    landcover_name = input_dict.get(InputParameters.LANDCOVER_NAME.value)
    weights_path = input_dict.get(InputParameters.WEIGHTS_PATH.value)
    epsg_code = input_dict.get(InputParameters.REF_EPSG.value)
    output_path = input_dict.get(InputParameters.OUTPUT_DIR.value)

    # -- Create temp_dir if not exist
    tmp_dir = os.path.join(output_path, "temp_dir")
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # -- Dict that store landcover name and landcover path
    landcover_path_dict = {
        LandcoverType.CLC.value: DataPath.CLC_PATH,
        LandcoverType.ESAWC.value: DataPath.ESAWC_PATH,
    }

    # -- Store landcover path in a variable
    lulc_path = landcover_path_dict[landcover_name]

    # -- Dict that store landcover name and landcover path
    landcover_path_dict = {
        LandcoverType.CLC.value: DataPath.CLC_PATH,
        LandcoverType.ESAWC.value: DataPath.ESAWC_PATH,
    }

    # -- Store landcover path in a variable
    lulc_path = landcover_path_dict[landcover_name]

    # -- Dict that store dem name and dem path
    dem_path_dict = {DemType.COPDEM_30.value: DataPath.COPDEM30_PATH,
                     DemType.FABDEM.value: DataPath.FABDEM_PATH,
                     DemType.SRTM.value: DataPath.SRTM30_PATH,
                     DemType.OTHER.value: other_dem_path}
    # Store DEM path in a variable
    dem_path = dem_path_dict[dem_name]

    # -- Read AOI
    aoi = vectors.read(aoi_path)
    if epsg_code:
        proj_crs = epsg_code
    else:
        proj_crs = aoi.crs

    raster_dict = {}

    # Define Weights dbfs paths:
    geology_dbf_path = os.path.join(weights_path, "Geology.dbf")
    slope_dbf_path = os.path.join(weights_path, "Slope.dbf")
    elevation_dbf_path = os.path.join(weights_path, "Elevation.dbf")
    aspect_dbf_path = os.path.join(weights_path, "Aspect.dbf")
    landuse_dbf_path = os.path.join(weights_path, "Land use.dbf")
    hydro_dbf_path = os.path.join(weights_path, "Hydro.dbf")
    final_weights_dbf_path = os.path.join(weights_path, "Final_weights.dbf")

    # Define folder for temporal files
    tmp_dir = os.path.join(output_path, "temp_dir")
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
 

    # -- Check location (EUROPE or OUTSIDE of Europe)
    if location == LocationType.EUROPE.value:
        print("do european stuff")

    if location == LocationType.GLOBAL.value:
        print("do gLoBaL CitIZeN stuff")

        #0. DEM
        aoi_b = geometry.buffer(aoi, 0.1) # buffer (for LandUse + Hydro rasters)
        dem_buff = rasters.crop(dem_path, aoi_b)
        dem = rasters.crop(dem_buff, aoi)

        # -- 1. Geology
        litho_shp = vectors.read(vector_path = lithology_path,
                         crs = proj_crs,
                         window = aoi)
        litho_shp = gpd.clip(litho_shp, aoi)

        geology_dbf = gpd.read_file(geology_dbf_path)
        # momentaneous line to add Not Applicable class
        geology_dbf.loc[len(geology_dbf)] = ["Not Applicable", 997, 0.0, 0.0, None]

        geology_layer = geology_raster(geology_dbf, litho_shp, dem, aoi, output_path)

        # -- 2. Slope
        slope_dbf = gpd.read_file(slope_dbf_path)
        slope_layer = slope_raster(slope_dbf, dem, aoi, output_path)
        # -- 3. Landcover
        # -- 4. Elevation
                

        # -- 5. Hydro
        hydro_dbf = gpd.read_file(hydro_dbf_path)
        hydro_layer = hydro_raster(hydro_dbf, dem_buff, aoi, tmp_dir, output_path)

        # -- 6. Aspect 
        
        # Final
        raster_dict = {
        }

def lsi_core(input_dict: dict) -> str:
    """
    TODO: Complete arguments and dosctring
    """

    # 0. Prepare Rasters and Database files of Weights
    make_raster_list(input_dict)


    # 1. GEOLOGY

    # 2. SLOPE (DEGREES)

    # 3. LANDCOVER

    # 4. ELEVATION

    # 5. DISTANCE TO HYDRO

    # 6. ASPECT

    # 7. FINAL WEIGTHS

    # 8. LSI RASTER CALCULATION





    return
    # raise NotImplementedError
