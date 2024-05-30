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

from lsi.src.utils import (initialize_whitebox_tools, compute_flow_accumulation
                       , RoutingAlgorithm, np_to_xr, xr_to_gdf, classify_raster,
                       aspect)

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

    # aoi_path = input_dict.get(InputParameters.AOI_PATH.value)
    # location = input_dict.get(InputParameters.LOCATION.value)
    dem_name = input_dict.get(InputParameters.DEM_NAME.value)
    other_dem_path = input_dict.get(InputParameters.OTHER_DEM_PATH.value)
    # lithology_path = input_dict.get(InputParameters.LITHOLOGY_PATH.value)
    # landcover_name = input_dict.get(InputParameters.LANDCOVER_NAME.value)
    # weights_path = input_dict.get(InputParameters.WEIGHTS_PATH.value)
    # epsg_code = input_dict.get(InputParameters.REF_EPSG.value)
    # output_path = input_dict.get(InputParameters.OUTPUT_DIR.value)

    # Check if other_dem_path is needed
    if (dem_name == DemType.OTHER.value) and (other_dem_path is None):
        raise ValueError(f"Dem path is needed !")
    return  

def geology_raster(geology_dbf, litho_shp, dem, aoi, output_path):
    """
    
    """
    LOGGER.info(
        "-- Produce the Geology/Lithology raster for the LSI model"
    )
    if not os.path.exists(os.path.join(output_path, "geology_weight.tif")):
        proj_crs = aoi.crs
        litho_raster = rasters.rasterize( path_or_ds = dem
                                        , vector = litho_shp
                                        , value_field = "Rating"
                                        )
        litho_raster = litho_raster.fillna(997)
        litho_raster = rasters.crop(litho_raster, aoi)

        litho_gdf = xr_to_gdf(litho_raster
                        , proj_crs
                        , column_name = litho_raster.name
                        , column_rename = "Value")
        
        # -- JOIN with Geology_dbf
        geology_tif = litho_gdf.merge(geology_dbf, on="Value")
        geology_tif = geology_tif.set_index(["y", "x"]).Weights.to_xarray()
        geology_tif = geology_tif.rio.write_crs(litho_raster.rio.crs)
        geology_tif = rasters.crop(geology_tif, aoi)

        rasters.write(geology_tif, os.path.join(output_path, "geology_weight.tif"))

        return geology_tif
    else:
        return rasters.read(os.path.join(output_path, "geology_weight.tif"))

def slope_raster(slope_dbf, dem, aoi, output_path):
    """
    """
    LOGGER.info(
        "-- Produce the Slope raster for the LSI model"
    )
    if not os.path.exists(os.path.join(output_path, "slope_weight.tif")):
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

        rasters.write(slope_tif, os.path.join(output_path, "slope_weight.tif"))

        return slope_tif
    else:
        return rasters.read(os.path.join(output_path, "slope_weight.tif"))

def landuse_raster(landuse_dbf, lulc_path, dem, aoi, output_path):
    """
    """
    LOGGER.info(
    "-- Produce the Land use raster for the LSI model"
    )

    if not os.path.exists(os.path.join(output_path, "landcover_weight.tif")):
        proj_crs = aoi.crs
        aoi_b = geometry.buffer(aoi, buffer_n) # buffer (for LandUse + Hydro rasters)
        lulc = rasters.crop(lulc_path, aoi_b)
        lulc = rasters.collocate(dem, lulc, Resampling.nearest)

        # Reclassification of LULC for LSI calculation
        landcover_reclass = xr.where(lulc == 10, 3, lulc) # Tree cover -> Forest
        landcover_reclass = xr.where(landcover_reclass == 20, 3, landcover_reclass) # Shrubland -> Forest
        landcover_reclass = xr.where(landcover_reclass == 30, 2, landcover_reclass) # Grassland -> Grassland
        landcover_reclass = xr.where(landcover_reclass == 40, 1, landcover_reclass) # Cropland -> Arable land
        landcover_reclass = xr.where(landcover_reclass == 50, 5, landcover_reclass) # Built-up -> Urban areas
        landcover_reclass = xr.where(landcover_reclass == 60, 4, landcover_reclass) # Bare/Sparse vegetation -> Bare
        landcover_reclass = xr.where(landcover_reclass == 70, 997, landcover_reclass) # Snow and ice -> Not applicable
        landcover_reclass = xr.where(landcover_reclass == 80, 6, landcover_reclass) # Permanent water bodies -> Water areas
        landcover_reclass = xr.where(landcover_reclass == 90, 997, landcover_reclass) # Herbaceous wetland -> Not applicable
        landcover_reclass = xr.where(landcover_reclass == 95, 997, landcover_reclass) # Mangroves -> Not applicable
        landcover_reclass = xr.where(landcover_reclass == 100, 2, landcover_reclass) # Moss and lichen -> Grassland
        landcover_reclass = landcover_reclass.rio.write_crs(proj_crs, inplace=True)

        landcover_gdf = xr_to_gdf(landcover_reclass, proj_crs, landcover_reclass.name, "Value")
        lulc_tif = landcover_gdf.merge(landuse_dbf, on="Value")
        lulc_tif = lulc_tif.set_index(["y", "x"]).Weights.to_xarray()
        lulc_tif = lulc_tif.rio.write_crs(landcover_reclass.rio.crs)
        lulc_tif = rasters.crop(lulc_tif, aoi)

        rasters.write(lulc_tif, os.path.join(output_path, "landcover_weight.tif"))
        return lulc_tif
    else:
        return rasters.read(os.path.join(output_path, "landcover_weight.tif"))

def elevation_raster(elevation_dbf, dem, aoi, output_path):
    """
    """
    LOGGER.info(
        "-- Produce the Elevation raster for the LSI model"
    )
    if not os.path.exists(os.path.join(output_path, "elevation_weight.tif")):

        proj_crs = aoi.crs

        # Reclassify
        ELEVATION_STEPS = [0, 500, 600, 700, 800, 9000]
        ELEVATION_CLASSES = {1: f"{ELEVATION_STEPS[0]} - {ELEVATION_STEPS[1]}", 
                            2: f"{ELEVATION_STEPS[1]} - {ELEVATION_STEPS[2]}", 
                            3: f"{ELEVATION_STEPS[2]} - {ELEVATION_STEPS[3]}", 
                            4: f"{ELEVATION_STEPS[3]} - {ELEVATION_STEPS[4]}", 
                            5: f"{ELEVATION_STEPS[4]} - {ELEVATION_STEPS[5]}",
                            6: f"{ELEVATION_STEPS[5]}"
                            }
        elevation_name = "Value"
        elevation_arr = classify_raster(dem, ELEVATION_STEPS, ELEVATION_CLASSES)
        elevation_d = dem.copy(data=elevation_arr).astype(np.float32).rename(elevation_name)
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

def hydro_raster(hydro_dbf, dem_buff, aoi, ref_raster, tmp_dir, output_path):
    """
    Make raster of hydro_weights
    """

    LOGGER.info(
        "-- Produce the Distance to Hydro raster for the LSI model"
    )
    if not os.path.exists(os.path.join(output_path, "hydro_weight.tif")):
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
            os.path.join(tmp_dir,"dem_d.tif"), 
            compress="deflate", 
            predictor=1, 
            nodata=FLOAT_NODATA,
            dtype=np.float32
        )
        # -- Hydro processing

        wbt = initialize_whitebox_tools()

        # -- Fill pits
        wbt.fill_single_cell_pits(
                                os.path.join(tmp_dir,"dem_d.tif"),
                                os.path.join(tmp_dir,"filled_pits.tif"),
                                )
        # -- Fill depressions
        wbt.fill_depressions(
                                os.path.join(tmp_dir,"filled_pits.tif"), 
                                os.path.join(tmp_dir,"filled_depressions.tif"),
                            )
        # -- Flow accumulation
        compute_flow_accumulation(
                                os.path.join(tmp_dir,"filled_depressions.tif"),
                                os.path.join(tmp_dir,"flow_acc.tif"),
                                wbt,
                                RoutingAlgorithm.D8,
                                )
        
        flow_acc = rasters.read(os.path.join(tmp_dir,"flow_acc.tif"))
        # Thresholding the flow accumulation
        elevation_threshold = (dem_buff.max() - abs(dem_buff.min())).values
        flow_acc_thresh = xr.where(flow_acc > elevation_threshold, 1, 0)
        rasters.write(
                        flow_acc_thresh, 
                        os.path.join(tmp_dir,"flow_acc_thresh.tif"), 
                        compress="deflate", 
                        predictor=1, 
                        dtype=np.float32
                    )
        # Flow_acc raster to polyline
        wbt.raster_to_vector_lines(
                                    os.path.join(tmp_dir,'flow_acc_thresh.tif'), 
                                    os.path.join(tmp_dir,'flowacc_thresh_polyline.shp'))
        flowacc_thresh_lines = vectors.read(os.path.join(tmp_dir,"flowacc_thresh_polyline.shp"))
        flowacc_thresh_lines = rasters.rasterize(  
                                                path_or_ds = dem_b
                                                , vector = flowacc_thresh_lines
                                                , value_field = "VALUE"
                                                )
        rasters.write(flowacc_thresh_lines
                        ,os.path.join(tmp_dir,"flowacc_thresh_lines.tif")
                        , compress="deflate"
                        , predictor=1)
        
        # -- Euclidean Distance to River
        with rio.open(os.path.join(tmp_dir,"flowacc_thresh_lines.tif")) as src:
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

        # From UTM to LatLon
        dst_transform = ref_raster.rio.transform
        dst_crs = ref_raster.rio.crs

        hydro_tif = hydro_tif.rio.reproject(dst_crs = dst_crs,
                                            dst_transfom = dst_transform,
                                            resampling=Resampling.bilinear)
        hydro_tif = rasters.collocate(ref_raster, hydro_tif, Resampling.bilinear)

        # Write in memory
        rasters.write(hydro_tif,  os.path.join(output_path, "hydro_weight.tif"))
        return hydro_tif
    else:
        return rasters.read(os.path.join(output_path, "hydro_weight.tif"))

def aspect_raster(aspect_dbf, dem, aoi, output_path):
    """
    """
    LOGGER.info(
        "-- Produce the Aspect raster for the LSI model"
    )
    if not os.path.exists(os.path.join(output_path, "aspect_weight.tif")):
        proj_crs = aoi.crs

        aspect_tif = aspect(dem, proj_crs)

        aspect_tif_deg = np.round(aspect_tif * (180/np.pi)) + 180

        # Taking the maximum Azimuth as the Flat (360 -> -1) for classification purposes
        aspect_tif_deg = xr.where(aspect_tif_deg == aspect_tif_deg.max(), -1, aspect_tif_deg)

        # -- Classify from degrees [0, 359] Flat/North/Northeast/etc...
        ASPECT_STEPS = [-1, 0, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5]
        ASPECT_CLASSES = {   1: f"{ASPECT_STEPS[0]} - {ASPECT_STEPS[1]}", # Flat
                            2: f"{ASPECT_STEPS[1]} - {ASPECT_STEPS[2]}", # North
                            3: f"{ASPECT_STEPS[2]} - {ASPECT_STEPS[3]}", # Northeast
                            4: f"{ASPECT_STEPS[3]} - {ASPECT_STEPS[4]}", # East
                            5: f"{ASPECT_STEPS[4]} - {ASPECT_STEPS[5]}", # Southeast
                            6: f"{ASPECT_STEPS[5]} - {ASPECT_STEPS[6]}", # South
                            7: f"{ASPECT_STEPS[6]} - {ASPECT_STEPS[7]}", # Southwest
                            8: f"{ASPECT_STEPS[7]} - {ASPECT_STEPS[8]}", # West
                            9: f"{ASPECT_STEPS[8]} - {ASPECT_STEPS[9]}", # Northwest
                            10: f"{ASPECT_STEPS[9]}" # North
                            }
        aspect_class = classify_raster(aspect_tif_deg, ASPECT_STEPS, ASPECT_CLASSES)

        # -- Classify following the aspect classes from dbf
        # the following codes are not considered = {1: Flat, 2: North, 3: Northeast, 4: East}
        # as they already are within the class defined in Aspect_dbf

        # Transform to Aspect_dbf scale
        aspect_reclass = xr.where(aspect_class == 5, 4, aspect_class) # Southeast
        aspect_reclass = xr.where(aspect_reclass == 6, 5, aspect_reclass) # South
        aspect_reclass = xr.where(aspect_reclass == 7, 5, aspect_reclass) # Southwest
        aspect_reclass = xr.where(aspect_reclass == 8, 3, aspect_reclass) # West
        aspect_reclass = xr.where(aspect_reclass == 9, 3, aspect_reclass) # Northeast
        aspect_reclass = xr.where(aspect_reclass == 10, 2, aspect_reclass) # North

        aspect_reclass_xr = np_to_xr(aspect_reclass, dem, proj_crs)
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



    # Store landcover path in a variable depending on LOCATION
        # -- Check location (EUROPE or OUTSIDE of Europe)
    if location == LocationType.EUROPE.value:
            # -- Dict that store landcover name and landcover path
        landcover_path_dict = {
            LandcoverType.CLC.value: DataPath.CLC_PATH,
            LandcoverType.ESAWC.value: DataPath.ESAWC_PATH,
        }

        lulc_path = landcover_path_dict[landcover_name]

    if location == LocationType.GLOBAL.value:
            # -- Dict that store landcover name and landcover path
        landcover_path_dict = {
            LandcoverType.ESAWC.value: DataPath.ESAWC_PATH,
        }

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

    # -- Define Weights dbfs paths:
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


    #0. DEM
    LOGGER.info(
        "-- Reading/cropping DEM"
    )
    aoi_b = geometry.buffer(aoi, buffer_n) # buffer (for LandUse + Hydro rasters)
    dem_buff = rasters.crop(dem_path, aoi_b)
    dem = rasters.crop(dem_buff, aoi)

    # -- 1. Geology
    # Reading geology database, clip to aoi and reproject to proj_crs
    litho_db = gpd.read_file(lithology_path, driver='FileGDB', mask = aoi)
    aoi_db = aoi.to_crs(litho_db.crs)
    litho_shp = gpd.clip(litho_db, aoi_db)
    litho_shp = litho_shp.to_crs(proj_crs)

    geology_dbf = gpd.read_file(geology_dbf_path)
    # momentaneous line to add Not Applicable class
    geology_dbf.loc[len(geology_dbf)] = ["Not Applicable", 997, 0.0, 0.0, None]
    geology_layer = geology_raster(geology_dbf, litho_shp, dem, aoi, output_path)

    # -- 2. Slope
    slope_dbf = gpd.read_file(slope_dbf_path)
    slope_layer = slope_raster(slope_dbf, dem, aoi, output_path)
    
    # -- 3. Landcover
    landuse_dbf = gpd.read_file(landuse_dbf_path)
    landuse_dbf.loc[len(landuse_dbf)] = ["Not Applicable", 997, 0.0, 0.0, None]
    # Landcover
    landuse_layer = landuse_raster(landuse_dbf, lulc_path, dem, aoi_b, output_path)
    
    # -- 4. Elevation
    elevation_dbf = gpd.read_file(elevation_dbf_path)
    elevation_layer = elevation_raster(elevation_dbf, dem, aoi, output_path)

    # -- 5. Hydro
    hydro_dbf = gpd.read_file(hydro_dbf_path)
    hydro_layer = hydro_raster(hydro_dbf, dem_buff, aoi, elevation_layer, tmp_dir, output_path)

    # -- 6. Aspect 
    aspect_dbf = gpd.read_file(aspect_dbf_path)
    
    aspect_layer = aspect_raster(aspect_dbf, dem, aoi, output_path)


    # -- Final weights computing
    fw_dbf = gpd.read_file(final_weights_dbf_path)

    # Extracting final weights
    slope_weights = fw_dbf[fw_dbf.Factors == "Slope"].Weights.iloc[0]
    geology_weights = fw_dbf[fw_dbf.Factors == "Geology"].Weights.iloc[0]
    aspect_weights = fw_dbf[fw_dbf.Factors == "Slope aspect"].Weights.iloc[0]
    elevation_weights = fw_dbf[fw_dbf.Factors == "Elevation"].Weights.iloc[0]
    hydro_weights = fw_dbf[fw_dbf.Factors == "Distance from river"].Weights.iloc[0]
    landuse_weights = fw_dbf[fw_dbf.Factors == "Land use"].Weights.iloc[0]

    # Final weight

    lsi_tif = (slope_layer*float(slope_weights) 
            + geology_layer*float(geology_weights) 
            + elevation_layer*float(elevation_weights) 
            + aspect_layer*float(aspect_weights)
            + landuse_layer*float(landuse_weights)
            + hydro_layer*float(hydro_weights)
            )
    
    # Write in memory
    rasters.write(lsi_tif, os.path.join(output_path, "lsi.tif"))

    return
    # raise NotImplementedError


# def lsi_core(input_dict: dict) -> str:
#     """
#     TODO: Complete arguments and dosctring
#     """

#     # 0. Prepare Rasters and Database files of Weights
#     lsi_compute(input_dict)

#     return
#     # raise NotImplementedError
