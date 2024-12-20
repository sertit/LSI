# This file is part of LSI.
# LSI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# LSI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with LSI. If not, see <https://www.gnu.org/licenses/>.
"""Lsi calculator"""

import logging
import os
import warnings

import geopandas as gpd
import numpy as np
import rasterio as rio
import whitebox_workflows as wbw
import xarray as xr
from rasterio.enums import Resampling
from scipy.ndimage import distance_transform_edt
from sertit import AnyPath, geometry, rasters, vectors
from sertit.rasters import FLOAT_NODATA

from lsi.src.reclass import classify_raster, reclass_landcover, reclass_landcover_elsus
from lsi.src.utils import aspect, np_to_xr, xr_to_gdf

LOGGING_FORMAT = "%(asctime)s - [%(levelname)s] - %(message)s"
LOGGER = logging.getLogger("LSI")

VERY_BIG_BUFFER = 5000
BIG_BUFFER = 1500
REGULAR_BUFFER = 1000
SMALL_BUFFER = 100


def join_weights(xr_raster, weights_dbf, out_crs, weight_column="Value"):
    """
    Args:
        xr_raster: Raster in Xarray format
        weights_dbf: weights database file in GeoDataFrame format
        out_crs: Coordinate reference system in String format
        weight_column:
    """

    raster_gdf = xr_to_gdf(
        xr_raster, out_crs, column_name=xr_raster.name, column_rename=weight_column
    )
    raster_weighted = raster_gdf.merge(weights_dbf, on=weight_column)
    try:
        raster_weighted = raster_weighted.set_index(["y", "x"]).Weights.to_xarray()
    except AttributeError:  # ELSUS method uses Weight instead of Weights
        raster_weighted = raster_weighted.set_index(["y", "x"]).Weight.to_xarray()

    return raster_weighted


# --- GLOBAL LSI method functions


def geology_raster(lithology_path, dem, aoi, proj_crs, output_path):
    """ """
    LOGGER.info("-- Produce the Geology/Lithology raster for the LSI model")
    if not os.path.exists(os.path.join(output_path, "geology_weight.tif")):
        # Reading geology database, clip to aoi and reproject to proj_crs
        with (
            warnings.catch_warnings()
        ):  # For cases of polygons with more than 100 parts
            warnings.simplefilter("ignore")
            litho_db = vectors.read(lithology_path, window=aoi)
        aoi_db = aoi.to_crs(litho_db.crs)
        litho_shp = gpd.clip(litho_db, aoi_db)
        litho_shp = litho_shp.to_crs(proj_crs)

        # Rasterize Lithology dataset
        litho_raster = rasters.rasterize(dem, litho_shp, value_field="Rating")

        # litho_raster = litho_raster.fillna(997)
        litho_raster = rasters.crop(litho_raster, aoi)

        def geology_transform(x):
            values = [1, 2, 3, 4, 5]
            weights = [0.062376, 0.098573, 0.16105, 0.261788, 0.416212]
            if x == float(values[0]):
                out_value = weights[0]
            elif x == float(values[1]):
                out_value = weights[1]
            elif x == float(values[2]):
                out_value = weights[2]
            elif x == float(values[3]):
                out_value = weights[3]
            elif x == float(values[4]):
                out_value = weights[4]
            else:
                out_value = float(0)
            return out_value

        geology_tif = xr.apply_ufunc(
            geology_transform, litho_raster, vectorize=True
        ).astype(np.float32)
        geology_tif = geology_tif.copy(data=geology_tif).rename("geology_weight")
        geology_tif = rasters.crop(geology_tif, aoi)

        rasters.write(geology_tif, os.path.join(output_path, "geology_weight.tif"))

        return geology_tif
    else:
        return rasters.read(os.path.join(output_path, "geology_weight.tif"))


def slope_raster(dem_b, aoi, output_path):
    """ """
    LOGGER.info("-- Produce the Slope raster for the LSI model")
    if not os.path.exists(os.path.join(output_path, "slope_weight.tif")):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            slope = rasters.slope(dem_b.astype(np.float32), in_rad=False)
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
        slope_arr = classify_raster(slope, SLOPE_STEPS, SLOPE_CLASSES).astype(
            np.float32
        )
        slope_d = slope.copy(data=slope_arr).rename(slope_name)
        slope_d.attrs["long_name"] = slope_name

        def slope_transform(x):
            values = [1, 2, 3, 4, 5]
            weights = [0.06, 0.1, 0.16, 0.26, 0.42]
            if x == float(values[0]):
                out_value = weights[0]
            elif x == float(values[1]):
                out_value = weights[1]
            elif x == float(values[2]):
                out_value = weights[2]
            elif x == float(values[3]):
                out_value = weights[3]
            elif x == float(values[4]):
                out_value = weights[4]
            else:
                out_value = float(0)
            return out_value

        slope_tif = xr.apply_ufunc(slope_transform, slope_d, vectorize=True).astype(
            np.float32
        )
        slope_tif = slope_tif.copy(data=slope_tif).rename("slope_weight")
        slope_tif = rasters.crop(slope_tif, aoi)

        rasters.write(slope_tif, os.path.join(output_path, "slope_weight.tif"))

        return slope_tif
    else:
        return rasters.read(os.path.join(output_path, "slope_weight.tif"))


def landcover_raster(
    lulc,
    landcover_name,
    aoi,
    proj_crs,
    output_resolution,
    output_path,
):
    """ """
    LOGGER.info("-- Produce the Land use raster for the LSI model")

    if not os.path.exists(os.path.join(output_path, "landcover_weight.tif")):
        # Reclassification of LULC for LSI calculation
        landcover_reclass = reclass_landcover(lulc, landcover_name)
        landcover_reclass = landcover_reclass.rio.write_crs(lulc.rio.crs, inplace=True)

        def landuse_transform(x):
            values = [1, 2, 3, 4, 5, 6]
            weights = [0.37, 0.25, 0.15, 0.09, 0.09, 0.05]
            if x == float(values[0]):
                out_value = weights[0]
            elif x == float(values[1]):
                out_value = weights[1]
            elif x == float(values[2]):
                out_value = weights[2]
            elif x == float(values[3]):
                out_value = weights[3]
            elif x == float(values[4]):
                out_value = weights[4]
            elif x == float(values[5]):
                out_value = weights[5]
            else:
                out_value = float(0)
            return out_value

        lulc_tif = xr.apply_ufunc(
            landuse_transform, landcover_reclass, vectorize=True
        ).astype(np.float32)
        lulc_tif = lulc_tif.copy(data=lulc_tif).rename("landcover_weight")

        # Reproject to proj_crs and resolution
        lulc_tif = lulc_tif.rio.reproject(
            proj_crs, resolution=output_resolution, resampling=Resampling.bilinear
        )
        lulc_tif = rasters.crop(lulc_tif, aoi)

        rasters.write(lulc_tif, os.path.join(output_path, "landcover_weight.tif"))
        return lulc_tif
    else:
        return rasters.read(os.path.join(output_path, "landcover_weight.tif"))


def elevation_raster(dem_b, aoi, output_path):
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
        elevation_arr = classify_raster(
            dem_b, ELEVATION_STEPS, ELEVATION_CLASSES
        ).astype(np.float32)
        elevation_d = (
            dem_b.copy(data=elevation_arr).astype(np.float32).rename(elevation_name)
        )
        elevation_d.attrs["long_name"] = elevation_name

        def elevation_transform(x):
            values = [1, 2, 3, 4, 5]
            weights = [0.06, 0.1, 0.16, 0.26, 0.42]
            if x == float(values[0]):
                out_value = weights[0]
            elif x == float(values[1]):
                out_value = weights[1]
            elif x == float(values[2]):
                out_value = weights[2]
            elif x == float(values[3]):
                out_value = weights[3]
            elif x == float(values[4]):
                out_value = weights[4]
            else:
                out_value = float(0)
            return out_value

        elevation_tif = xr.apply_ufunc(
            elevation_transform, elevation_d, vectorize=True
        ).astype(np.float32)
        elevation_tif = elevation_tif.copy(data=elevation_tif).rename(
            "elevation_weight"
        )

        elevation_tif = rasters.crop(elevation_tif, aoi)

        rasters.write(elevation_tif, os.path.join(output_path, "elevation_weight.tif"))

        return elevation_tif
    else:
        return rasters.read(os.path.join(output_path, "elevation_weight.tif"))


def hydro_raster_wbw(
    dem_buff, aoi, proj_crs, dem_max, dem_min, output_resolution, tmp_dir
):
    """
    Make raster of hydro_weights
    """
    LOGGER.info("-- Produce the Distance to Hydro raster for the LSI model")
    if not os.path.exists(os.path.join(tmp_dir, "hydro_weight.tif")):
        LOGGER.info("-- -- Preparing the DEM for hydro analysis")

        # Prepare DEM
        aoi_b = aoi.geometry.buffer(REGULAR_BUFFER)
        # reproject
        with (
            warnings.catch_warnings()
        ):  # Catching -> The nodata value (3.402823466e+38)
            warnings.simplefilter("ignore")
            dem_b = dem_buff.astype(np.float32).rio.reproject(
                dst_crs=proj_crs,
                nodata=rasters.FLOAT_NODATA,
                resampling=Resampling.bilinear,
            )
            # no data
            dem_b = xr.where(dem_buff <= -700, FLOAT_NODATA, dem_buff)
            dem_b = dem_b.rio.write_crs(proj_crs)

            # reproject
            dem_b = dem_b.rio.reproject(
                dst_crs=proj_crs,
                resampling=Resampling.bilinear,
            )
            dem_b = rasters.crop(dem_b, aoi_b)

        # Write DEM in memory for Whitebox to work
        rasters.write(
            dem_b,
            os.path.join(tmp_dir, "dem_d.tif"),
            compress="deflate",
            predictor=1,
            nodata=FLOAT_NODATA,
            dtype=np.float32,
        )
        # -- Hydro processing
        wbe = wbw.WbEnvironment()

        LOGGER.info("-- -- Preparing the DEM: Filling Pits")
        # -- Fill pits
        filled_pits = wbe.fill_pits(wbe.read_raster(os.path.join(tmp_dir, "dem_d.tif")))

        LOGGER.info("-- -- Preparing the DEM: Filling Depressions")
        # Write and read using WbW
        wbe.write_raster(filled_pits, os.path.join(tmp_dir, "filled_pits.tif"))
        wbe.read_raster(os.path.join(tmp_dir, "filled_pits.tif"))
        # -- Fill depressions
        filled_depressions = wbe.fill_depressions(filled_pits)

        # -- Flow accumulation
        LOGGER.info("-- -- Compute Flow Accumulation")

        flow_acc = wbe.fd8_flow_accum(filled_depressions, out_type="cells")

        # Write flow_acc so it can be read with sertit_utils for the following steps
        wbe.write_raster(
            flow_acc, os.path.join(tmp_dir, "flow_acc.tif"), compress=False
        )

        flow_acc = rasters.read(os.path.join(tmp_dir, "flow_acc.tif"))

        # Thresholding the flow accumulation
        elevation_threshold = (dem_max - abs(dem_min)).values
        flow_acc_thresh = xr.where(flow_acc > elevation_threshold, 1, 0)

        # Write to have a Sertit Raster for next step
        rasters.write(
            flow_acc_thresh,
            os.path.join(tmp_dir, "flow_acc_thresh.tif"),
            compress="deflate",
            predictor=1,
            dtype=np.float32,
        )

        flow_acc_thresh = wbe.read_raster(os.path.join(tmp_dir, "flow_acc_thresh.tif"))

        # Flow_acc raster to polyline
        flowacc_thresh_lines = wbe.raster_to_vector_lines(flow_acc_thresh)

        # Write vector in temporal directory
        wbe.write_vector(
            flowacc_thresh_lines, os.path.join(tmp_dir, "flowacc_thresh_lines.shp")
        )
        # Read with Sertit
        flowacc_thresh_lines = vectors.read(
            os.path.join(tmp_dir, "flowacc_thresh_lines.shp"), crs=proj_crs
        )

        # Rasterization
        flowacc_thresh_lines = flowacc_thresh_lines.set_crs(proj_crs)
        flowacc_thresh_lines = flowacc_thresh_lines.to_crs(aoi.crs)

        flowacc_thresh_lines = rasters.rasterize(
            dem_b, flowacc_thresh_lines, value_field="VALUE"
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
        euclidean_distance = distance_transform_edt(river_streams_inverted)

        euclidean_distance_xr = np_to_xr(euclidean_distance, dem_b)
        euclidean_distance_xr = (
            euclidean_distance_xr.copy(data=euclidean_distance_xr)
            .astype(np.float32)
            .rename("euclidean_distance")
        )

        # transform from pixel to meters

        LOGGER.info("-- -- Distance to rivers classification")
        euclidean_distance_xr = euclidean_distance_xr * int(output_resolution)

        # -- Reclassify
        ED_STEPS = [0, 100, 200, 300, 400, 20000]
        ED_CLASSES = {
            1: f"{ED_STEPS[0]} - {ED_STEPS[1]}",
            2: f"{ED_STEPS[1]} - {ED_STEPS[2]}",
            3: f"{ED_STEPS[2]} - {ED_STEPS[3]}",
            4: f"{ED_STEPS[3]} - {ED_STEPS[4]}",
            5: f"{ED_STEPS[4]} - {ED_STEPS[5]}",
            6: f"{ED_STEPS[5]}",
        }
        ed_class = classify_raster(euclidean_distance_xr, ED_STEPS, ED_CLASSES)
        ed_class = np_to_xr(ed_class, euclidean_distance_xr)
        ed_class = (
            ed_class.copy(data=ed_class).astype(np.int32).rename("euclidean_distance")
        )

        ed_reclass = rasters.crop(ed_class, aoi)
        # -- JOIN with Hydro.dbf

        def hydro_transform(x):
            values = [1, 2, 3, 4, 5]
            weights = [0.42, 0.26, 0.16, 0.1, 0.06]
            if x == float(values[0]):
                out_value = weights[0]
            elif x == float(values[1]):
                out_value = weights[1]
            elif x == float(values[2]):
                out_value = weights[2]
            elif x == float(values[3]):
                out_value = weights[3]
            elif x == float(values[4]):
                out_value = weights[4]
            else:
                out_value = float(0)
            return out_value

        hydro_tif = xr.apply_ufunc(hydro_transform, ed_reclass, vectorize=True).astype(
            np.float32
        )
        hydro_tif = hydro_tif.copy(data=hydro_tif).rename("hydro_weight")

        hydro_tif = hydro_tif.rio.reproject(
            proj_crs, resolution=output_resolution, resampling=Resampling.bilinear
        )
        rasters.write(hydro_tif, os.path.join(tmp_dir, "hydro_weight.tif"))
        return hydro_tif
    else:
        return rasters.read(os.path.join(tmp_dir, "hydro_weight.tif"))


def aspect_raster(dem_b, aoi, output_path):
    """ """
    LOGGER.info("-- Produce the Aspect raster for the LSI model")
    if not os.path.exists(os.path.join(output_path, "aspect_weight.tif")):
        aspect_tif = aspect(dem_b)

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
        #    the following codes are not considered = {1: Flat, 2: North, 3: Northeast, 4: East}
        #    as they already are within the class defined in Aspect_dbf

        # Transform to Aspect_dbf scale
        aspect_reclass = xr.where(aspect_class == 5, 4, aspect_class)  # Southeast
        aspect_reclass = xr.where(aspect_reclass == 6, 5, aspect_reclass)  # South
        aspect_reclass = xr.where(aspect_reclass == 7, 5, aspect_reclass)  # Southwest
        aspect_reclass = xr.where(aspect_reclass == 8, 3, aspect_reclass)  # West
        aspect_reclass = xr.where(aspect_reclass == 9, 3, aspect_reclass)  # Northeast
        aspect_reclass = xr.where(aspect_reclass == 10, 2, aspect_reclass)  # North

        aspect_reclass_xr = (
            aspect_tif.copy(data=aspect_reclass)
            .astype(np.float32)
            .rename("aspect_weight")
        )

        def aspect_transform(x):
            values = [1, 2, 3, 4, 5]
            weights = [0.03118, 0.046272, 0.075935, 0.130972, 0.2164]
            if x == float(values[0]):
                out_value = weights[0]
            elif x == float(values[1]):
                out_value = weights[1]
            elif x == float(values[2]):
                out_value = weights[2]
            elif x == float(values[3]):
                out_value = weights[3]
            elif x == float(values[4]):
                out_value = weights[4]
            else:
                out_value = float(0)
            return out_value

        aspect_tif = xr.apply_ufunc(
            aspect_transform, aspect_reclass_xr, vectorize=True
        ).astype(np.float32)
        aspect_tif = aspect_tif.copy(data=aspect_tif).rename("aspect_weight")
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
    landcover_weight_dbf = gpd.read_file(landcover_weight_path)
    aoi_b = geometry.buffer(aoi_zone, REGULAR_BUFFER)
    aoi_m = geometry.buffer(aoi_zone, SMALL_BUFFER)
    try:
        landcover = rasters.crop(landcover_path, aoi_b)
    except ValueError:
        raise ValueError("Your AOI doesn't cover your Landcover layer.") from ValueError

    # Reclassification based on ELSUS
    landcover = rasters.collocate(
        reference_raster, landcover.astype(np.float32), Resampling.nearest
    )
    landcover_reclass = reclass_landcover_elsus(landcover, proj_crs, landcover_name)
    landcover_reclass = rasters.crop(landcover_reclass, aoi_m)

    # JOIN with LULC_dbf
    def landuse_transform(x):
        values = list(landcover_weight_dbf["Value"])
        weights = list(landcover_weight_dbf["Weight"])
        if x == float(values[0]):
            out_value = weights[0]
        elif x == float(values[1]):
            out_value = weights[1]
        elif x == float(values[2]):
            out_value = weights[2]
        elif x == float(values[3]):
            out_value = weights[3]
        elif x == float(values[4]):
            out_value = weights[4]
        elif x == float(values[5]):
            out_value = weights[5]
        elif x == float(values[6]):
            out_value = weights[6]
        else:
            out_value = float(0)
        return out_value

    landcover_tif = xr.apply_ufunc(
        landuse_transform, landcover_reclass, vectorize=True
    ).astype(np.float32)
    landcover_tif = landcover_tif.copy(data=landcover_tif).rename("landcover_weight")

    landcover_tif = rasters.crop(landcover_tif, aoi_m)
    # Calculating final Weights
    zone_class = "Z" + str(zone_id)  # Class for the zone
    final_weight_factor = final_weight_dbf[final_weight_dbf.Factor == "Landcover"][
        zone_class
    ].iloc[0]
    landcover_tif = landcover_tif * final_weight_factor

    landcover_tif = xr.where(
        landcover_tif > 10, np.nan, landcover_tif
    )  # There should not be values greater than 10 (border effect)
    landcover_tif = xr.where(
        landcover_tif < 0, np.nan, landcover_tif
    )  # There should not be negative values (border effect)
    landcover_tif = landcover_tif.rio.write_crs(proj_crs, inplace=True)

    # Write in memory
    temp_dir = os.path.join(
        tmp_dir, AnyPath("landcover_" + str(zone_id) + "_" + str(counter) + ".tif")
    )

    rasters.write(
        landcover_tif,
        temp_dir,
        predictor=1,
    )

    return temp_dir


def slope_raster_eu(
    slope_weight_path,
    dem,
    aoi_zone,
    proj_crs,
    zone_id,
    counter,
    final_weight_dbf,
    tmp_dir,
):
    # -- Slope raster computation
    slope_dbf = gpd.read_file(slope_weight_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        slope_degrees = rasters.slope(dem.astype(np.float32), in_rad=False)
    aoi_b = geometry.buffer(aoi_zone, REGULAR_BUFFER)
    aoi_m = geometry.buffer(aoi_zone, SMALL_BUFFER)
    # -- Slope Reclassification:
    # Define steps of classification depending on the zone
    if int(zone_id) == 0:
        SLOPE_STEPS = [0, 1, 9, 13, 21, 27, 35, 42, 90]
    elif int(zone_id) in range(1, 5):
        SLOPE_STEPS = [0, 1, 5, 9, 13, 17, 21, 31, 90]
    elif int(zone_id) in range(5, 7):
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

    # -- JOIN Slope with Weights

    def slope_transform(x):
        values = list(slope_dbf["Value"])
        weights = list(slope_dbf["Weight"])
        if x == float(values[0]):
            out_value = weights[0]
        elif x == float(values[1]):
            out_value = weights[1]
        elif x == float(values[2]):
            out_value = weights[2]
        elif x == float(values[3]):
            out_value = weights[3]
        elif x == float(values[4]):
            out_value = weights[4]
        elif x == float(values[5]):
            out_value = weights[5]
        elif x == float(values[6]):
            out_value = weights[6]
        elif x == float(values[7]):
            out_value = weights[7]
        else:
            out_value = float(0)
        return out_value

    slope_tif = xr.apply_ufunc(slope_transform, slope_d, vectorize=True).astype(
        np.float32
    )
    slope_tif = slope_tif.copy(data=slope_tif).rename("slope_weight")
    slope_tif = rasters.crop(slope_tif, aoi_m)

    # -- Apply Final Weights
    zone_class = "Z" + str(zone_id)  # Class for the zone
    final_weight_factor = final_weight_dbf[final_weight_dbf.Factor == "Slope"][
        zone_class
    ].iloc[0]

    # Apply factor
    slope_tif = slope_tif * final_weight_factor

    slope_tif = xr.where(
        slope_tif > 10, np.nan, slope_tif
    )  # There should not be values greater than 10
    slope_tif = xr.where(
        slope_tif < 0, np.nan, slope_tif
    )  # There should not be negative values
    slope_tif = slope_tif.rio.write_crs(proj_crs, inplace=True)

    # Write in memory
    temp_dir = os.path.join(
        tmp_dir, AnyPath("slope_" + str(zone_id) + "_" + str(counter) + ".tif")
    )

    rasters.write(
        slope_tif,
        temp_dir,
        predictor=1,
    )

    return temp_dir
