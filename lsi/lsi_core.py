""" lsi """
import logging
from enum import unique

from sertit.misc import ListEnum
from sertit import AnyPath
from sertit.unistra import get_geodatastore, s3_env

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
        cls.WORLD_COUNTRIES_PATH = (
            cls.GLOBAL_DIR / "World_countries_poly" / "world_countries_poly.shp"
        )


@unique
class InputParameters(ListEnum):
    """
    List of the input parameters
    """

    AOI_PATH = "aoi_path"
    LOCATION = "location"
    LITHOLOGY_PATH = "lithology_path"
    LANDCOVER_NAME = "landcover_name"
    FINAL_WEIGHTS_PATH = "final_weights_path"
    REF_EPSG = "epsg_code"
    OUTPUT_DIR = "output_path"


def lsi_core() -> str:
    """
    TODO: Complete arguments and dosctring
    """

    # 1. GEOLOGY

    # 2. SLOPE (DEGREES)

    # 3. LANDCOVER

    # 4. ELEVATION

    # 5. DISTANCE TO HYDRO

    # 6. ASPECT

    # 7. FINAL WEIGTHS

    # 8. LSI RASTER CALCULATION






    raise NotImplementedError
