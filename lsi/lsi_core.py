""" lsi """
import logging
from enum import unique

from sertit.misc import ListEnum

# from lsi import LOGGER_NAME # on the meantime to solve the acces to lsi.py

DEBUG = False
LOGGING_FORMAT = "%(asctime)s - [%(levelname)s] - %(message)s"
LOGGER = logging.getLogger("LSI")  # on the meantime to solve the acces to lsi.py


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
    raise NotImplementedError
