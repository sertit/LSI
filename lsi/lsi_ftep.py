# -*- coding: utf-8 -*-
# This file is part of LSI.
# LSI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# LSI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with LSI. If not, see <https://www.gnu.org/licenses/>.

import logging
import logging.handlers
import os
import sys
from functools import wraps

import ftep_util as ftep
from sertit.logs import SU_NAME

# from sertit import s3
from sertit.s3 import define_s3_client

LOGGER = logging.getLogger(SU_NAME)

DEBUG = False
LOGGING_FORMAT = "%(asctime)s - [%(levelname)s] - %(message)s"
LOGGER = logging.getLogger("OSM Charter")

FTEP_S3_ENDPOINT = "s3.waw2-1.cloudferro.com"


AWS_ACCESS_KEY_ID = "AWS_ACCESS_KEY_ID"
"""
Environment variable linked to AWS Access Key ID.
"""

AWS_SECRET_ACCESS_KEY = "AWS_SECRET_ACCESS_KEY"
"""
Environment variable linked to AWS Secret Access Key.
"""

AWS_S3_ENDPOINT = "AWS_S3_ENDPOINT"
"""
Environment variable linked to AWS endpoint.
"""

USE_S3_STORAGE = "USE_S3_STORAGE"
"""
Environment variable created to use Unistra's S3 bucket.
"""


def s3_env(*args, **kwargs):
    """
    Create S3 compatible storage environment
    You need to set endpoint url if you use s3 compatible storage
    since GDAL/Rasterio does not read endpoint url from config file.

    This function searches for S3 configuration in many places.
    It does apply configuration variables precedence, and you might have a use for it.
    Here is the order of precedence from least to greatest
    (the last listed configuration variables override all other variables):
        1. AWS profile
        2. Given endpoint_url as function argument
        3. AWS environment variable

    Returns:
        Callable: decorated function

    Example:
        >>> from sertit.s3 import s3_env
        >>> from sertit import AnyPath
        >>> @s3_env(endpoint="s3.unistra.fr")
        >>> def file_exists(path: str):
        >>>     pth = AnyPath(path)
        >>>     print(pth.exists())
        >>> file_exists("s3://sertit-geodatastore/GLOBAL/COPDEM_30m/COPDEM_30m.vrt")
        True
    """
    import rasterio

    use_s3 = kwargs.get("use_s3_env_var", USE_S3_STORAGE)
    requester_pays = kwargs.get("requester_pays")
    no_sign_request = kwargs.get("no_sign_request")
    endpoint = os.getenv(AWS_S3_ENDPOINT, kwargs.get("endpoint"))
    profile_name = kwargs.get("profile_name", None)

    def decorator(function):
        @wraps(function)
        def s3_env_wrapper(*_args, **_kwargs):
            """S3 environment wrapper"""
            if int(os.getenv(use_s3, 1)):
                args_rasterio = {
                    "profile_name": profile_name,
                    "CPL_CURL_VERBOSE": False,
                    "AWS_VIRTUAL_HOSTING": False,
                    "GDAL_DISABLE_READDIR_ON_OPEN": False,
                    "AWS_NO_SIGN_REQUEST": "YES" if no_sign_request else "NO",
                    "AWS_REQUEST_PAYER": "requester" if requester_pays else None,
                }
                args_s3_client = {
                    "profile_name": profile_name,
                    "requester_pays": requester_pays,
                    "no_sign_request": no_sign_request,
                }
                args_s3_client.update(kwargs)

                if endpoint is not None:
                    args_rasterio["AWS_S3_ENDPOINT"] = endpoint
                    args_s3_client["endpoint_url"] = (
                        f"https://{endpoint}"  # cloudpathlib can read endpoint from config file
                    )

                # Define S3 client for S3 paths
                define_s3_client(**args_s3_client)
                os.environ[use_s3] = "1"
                LOGGER.info("Using S3 files")
                with rasterio.Env(**args_rasterio):
                    return function(*_args, **_kwargs)

            else:
                os.environ[use_s3] = "0"
                LOGGER.info("Using on disk files")
                return function(*_args, **_kwargs)

        return s3_env_wrapper

    return decorator


def ftep_s3_env(*args, **kwargs):
    return s3_env(endpoint=FTEP_S3_ENDPOINT)(*args, **kwargs)


@ftep_s3_env
def compute_lsi():
    parameters_file_path = "/home/worker/workDir/FTEP-WPS-INPUT.properties"
    # Default parameter values
    params = ftep.Params(parameters_file_path)

    # Add parameters from file
    params.readFile(parameters_file_path)

    # --- Parameters ---
    # Load inputs

    from sertit import logs

    from lsi.lsi_core import LOGGER, LOGGING_FORMAT, DataPath, InputParameters, lsi_core

    logs.init_logger(LOGGER, logging.INFO, LOGGING_FORMAT)
    LOGGER.info("--- LSI ---")

    input_dict = {
        InputParameters.AOI_PATH.value: params.getString("aoi"),
        InputParameters.LOCATION.value: params.getString("location"),
        InputParameters.DEM_NAME.value: params.getString("dem_name"),
        InputParameters.OTHER_DEM_PATH.value: None,
        InputParameters.LANDCOVER_NAME.value: params.getString("lulc"),
        InputParameters.EUROPE_METHOD.value: "Refined",
        InputParameters.OUTPUT_RESOLUTION.value: 10,
        InputParameters.REF_EPSG.value: None,
        InputParameters.OUTPUT_DIR.value: "/home/worker/workDir/outDir/output",
        InputParameters.TEMP.value: False,
        InputParameters.JENKS.value: True,
    }
    DataPath.load_paths(ftep=True)

    try:
        # Compute LSI charter
        lsi_core(input_dict)
        LOGGER.info("--- LSI was a success.")
        sys.exit(0)

    except Exception as ex:  # noqa
        LOGGER.error("LSI has failed:", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    compute_lsi()
