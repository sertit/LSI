"""

This file is part of LSI.

LSI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

LSI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with LSI. If not, see <https://www.gnu.org/licenses/>.

"""

import logging.handlers
import sys

import ftep_util as ftep
from sertit import s3

DEBUG = False
LOGGING_FORMAT = "%(asctime)s - [%(levelname)s] - %(message)s"
LOGGER = logging.getLogger("OSM Charter")

FTEP_S3_ENDPOINT = "s3.waw2-1.cloudferro.com"


def ftep_s3_env(*args, **kwargs):
    return s3.s3_env(endpoint=FTEP_S3_ENDPOINT)(*args, **kwargs)


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
        InputParameters.LOCATION.value: "Global",  # params.getString("location"),
        InputParameters.DEM_NAME.value: params.getString("dem_name"),
        InputParameters.OTHER_DEM_PATH.value: None,
        InputParameters.LANDCOVER_NAME.value: "ESA WorldCover - 2021 (10m)",  # params.getString("lulc"),
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
