"""
This file is part of LSI.
LSI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
LSI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with LSI. If not, see <https://www.gnu.org/licenses/>.
"""

"""Tests"""

import logging
import os
import tempfile

# import pytest
from sertit import AnyPath, ci  # noqa
from sertit.types import AnyPathType
from sertit.unistra import s3_env

from lsi.lsi_core import DataPath, InputParameters, lsi_core

logging.getLogger("sertit").setLevel(logging.INFO)
ci.reduce_verbosity()


# @pytest.fixture(scope="session", autouse=True)
# def set_env():
# os.environ["NUMBA_DEBUG"] = "0"


def get_ci_path() -> AnyPathType:
    """
    Get the path to the CI folder.
    """
    s3_path = AnyPath("s3://sertit-ci")

    return s3_path / "lsi" / "Test_Austria"


def get_directories() -> AnyPathType:
    """
    Get the path to the CI folder and its subdirectories.
    """
    ci_directory = get_ci_path()
    aoi_directory = ci_directory / "AOI"
    output_directory = ci_directory / "out_expected"

    return ci_directory, aoi_directory, output_directory


@s3_env
def test_lsi_global():
    """ """
    # Make directory paths
    _, aoi_directory, output_directory = get_directories()

    aoi_path = aoi_directory / "ci_aoi_salzburg.shp"
    expected_output_tif = output_directory / "GLOBAL" / "LandslideSusceptibility.tif"
    expected_output_shp = output_directory / "GLOBAL" / "FER_LR_ave.shp"

    assert expected_output_tif.exists()
    assert expected_output_shp.exists()

    with tempfile.TemporaryDirectory() as output:
        input_dict = {
            InputParameters.AOI_PATH.value: str(aoi_path),
            InputParameters.LOCATION.value: "Global",
            InputParameters.DEM_NAME.value: "COPDEM 30m",
            InputParameters.LANDCOVER_NAME.value: "Global Land Cover - Copernicus 2019 (100m)",
            InputParameters.EUROPE_METHOD.value: "Refined",  # Not relevant for this test
            InputParameters.OUTPUT_RESOLUTION.value: 50,
            InputParameters.TEMP.value: False,  # Don't keep temporary files
            InputParameters.JENKS.value: False,  # Apply jenks breaks to compute LandslideRisk
            InputParameters.OUTPUT_DIR.value: str(output),
        }
        DataPath.load_paths()
        lsi_core(input_dict=input_dict, ftep=False)

        output_classification_tif = os.path.join(output, "LandslideSusceptibility.tif")
        output_classification_shp = os.path.join(output, "FER_LR_ave.shp")
        ci.assert_geom_almost_equal(expected_output_shp, output_classification_shp)
        ci.assert_raster_max_mismatch(
            expected_output_tif, output_classification_tif, max_mismatch_pct=3
        )


@s3_env
def test_lsi_europe():
    """ """
    # Make directory paths
    _, aoi_directory, output_directory = get_directories()

    aoi_path = aoi_directory / "ci_aoi_salzburg.shp"
    expected_output_tif = output_directory / "EUROPE" / "LandslideSusceptibility.tif"
    expected_output_shp = output_directory / "EUROPE" / "FER_LR_ave.shp"

    assert expected_output_tif.exists()
    assert expected_output_shp.exists()

    with tempfile.TemporaryDirectory() as output:
        input_dict = {
            InputParameters.AOI_PATH.value: str(aoi_path),
            InputParameters.LOCATION.value: "Europe",
            InputParameters.DEM_NAME.value: "COPDEM 30m",
            InputParameters.LANDCOVER_NAME.value: "Corine Land Cover - 2018 (100m)",
            InputParameters.EUROPE_METHOD.value: "Refined",  # Not relevant for this test
            InputParameters.OUTPUT_RESOLUTION.value: 50,
            InputParameters.TEMP.value: False,  # Don't keep temporary files
            InputParameters.JENKS.value: False,  # Apply jenks breaks to compute LandslideRisk
            InputParameters.OUTPUT_DIR.value: str(output),
        }
        DataPath.load_paths()
        lsi_core(input_dict=input_dict, ftep=False)

        output_classification_tif = os.path.join(output, "LandslideSusceptibility.tif")
        output_classification_shp = os.path.join(output, "FER_LR_ave.shp")
        ci.assert_geom_almost_equal(expected_output_shp, output_classification_shp)
        ci.assert_raster_max_mismatch(
            expected_output_tif, output_classification_tif, max_mismatch_pct=3
        )
