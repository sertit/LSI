"""
This file is part of LSI.
LSI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
LSI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with LSI. If not, see <https://www.gnu.org/licenses/>.
"""

"""Tests"""

import os
import tempfile

from sertit import AnyPath, ci  # noqa
from sertit.unistra import s3_env

from lsi.lsi_core import DataPath, lsi_core

ci.reduce_verbosity()

os.environ["USE_S3_STORAGE"] = "1"
os.environ["AWS_S3_ENDPOINT"] = "s3.unistra.fr"


@s3_env
def test_lsi_global():
    ci_path = AnyPath("s3://sertit-ci") / "lsi" / "Test_Austria"
    aoi_path = ci_path / "AOI" / "ci_aoi_salzburg.shp"
    expected_output_tif = (
        ci_path / "out_expected" / "GLOBAL" / "LandslideSusceptibility.tif"
    )
    expected_output_shp = ci_path / "out_expected" / "GLOBAL" / "FER_LR_ave.shp"

    assert expected_output_tif.exists()
    assert expected_output_shp.exists()

    with tempfile.TemporaryDirectory() as output:
        input_dict = {
            "aoi": str(aoi_path),
            "location": "Global",
            "dem_name": "COPDEM 30m",
            "landcover_name": "ESA WorldCover - 2021 (10m)",
            "europe_method": "Refined",  # Not relevant for this test
            "output_resolution": 30,
            "temp": False,  # Don't keep temporary files
            "jenks": True,  # Apply jenks breaks to compute LandslideRisk
            "output_path": output,
        }
        DataPath.load_paths()
        lsi_core(input_dict=input_dict)

        output_classification_tif = os.path.join(output, "LandslideSusceptibility.tif")
        output_classification_shp = os.path.join(output, "FER_LR_ave.shp")
        ci.assert_raster_max_mismatch(
            expected_output_tif, output_classification_tif, max_mismatch_pct=3
        )
        ci.assert_raster_max_mismatch(
            expected_output_shp, output_classification_shp, max_mismatch_pct=3
        )
