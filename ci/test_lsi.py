# -*- coding: utf-8 -*-
# This file is part of LSI.
# LSI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# LSI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with LSI. If not, see <https://www.gnu.org/licenses/>.
""" Tests """

import os
import tempfile

from sertit import ci  # noqa


def get_ci_path():
    """Get ci DATA path"""
    return os.path.dirname(os.path.realpath(__file__))


def test_xxx():
    """Test WGS84 + core"""
    xxx_path = os.path.join(get_ci_path(), "xxx")  # noqa
    with tempfile.TemporaryDirectory() as tmp:  # noqa
        pass
