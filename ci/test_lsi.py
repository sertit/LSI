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
