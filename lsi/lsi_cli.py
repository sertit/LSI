""" lsi : main script with CLI """
import argparse
import logging
import sys
try:
    import rich_click as click
except:
    import click


# from lsi import LOGGER_NAME # on the meantime to solve the acces to lsi.py
# from lsi.lsi_core import lsi_core # on the meantime to solve the acces to lsi.py
from lsi.lsi_core import LOGGER, LOGGING_FORMAT, DataPath, InputParameters, lsi_core
from sertit import logs
from sertit.files import to_abspath
from sertit.logs import LOGGING_FORMAT
from sertit.unistra import unistra_s3

# from lsi_core import LOGGER # on the meantime to solve the acces to lsi.py
# LOGGER = logging.getLogger(LOGGER_NAME) # on the meantime to solve the acces to lsi.py


@click.command()
@click.option(
    "-aoi",
    "--aoi_path",
    help="Path to the AOI (shp, geojson) or WKT string",
    type=click.Path(exists=True, resolve_path=True),
    required=True,
)
@click.option(
    "-loc",
    "--location",
    help="Location of the AOI",
    type=click.Choice(["Europe", "Global"]),  # to be added maybe Europe/Global_Legacy
    required=True,
)
@click.option(
    "-dem",
    "--dem_name",
    help="DEM Name needed",
    type=click.Choice(["COPDEM 30m", "FABDEM", "SRTM 30m", "Other"]), #"EUDEM 25m" ,  "MERIT 5 deg"
    default="COPDEM 30m",
    show_default=True,
)
@click.option(
    "-demp",
    "--other_dem_path",
    help="DEM path if dem = Other",
    type=click.Path(exists=True, resolve_path=True),
)
@click.option(
    "-litho",
    "--lithology_path",
    help="GDB of lithologies.",
    type=click.Path(exists=True, resolve_path=True),
)
@click.option(
    "-lulc",
    "--landcover_name",
    help="Land Cover Name",
    type=click.Choice(
        [
            "ESA WorldCover - 2021 (10m)"
            # "Corine Land Cover - 2018 (100m)",
            # "Global Land Cover - Copernicus 2019 (100m)",
            # "P03",
        ]
    ),
    default="ESA Worldcover - 2021 (10m)",
    show_default=True,
)
@click.option(
    "-weights",
    "--weights_path",
    help="Geotadabase with the weights for the LSI computation.",
    type=click.Path(exists=True, resolve_path=True),
)
@click.option(
    "-epsg",
    "--epsg_code",
    help="EPSG code, 4326 is not accepted. By default, it is the EPSG code of the AOI UTM zone.",
    type=click.IntRange(min=1024, max=32767),
    show_default=True,
)
@click.option(
    "-out",
    "--output_path",
    help="Output directory. ",
    type=click.Path(file_okay=False, resolve_path=True, writable=True),
    required=True,
)
@click.option(
    "--ftep",
    help="Set this flag if the command line is run on the ftep platform. ",
    default=False,
)


def compute_lsi(
    aoi_path,
    location,
    dem_name,
    other_dem_path,
    lithology_path,
    landcover_name,
    weights_path,
    epsg_code,
    output_path,
    ftep,
):
    logs.init_logger(LOGGER, logging.INFO, LOGGING_FORMAT)
    LOGGER.info("--- LSI ---")

    with unistra_s3():
        # Insert args in a dict
        input_dict = {
            InputParameters.AOI_PATH.value: aoi_path,
            InputParameters.LOCATION.value: location,
            InputParameters.DEM_NAME.value: dem_name,
            InputParameters.OTHER_DEM_PATH.value: other_dem_path,
            InputParameters.LITHOLOGY_PATH.value: lithology_path,
            InputParameters.LANDCOVER_NAME.value: landcover_name,
            InputParameters.WEIGHTS_PATH.value: weights_path,
            InputParameters.REF_EPSG.value: epsg_code,
            InputParameters.OUTPUT_DIR.value: output_path,
        }
        DataPath.load_paths(ftep)

    try:
        lsi_core(input_dict)
        LOGGER.info("lsi was a success.")
        sys.exit(0)

        print("success")

    # pylint: disable=W0703
    except Exception as ex:
        LOGGER.error("lsi has failed:", exc_info=True)
        print("not sucess")
        sys.exit(1)


if __name__ == "__main__":
    compute_lsi()
