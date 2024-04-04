""" lsi : main script with CLI """
import argparse
import logging
import sys

import click

# from lsi import LOGGER_NAME # on the meantime to solve the acces to lsi.py
# from lsi.lsi_core import lsi_core # on the meantime to solve the acces to lsi.py
from lsi_core import LOGGER, LOGGING_FORMAT, DataPath, InputParameters, lsi_core
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
    required=True,
)
@click.option(
    "-loc",
    "--location",
    help="Location of the AOI",
    type=click.Choice(["Europe", "Global"]),  # to be added maybe Europe/Global_Legacy
    required=True,
    show_default=True,
)
@click.option(
    "-litho",
    "--lithology_path",
    help="Geodatabases of lithologies.",
    type=click.Path(exists=True, resolve_path=True),
)
@click.option(
    "-lulc",
    "--landcover_name",
    help="Land Cover Name",
    type=click.Choice(
        [
            "ESA Worldcover - 2021 (10m)"
            # "Corine Land Cover - 2018 (100m)",
            # "Global Land Cover - Copernicus 2019 (100m)",
            # "P03",
        ]
    ),
    default="ESA Worldcover - 2021 (10m)",
    show_default=True,
)
@click.option(
    "-final_weights",
    "--final_weights_path",
    help="Geotadabase with the final weights for the LSI computation.",
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
    "-o",
    "--output",
    help="Output directory. ",
    type=click.Path(file_okay=False, resolve_path=True, writable=True),
    required=True,
)
@click.option(
    "--ftep",
    help="Set this flag if the command line is run on the ftep platform. ",
    default=False,
)
def main(
    aoi_path,
    location,
    lithology_path,
    landcover_name,
    final_weights_path,
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
            InputParameters.LITHOLOGY_PATH.value: lithology_path,
            InputParameters.LANDCOVER_NAME.value: landcover_name,
            InputParameters.FINAL_WEIGHTS_PATH.value: final_weights_path,
            InputParameters.REF_EPSG.value: epsg_code,
            InputParameters.OUTPUT_DIR.value: output_path,
        }
        DataPath.load_paths(ftep)

    try:
        lsi()
        LOGGER.info("lsi was a success.")
        sys.exit(0)

        print("success")

    # pylint: disable=W0703
    except Exception:
        LOGGER.error("lsi has failed:", exc_info=True)
        sys.exit(1)
        print("not sucess")


# def lsi():
#    """rrm_csv_coordinates with the command line"""
#    parser = argparse.ArgumentParser()

# parser.add_argument(
#     "-a",
#     "--aoi",
#     help="AOI path as shapefile.",
#     type=to_abspath,
#     required=True,
# )

#     parser.add_argument(
#     "-loc",
#     "--location",
#     help="Location",
#     choices=["Global"],
#     type=str,
#     required=True,
# ),

# parser.add_argument(
#     "-dem", "--dem", help="Please specify a DEM (30m).", required=True, type=to_abspath,
# )

# parser.add_argument(
#     "-odem", "--other_dem", help="Provide any other DEM", required=False, type=to_abspath,
# )

# parser.add_argument(
#     "-lc", "--landcover", help="Please specify a Landcover.", required=True, type=to_abspath,
# )

# parser.add_argument(
#     "-litho", "--litho_db", help="Lithology and permeability geodatabase.", required=True, type=to_abspath,
# )

# parser.add_argument(
#     "-field", "--field_geology", help="=Field geology to be used.", required=True
#     , type=str, choices=["Rating"]
# )

# parser.add_argument(
#     "-geom", "--geomor_db", help="=Database with Aspect, Elevation, Final weights, Geology, Hydrology, Land use and Slope data.", required=True
# )

# parser.add_argument(
#     "-epsg", "--epsg_code", help="EPSG code", type=int, required=True
# )

# parser.add_argument(
#     "-out", "--output", help="=Output directory.", required=True, type=to_abspath,
# )

# Parse args
#   args = parser.parse_args()  # noqa


# Process
# (Here to place the main function with all the LSI process from lsi_core.py)
#   lsi_core()


if __name__ == "__main__":
    main()
