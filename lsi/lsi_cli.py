""" lsi : main script with CLI """

import logging
import sys

import rich_click as click
from sertit import logs
from sertit.logs import LOGGING_FORMAT
from sertit.unistra import unistra_s3

from lsi.lsi_core import LOGGER, DataPath, InputParameters, lsi_core


@click.command(
    context_settings=dict(
        help_option_names=["-h", "--help"],
        max_content_width=300,
        show_default=True,
    )
)
@click.option(
    "-aoi",
    "--aoi",
    help="AOI (shp, geojson) or WKT string",
    type=click.Path(exists=True, resolve_path=True),
    required=True,
)
@click.option(
    "-loc",
    "--location",
    help="Location of the AOI",
    type=click.Choice(["Europe", "Global"]),  # to be added maybe Europe/Global_Legacy
    default="Global",
)
@click.option(
    "-dem",
    "--dem_name",
    help="DEM Name needed",
    type=click.Choice(
        ["COPDEM 30m", "FABDEM", "SRTM 30m", "Other"]
    ),  # "EUDEM 25m" ,  "MERIT 5 deg"
    default="COPDEM 30m",
    show_default=True,
)
@click.option(
    "-demp",
    "--other_dem",
    help="DEM path if dem = Other",
    type=click.Path(exists=True, resolve_path=True),
)
@click.option(
    "-lulc",
    "--landcover_name",
    help="Land Cover Name",
    type=click.Choice(
        [
            "ESA WorldCover - 2021 (10m)",
            "Corine Land Cover - 2018 (100m)",
            # "Global Land Cover - Copernicus 2019 (100m)",
            # "P03",
        ]
    ),
    default="ESA WorldCover - 2021 (10m)",
    show_default=True,
)
@click.option(
    "-eu_method",
    "--europe_method",
    help="if LOCATION = EUROPE, choose whether you want a fast computation with lower resolution (based on the pre-existent ELSUS layer) or a refined LSI computation",
    type=click.Choice(
        [
            "Refined",
            "Fast",  # ELSUS layer
        ]
    ),
    default="Refined",
)
@click.option(
    "-res",
    "--output_resolution",
    help="Output resolution. Taking from DEM if not provided",
    type=click.IntRange(min=1, max=1000),
    #    default=30,
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
    help="Output directory.",
    type=click.Path(file_okay=False, resolve_path=True, writable=True),
    required=True,
)
@click.option(
    "--ftep",
    help="Set this flag if the command line is run on the ftep platform. ",
    default=False,
)
def compute_lsi(
    aoi,
    location,
    dem_name,
    other_dem,
    landcover_name,
    europe_method,
    output_resolution,
    epsg_code,
    output_path,
    ftep,
):
    logs.init_logger(LOGGER, logging.INFO, LOGGING_FORMAT)
    LOGGER.info("--- LSI ---")

    with unistra_s3():
        # Insert args in a dict
        input_dict = {
            InputParameters.AOI_PATH.value: aoi,
            InputParameters.LOCATION.value: location,
            InputParameters.DEM_NAME.value: dem_name,
            InputParameters.OTHER_DEM_PATH.value: other_dem,
            InputParameters.LANDCOVER_NAME.value: landcover_name,
            InputParameters.EUROPE_METHOD.value: europe_method,
            InputParameters.OUTPUT_RESOLUTION.value: output_resolution,
            InputParameters.REF_EPSG.value: epsg_code,
            InputParameters.OUTPUT_DIR.value: output_path,
        }
        DataPath.load_paths(ftep)

        try:
            lsi_core(input_dict)
            LOGGER.info("lsi was a success.")
            sys.exit(0)

        # pylint: disable=W0703
        except Exception:  # noqa
            LOGGER.error("lsi has failed:", exc_info=True)
            sys.exit(1)


if __name__ == "__main__":
    compute_lsi()
