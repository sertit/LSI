import arcpy


class Toolbox:
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Sertit"
        self.alias = "Sertit"

        # List of tool classes associated with this toolbox
        self.tools = [Lsi]


class Lsi:
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lsi"
        self.description = ""
        self.canRunInBackground = False
        self.category = "RRM"

    def getParameterInfo(self):
        """Define parameter definitions"""
        # Define parameter definitions

        # 0. AOI
        aoi = arcpy.Parameter(
            displayName="Aoi",
            name="aoi",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input",
        )
        # 1. Location
        location = arcpy.Parameter(
            displayName="Location",
            name="location",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
        )
        location.filter.type = "ValueList"
        location.filter.list = ["Europe", "Global"]
        location.value = "Global"
        # 2. Method for Location=Europe
        method = arcpy.Parameter(
            displayName="Method",
            name="method",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
        )
        method.filter.type = "ValueList"
        method.filter.list = ["Refined", "Fast"]
        method.value = "Refined"
        # 3. Landcover
        landcover = arcpy.Parameter(
            displayName="Landcover",
            name="landcover",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
            category="Advanced",
        )

        landcover.filter.type = "ValueList"
        landcover.filter.list = [
            "ESA WorldCover - 2021 (10m)",
            "Corine Land Cover - 2018 (100m)",
            "Global Land Cover - Copernicus 2019 (100m)",
        ]
        landcover.value = "ESA WorldCover - 2021 (10m)"
        # 4. DEM
        dem = arcpy.Parameter(
            displayName="DEM",
            name="dem",
            datatype="GPString",
            parameterType="Optional",
            direction="Input",
            category="Advanced",
        )

        dem.filter.type = "ValueList"
        dem.filter.list = ["COPDEM 30m", "FABDEM", "Other"]
        dem.value = "COPDEM 30m"
        # 5. Dem Raster path
        dem_raster_path = arcpy.Parameter(
            displayName="Dem raster",
            name="dem_raster",
            datatype="DEFile",
            parameterType="Optional",
            direction="Input",
            category="Advanced",
        )
        # 6. Output resolution
        output_resolution = arcpy.Parameter(
            displayName="Output resolution",
            name="output_resolution",
            datatype="GPDouble",
            direction="Input",
            category="Advanced",
        )
        output_resolution.value = 30
        # 7. Output folder
        output_folder = arcpy.Parameter(
            displayName="Output folder",
            name="output_folder",
            datatype="DEFolder",
            parameterType="Required",
            direction="Input",
        )

        # 8. Temp folder
        temp_folder = arcpy.Parameter(
            displayName="Keep temporal folder",
            name="temp_folder",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
            category="Advanced",
        )

        temp_folder.filter.type = "ValueList"
        temp_folder.filter.list = ["Yes", "No"]
        temp_folder.value = "Yes"

        # 9. Jenks
        jenks_class = arcpy.Parameter(
            displayName="Apply Jenks breaks classification",
            name="jenks_class",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
            category="Advanced",
        )

        jenks_class.filter.type = "ValueList"
        jenks_class.filter.list = ["Yes", "No"]
        jenks_class.value = "Yes"

        params = [
            aoi,
            location,
            method,
            landcover,
            dem,
            dem_raster_path,
            output_resolution,
            output_folder,
            temp_folder,
            jenks_class,
        ]

        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""

        if parameters[1].value == "Europe":  # Location
            parameters[2].enabled = True  # ELSUS method
        else:
            parameters[2].enabled = False

        if parameters[4].value == "Other":  # DEM
            parameters[5].enabled = True  # Other DEM path enabling
        else:
            parameters[5].enabled = False

        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        import pathlib
        import sys

        # Don't remove these lines
        tools_path = pathlib.Path(__file__).parent

        # The tool is run from sertit_atools so add sertit_atools to python path
        if tools_path.name == "sertit_atools":
            tools_path = str(tools_path.absolute())
        # The tool is run from this project so only add the root folder to python path
        else:
            tools_path = str(tools_path.parent.absolute())
        if tools_path not in sys.path:
            sys.path.append(tools_path)

        main_arcgis(parameters, messages)
        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return


def epsg_from_arcgis_proj(arcgis_proj):
    """
    Extract espg code from arcgis proj
    Args:
        arcgis_proj () : Arcgis proj

    Returns:
        epsg_code : ndarray of the reclassified a raster
    """
    try:
        sr = arcpy.SpatialReference()
        sr.loadFromString(arcgis_proj)
        epsg_code = sr.factoryCode

    except:
        raise ValueError(
            "Input coordinate system is not from Arcgis coordinate system tools"
        )

    return epsg_code


def main_arcgis(parameters, messages):
    """
    Main function of your arcgis tool .

    To be completed.
    """
    import logging.handlers

    from sertit.arcpy import (
        ArcPyLogger,
        ArcPyLogHandler,
        feature_layer_to_path,
        init_conda_arcpy_env,
    )

    init_conda_arcpy_env()

    from sertit.arcpy import ArcPyLogHandler

    from lsi import LOGGER_NAME
    from lsi.lsi_core import DataPath, InputParameters, lsi_core  # noqa

    arcpy_logger = ArcPyLogger("LSI")
    logger = logging.getLogger(LOGGER_NAME)
    # handler = ArcPyLogHandler(
    #     "output_log.log", maxBytes=1024 * 1024 * 2, backupCount=10  # 2MB log files
    # )
    # formatter = logging.Formatter("%(levelname)-8s %(message)s")
    # handler.setFormatter(formatter)
    # logger.addHandler(handler)
    logger = arcpy_logger.logger
    logger.setLevel(logging.DEBUG)

    # --- ENV VAR ---
    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension("Spatial")

    # -> Check extensions in used in this tool
    # arcpy.CheckOutExtension("Spatial")

    logger.info("lsi_core")

    # -> Load inputs
    # aoi = parameters[0].valueAsText
    # output = parameters[1].valueAsText

    aoi_path = feature_layer_to_path(parameters[0].value)

    if parameters[8].value == "Yes":  # Keep temporal folders?
        temp = True
    else:
        temp = False

    if parameters[9].value == "Yes":  # Jenks
        jenks = True
    else:
        jenks = False

    # --- Parameters ---
    # Load inputs

    input_dict = {
        InputParameters.AOI_PATH.value: aoi_path,
        InputParameters.LOCATION.value: parameters[1].valueAsText,
        InputParameters.DEM_NAME.value: parameters[4].valueAsText,
        InputParameters.OTHER_DEM_PATH.value: parameters[5].valueAsText,
        InputParameters.LANDCOVER_NAME.value: parameters[3].valueAsText,
        InputParameters.EUROPE_METHOD.value: parameters[2].valueAsText,
        InputParameters.OUTPUT_RESOLUTION.value: parameters[6].valueAsText,
        InputParameters.OUTPUT_DIR.value: parameters[7].valueAsText,
        InputParameters.TEMP.value: temp,
        InputParameters.JENKS.value: jenks,
    }
    DataPath.load_paths()
    try:
        lsi_core(input_dict)
        logger.info("lsi is a success.")

    except Exception:
        import traceback

        logger.error("lsi has failed: %s", traceback.format_exc())

    # finally:
    #     logger.removeHandler(handler)
