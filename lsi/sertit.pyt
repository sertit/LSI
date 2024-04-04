import arcpy


class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Sertit"
        self.alias = "Sertit"

        # List of tool classes associated with this toolbox
        self.tools = [Template]


class Template(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Template"
        self.description = ""
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        # Define parameter definitions

        # AOI
        aoi = arcpy.Parameter(
            displayName="Aoi",
            name="aoi",
            datatype="DEFile",
            parameterType="Optional",
            direction="Input",
        )

        # Third parameter
        output_folder = arcpy.Parameter(
            displayName="Output Folder",
            name="output_folder",
            datatype="DEFolder",
            parameterType="Optional",
            direction="Input",
        )

        params = [aoi, output_folder]

        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
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
        tools_path = pathlib.Path(__file__).parent.parent

        # The tool is run from sertit_atools so add sertit_atools to python path
        if tools_path.parent.name == "sertit_atools":
            tools_path = str(tools_path.parent.absolute())
        # The tool is run from this project so only add the root folder to python path
        else:
            tools_path = str(tools_path.absolute())
        if tools_path not in sys.path:
            sys.path.append(tools_path)

        main_arcgis(parameters, messages)
        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return


def main_arcgis(parameters, messages):
    """
    Main function of your arcgis tool .

    To be completed.
    """
    import logging.handlers

    from sertit.arcpy import ArcPyLogHandler
    from xxx_to_be_renamed import LOGGER_NAME
    from xxx_to_be_renamed.xxx_to_be_renamed_core import xxx_to_be_renamed_core  # noqa

    logger = logging.getLogger(LOGGER_NAME)
    handler = ArcPyLogHandler(
        "output_log.log", maxBytes=1024 * 1024 * 2, backupCount=10  # 2MB log files
    )
    formatter = logging.Formatter("%(levelname)-8s %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)

    # -> Check extensions in used in this tool
    # arcpy.CheckOutExtension("Spatial")

    logger.info("xxx_to_be_renamed_core")

    # -> Load inputs
    # aoi = parameters[0].valueAsText
    # output = parameters[1].valueAsText

    try:
        xxx_to_be_renamed_core()
        logger.info("xxx_to_be_renamed is a success.")

    except Exception:
        import traceback

        logger.error("xxx_to_be_renamed has failed: %s", traceback.format_exc())
    finally:
        logger.removeHandler(handler)
