# LSI

> :warning: **CONDA ENVIRONMENT**  
> **Make sure your ArcGIS Pro conda environment is properly installed.**  
> For more information: [ArcGIS Pro EO](https://lab.egeos-services.it/bitbucket/projects/CPP/repos/arcgis-pro-eo/browse).

## Description

Computes de Landslide Susceptibility Index (LSI) based on the P17-Post-disaster landslide methodology. 
This Index is estimated by a combination of data sources that includes the assesment of terrain, geomorphology and lithotypes, such as the following:

* Geology/Lithology
* Terrain Slope
* Slope Aspect
* Elevation
* Land use and land cover (LULC)
* Distance from drainage network


# ArcGIS Pro inputs:

(Layers)

TODO

## Run from the command line

```text
Usage: lsi.py [OPTIONS]

   TODO
  [Description ]

Options:
  -aoi, --aoi_path FILE                  Area of Interest for the calculation
                                         of the LSI  [required]
  -loc, --location STRING                Location of the AOI studied (Europe 
                                         or Global).  [required]
  -dem, --dem_name STRING                Name of the DEM to be used.
  -demp, --other_dem_path FILE           External DEM to be used in the calculations
                                         for the LSI. Provided by the user.
  -litho, --lithology_path FILE          Geodatabase of Lithologies for the
                                         Geology calculations [required]
  -lulc, --landcover_name STRING         Name of the landcover to be used
                                         [required]
  -weights, --weights_path DIRECTORY      Path to the weights for each raster layer
                                          Geology, Slope, Aspect, Elevation, Landuse
                                          and Hydro to be used for the calculation of
                                          the LSI. [required]
  -epsg, --epsg_code INTEGER              EPSG Code
  -out, --output_path DIRECTORY           Path to output directory.  [required]
  --ftep BOOLEAN


  -h, --help                              Show this message and exit.
```


Example for running the tool from the command line:

Global

```shell
python lsi.py -aoi "\\path\to\aoi_brazil.shp" -loc "Global" -dem  "COPDEM 30m" -litho "\\path\to\lithology.gdb" -lulc "ESA WorldCover - 2021 (10m)" -weights "\\path\to\weights_outside" -out ""\\path\to\output""
```

Europe

```shell
python lsi.py -aoi "\\path\to\aoi_austria.shp" -loc "Europe" -dem  "COPDEM 30m" -litho "\\path\to\lithology.gdb" -lulc "Corine Land Cover - 2018 (100m)" -weights "\\path\to\weights_europe" -out ""\\path\to\output""
```

