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

There are two available methods (Global and Europe). The Global method is based on the use of all previous data sources, the Europe method on the other side only considers the following data:
* Geology/Lithology
* Terrain Slope
* Land use and land cover (LULC)

> The Europe method uses the division by Climate-physiographically differentiated zones used for the susceptibility ELSUS map.
> For more information: [ELSUS susceptibility map](https://www.sciencedirect.com/science/article/abs/pii/S0169555X14003675).

The minimum requirement for the LSI calculation is:
* AOI
* DEM
* Landcover
* Output path

# ArcGIS Pro inputs:

TODO

## Run from the command line

```text
Usage: lsi.py [OPTIONS]

  [Description TODO]

+- Options -------------------------------------------------------------------+
| *  --aoi               -aoi        PATH                 AOI (shp, geojson)  |
|                                                         or WKT string       |
|                                                         [required]          |
|    --location          -loc        [Europe|Global]      Location of the AOI |
|                                                         [default: Global]   |
|    --dem_name          -dem        [COPDEM              DEM Name needed     |
|                                    30m|FABDEM|SRTM      [default: COPDEM    |
|                                    30m|Other]           30m]                |
|    --other_dem         -demp       PATH                 DEM path if dem =   |
|                                                         Other               |
|    --landcover_name    -lulc       [ESA WorldCover -    Land Cover Name     |
|                                    2021 (10m)|Corine    [default: ESA       |
|                                    Land Cover - 2018    WorldCover - 2021   |
|                                    (100m)]              (10m)]              |
|    --europe_method     -eu_method  [Refined|Fast]       if LOCATION =       |
|                                                         EUROPE, choose      |
|                                                         whether you want a  |
|                                                         fast computation    |
|                                                         with lower          |
|                                                         resolution (based   |
|                                                         on the pre-existent |
|                                                         ELSUS layer) or a   |
|                                                         refined LSI         |
|                                                         computation         |
|                                                         [default: Refined]  |
|    --output_resoluti▒  -res        INTEGER RANGE        Output resolution.  |
|                                    [1<=x<=1000]         Taking from DEM if  |
|                                                         not provided        |
|    --epsg_code         -epsg       INTEGER RANGE        EPSG code, 4326 is  |
|                                    [1024<=x<=32767]     not accepted. By    |
|                                                         default, it is the  |
|                                                         EPSG code of the    |
|                                                         AOI UTM zone.       |
| *  --output_path       -out        DIRECTORY            Output directory.   |
|                                                         [required]          |
|    --ftep                                               Set this flag if    |
|                                                         the command line is |
|                                                         run on the ftep     |
|                                                         platform.           |
|                                                         [default: False]    |
|    --help              -h                               Show this message   |
|                                                         and exit.           |
+-----------------------------------------------------------------------------+


```


Example for running the tool from the command line:

Global

```shell
python lsi.py -aoi "\\path\to\aoi_venezuela.shp" -loc "Global" -dem "COPDEM 30m" -lulc "ESA Worldcover - 2021 (10m)" -res 30 -out "\\path\to\output_folder"
```

Europe

```shell
python lsi.py -aoi "\\path\to\aoi_grenoble_france.shp" -loc "Europe" -dem "COPDEM 30m" -lulc "Corine Land Cover - 2018 (100m)" -eu_method "Refined" -res 30 -out "\\path\to\output_folder"
```


