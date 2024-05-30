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

The minimum requirement for the LSI calculation is:
* AOI
* DEM
* Lithology geodatabase
* Weights geodatabases
* Landcover

# ArcGIS Pro inputs:

TODO

## Run from the command line

```text
Usage: lsi.py [OPTIONS]

  [Description TODO]

Options:
--aoi             -aoi      PATH                   Path to the AOI (shp, 
                                                    geojson) or WKT       
                                                    string                
                                                    [required]            
--location        -loc      [Europe|Global]        Location of the AOI   
                                                    [required]            
--dem_name        -dem      [COPDEM                DEM Name needed       
                            30m|FABDEM|SRTM        [default: COPDEM 30m] 
                            30m|Other]                                   
--other_dem       -demp     PATH                   DEM path if dem =     
                                                    Other                 
--lithology_gdb   -litho    PATH                   GDB of lithologies.
                                                   [required]   
--landcover_name  -lc       [ESA WorldCover -      Land Cover Name       
                            2021 (10m)|Corine      [default: ESA         
                            Land Cover - 2018      Worldcover - 2021     
                            (100m)]                (10m)]                
--weights_path    -weights  PATH                   Geotadabase with the  
                                                    weights for the LSI   
                                                    computation.
                                                    [required]          
--epsg_code       -epsg     INTEGER RANGE          EPSG code, 4326 is    
                            [1024<=x<=32767]       not accepted. By      
                                                    default, it is the    
                                                    EPSG code of the AOI  
                                                    UTM zone.             
--output_path     -out      DIRECTORY              Output directory.     
                                                    [required]            
--ftep                                             Set this flag if the  
                                                    command line is run   
                                                    on the ftep platform. 
                                                    [default: False]      
--help            -h                               Show this message and 
                                                    exit.

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


