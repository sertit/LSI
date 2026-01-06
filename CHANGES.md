# Release History

## 1.4.2 (2026-MM-DD)
- CI: Use new groupware CI template for triggering sertit_atools update
- CI: Add pre-commit bot
- FIX: Float Nodata to Int for DEM at hydro layer

## 1.4.1 (2025-05-23)
- ENH: Change the print errors to LOGGER errors to catch and explain the FTEP S3 reading errors

## 1.4.0 (2025-05-06)
- ENH: Add ftep=False as input to arcgispro toolbox
- ENH: Delete clip of AOI by coastline, instead clip by GADM layer which is better aligned with coasts.
- ENH: Add exception control for error in reading rasters with rasterio and vectors with pyogrio at the FTEP, assumed to be related with problems of timeouts, failed reads, networking issues, rate limits, etc.
- CI: Change GLC to CLC for european test.
- DOC: Modify in CHANGES.md organization.

## 1.3.7 (2025-04-24)
- Add the needed ftep as input for lsi_core in lsi_ftep.py

## 1.3.6 (2025-04-23)
- Fix typo

## 1.3.5 (2025-04-23)
- Fix typo
- Fix pyproject

## 1.3.4 (2025-04-23)
- Fix numba debugging in CI.
- Fix typing error in dependencies for whitebox workflows.

## 1.3.3 (2025-04-22)
- Add ftep option for hydro layer computation at GLOBAL method. This method uses pysheds for filling depressions to avoid panicking error at the FTEP.
- Delete the value_field input at the flow accumulation polyline rasterization to keep it binary

## 1.3.1 (2025-01-27)
- Change in reclass_landcover and reclass_landcover_elsus from continous xr.where to dict comprehension for better readeability.

## 1.3.0 (2025-01-17)
- CI added

## 1.2.1 (2024-06-01)
- Replace flake8, isort and black to ruff
- The tag for sertit_atools changed from "RRM" to "CEMS RRM"

## 0.1.0 (2024-MM-DD)

- :rocket: First release