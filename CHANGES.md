# Release History

## 0.1.0 (2024-MM-DD)

- :rocket: First release

## 1.2.1 (2024-06-01)
- Replace flake8, isort and black to ruff
- The tag for sertit_atools changed from "RRM" to "CEMS RRM"

## 1.3.0 (2025-01-17)
- CI added

## 1.3.1 (2025-01-27)
- Change in reclass_landcover and reclass_landcover_elsus from continous xr.where to dict comprehension for better readeability.


## 1.3.3 (2025-04-22)
- Add ftep option for hydro layer computation at GLOBAL method. This method uses pysheds for filling depressions to avoid panicking error at the FTEP.
- Delete the value_field input at the flow accumulation polyline rasterization to keep it binary

## 1.3.4 (2025-04-23)
- Fix numba debugging in CI.
- Fix typing error in dependencies for whitebox workflows.
