# Multi-Dimensional Binned Ancestral Reconstruction

This program performs ancestral state reconstruction for multiple environmental variables (e.g., temperature and precipitation) using a binned approach on a phylogenetic tree. It is designed for ecological and evolutionary studies where understanding the historical niche evolution of species is important.

## Features
- Supports multiple environmental variables (currently temperature and precipitation)
- Uses maximum likelihood for ancestral state reconstruction
- Generates both pie chart and continuous trait visualizations for each variable
- Handles raster data and species occurrence data

## Dependencies
- R (>= 4.0)
- R packages: `nichevol`, `terra`, `ape`, `geiger`, `phytools`, `plotrix`, `viridis`, `ggplot2`
- For downloading climate data: `geodata`

## Input Data
- Species occurrence data (from the `nichevol` package example)
- Niche model files (from the `nichevol` package example)
- Phylogenetic tree (from the `nichevol` package example)
- Environmental rasters:
  - `temp.tif` (from `nichevol` package)
  - `precip.tif` (downloaded using the provided script in `newfolder/climate_data/`)

## How to Run
1. **Download Precipitation Data**
   - Run `download_climate_data.R` to download and prepare the `precip.tif` file:
     ```sh
     Rscript newfolder/download_climate_data.R
     ```
2. **Run the Ancestral Reconstruction**
   - Execute the main script:
     ```sh
     Rscript newfolder/binnedancestralreconstruction_multi.R
     ```

## Outputs
- `ancestral_reconstruction_visualization_temp.png`: Pie chart visualization of temperature bins on the tree
- `ancestral_reconstruction_visualization_precip.png`: Pie chart visualization of precipitation bins on the tree
- `ancestral_state_contmap_temp.png`: Continuous trait map for temperature
- `ancestral_state_contmap_precip.png`: Continuous trait map for precipitation
- `precipitation_map.png`: Map of mean annual precipitation (from the download script)

## Notes
- The script expects the precipitation raster at `newfolder/climate_data/precip.tif`. Adjust the path in the script if you move the file.
- Warnings about NAs may appear if some species or bins have missing data; these do not affect the main outputs.
- You can extend the script to include additional environmental variables by following the same pattern.

## Citation
If you use this workflow, please cite the `nichevol` package and any climate data sources (e.g., WorldClim).

---
For questions or suggestions, please contact the repository maintainer. 
