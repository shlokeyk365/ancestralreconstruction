# Script to download and prepare WorldClim precipitation data
library(geodata)
library(terra)

# Set the working directory to the climate_data folder
setwd("newfolder/climate_data")

# Download WorldClim data (precipitation)
precip_data <- worldclim_global(var = "prec", res = 10, path = ".")

# Calculate mean annual precipitation
mean_precip <- mean(precip_data)

# Save the mean precipitation raster
writeRaster(mean_precip, "precip.tif", overwrite = TRUE)

# Print information about the downloaded data
print("Precipitation data downloaded and processed successfully!")
print(paste("Resolution:", res(mean_precip)))
print(paste("Extent:", ext(mean_precip)))
print(paste("CRS:", crs(mean_precip)))

# Create a simple visualization
png("precipitation_map.png", width = 1200, height = 800)
plot(mean_precip, main = "Mean Annual Precipitation (mm)")
dev.off()

print("Visualization saved as 'precipitation_map.png'") 
