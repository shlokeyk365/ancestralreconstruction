# Load necessary libraries for different functionalities
if (!requireNamespace("plotrix", quietly = TRUE)) {
  install.packages("plotrix")
}
library(nichevol)  # Provides tools for ecological niche evolution analysis
library(terra)     # Handles spatial and environmental data, such as temperature maps
library(ape)       # Used for working with phylogenetic trees
library(geiger)    # Provides functions for evolutionary biology analysis
library(phytools)  # Additional tools for manipulating and analyzing phylogenetic trees
library(plotrix)   # Required for floating.pie function

# Load species occurrence data from the 'nichevol' package
data("occ_list", package = "nichevol")

# Find all model files in the package's extdata folder that describe ecological niches
m_files <- list.files(system.file("extdata", package = "nichevol"), 
                      pattern="m\\d.gpkg", full.names=TRUE)

# Read these niche model files as spatial vector objects
# These models define the regions accessible to species
m_list <- lapply(m_files, terra::vect)

# Load a raster file containing temperature data for the study area
# Raster files store environmental data across geographical areas
temp <- rast(system.file("extdata", "temp.tif", package = "nichevol"))

# Load a phylogenetic tree dataset that represents evolutionary relationships between species
data("tree", package = "nichevol")

# Generate a table dividing temperature data into bins (ranges) for each species
# This allows for a more structured comparison of environmental conditions
bin_tabl <- bin_table(Ms = m_list, occurrences = occ_list, species = "species", 
                      longitude = "x", latitude = "y", variable = temp, 
                      percentage_out = 5, n_bins = 20)  # Divides temperature values into 20 bins

# Ensure that the species names in the phylogenetic tree match the data table
tree$tip.label <- rownames(bin_tabl)

# Prepare data for ancestral state reconstruction by combining tree and environmental data
tree_data <- treedata(tree, bin_tabl)

# Perform ancestral reconstruction using the Maximum Likelihood method
# This method estimates the most probable evolutionary changes
table_ml_rec <- bin_ml_rec(tree_data)

# Apply smoothing to Maximum Likelihood results to make trends clearer
s_ml_rec_table <- smooth_rec(table_ml_rec)

# Create a color palette for the pie charts
# We'll use a color gradient from cool to warm temperatures
color_palette <- colorRampPalette(c("blue", "green", "yellow", "red"))(20)

# Function to convert bin probabilities to pie chart data
create_pie_data <- function(node_data) {
  # Convert to numeric if not already, handling NAs
  node_data <- as.numeric(node_data)
  # Replace NAs with 0
  node_data[is.na(node_data)] <- 0
  # Normalize the probabilities to sum to 1
  probs <- node_data / sum(node_data)
  return(probs)
}

# Create pie chart data for tips and nodes
bin_tabl_matrix <- as.matrix(bin_tabl)
s_ml_rec_matrix <- as.matrix(s_ml_rec_table)

# Handle any remaining NAs in the matrices
bin_tabl_matrix[is.na(bin_tabl_matrix)] <- 0
s_ml_rec_matrix[is.na(s_ml_rec_matrix)] <- 0

tip_pies <- t(apply(bin_tabl_matrix, 1, create_pie_data))
node_pies <- t(apply(s_ml_rec_matrix, 1, create_pie_data))

# Create the first visualization (pie chart tree)
png("ancestral_reconstruction_visualization.png", width = 1200, height = 800, res = 100)

# Plot the tree and get the coordinates
tree_plot <- plot.phylo(tree, label.offset = 0.04, type = "phylogram", 
                        direction = "rightwards", show.tip.label = TRUE, cex = 0.7)
# Get the coordinates of the tips and nodes
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)

# Add pie charts at tips
for(i in 1:nrow(tip_pies)) {
  x <- lastPP$xx[i]
  y <- lastPP$yy[i]
  plotrix::floating.pie(x, y, tip_pies[i,], radius = 0.05, col = color_palette, border = "white")
}

# Add pie charts at nodes (internal nodes)
for(i in 1:nrow(node_pies)) {
  x <- lastPP$xx[length(tree$tip.label) + i]
  y <- lastPP$yy[length(tree$tip.label) + i]
  plotrix::floating.pie(x, y, node_pies[i,], radius = 0.05, col = color_palette, border = "white")
}

# Add a legend
legend("bottomright", 
       legend = paste("Bin", 1:20),
       fill = color_palette,
       cex = 0.6,
       ncol = 4,
       title = "Temperature Bins")

# Add a title
title("Ancestral Reconstruction of Temperature Niches\nwith Maximum Likelihood Estimation")

dev.off()

# Print confirmation message
cat("\nVisualization has been created:\n")
cat("ancestral_reconstruction_visualization.png\n")

# Calculate a continuous trait value for each tip and node (e.g., weighted mean of bins)
get_trait_value <- function(row) {
  bins <- as.numeric(row)
  bins[is.na(bins)] <- 0
  # Avoid division by zero
  if(sum(bins) == 0) return(NA)
  # Weighted mean: bin index * probability
  sum(bins * seq_along(bins)) / sum(bins)
}

# For tips
tip_trait <- apply(bin_tabl_matrix, 1, get_trait_value)

# Check for NAs in tip_trait and handle them
if(any(is.na(tip_trait))) {
  cat("Warning: Some tip trait values are NA. Replacing with mean value.\n")
  tip_trait[is.na(tip_trait)] <- mean(tip_trait, na.rm = TRUE)
}

# Make sure tip_trait is a named vector
names(tip_trait) <- tree$tip.label

# Debug: Print some information
cat("\nDebugging information:\n")
cat("Number of tips:", length(tree$tip.label), "\n")
cat("Length of tip_trait:", length(tip_trait), "\n")
cat("Range of tip_trait:", range(tip_trait, na.rm = TRUE), "\n")
cat("Any NAs in tip_trait:", any(is.na(tip_trait)), "\n")

# Try to create the contMap object with error handling
tryCatch({
  cat("\nAttempting to create contMap object...\n")
  contmap_obj <- phytools::contMap(tree, tip_trait, plot=FALSE)
  cat("contMap object created successfully!\n")
  
  # Save the plot as a PNG
  png("ancestral_state_contmap.png", width = 1200, height = 800, res = 100)
  plot(contmap_obj, legend=0.7*max(nodeHeights(tree)), fsize=0.5, outline=FALSE)
  title("Ancestral State Reconstruction: Temperature")
  dev.off()
  
  cat("\ncontMap-style visualization has been created:\n")
  cat("ancestral_state_contmap.png\n")
  
}, error = function(e) {
  cat("Error creating contMap:\n")
  cat(e$message, "\n")
  
  # Alternative approach using a simpler method
  cat("\nTrying alternative approach with ace() function...\n")
  
  # Use ace() for ancestral state estimation
  ace_result <- ape::ace(tip_trait, tree, type = "continuous")
  
  # Create a simple color-coded tree plot
  png("ancestral_state_alternative.png", width = 1200, height = 800, res = 100)
  
  # Plot tree with branch colors based on trait values
  plot(tree, show.tip.label = TRUE, cex = 0.7)
  
  # Add colored points at tips based on trait values
  tip_colors <- heat.colors(100)[as.numeric(cut(tip_trait, breaks = 100))]
  tiplabels(pch = 21, bg = tip_colors, cex = 1.5)
  
  # Add colored points at nodes based on ancestral estimates
  node_colors <- heat.colors(100)[as.numeric(cut(ace_result$ace, breaks = 100))]
  nodelabels(pch = 21, bg = node_colors, cex = 1.2)
  
  title("Ancestral State Reconstruction: Temperature (Alternative Method)")
  
  # Add a simple legend
  legend("bottomright", 
         legend = c("Low Temp", "High Temp"),
         fill = c("red", "yellow"),
         title = "Temperature")
  
  dev.off()
  
  cat("Alternative visualization created: ancestral_state_alternative.png\n")
})
