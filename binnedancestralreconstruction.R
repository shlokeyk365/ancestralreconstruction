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
library(viridis)   # For better color palettes
library(ggplot2)   # For enhanced visualization capabilities

# Enhanced error handling function (commented out for now)
# validate_input_data <- function(occ_list, m_list, temp, tree) {
#   # Check if all required data is present
#   if (is.null(occ_list)) stop("Species occurrence data is missing")
#   if (is.null(m_list)) stop("Niche model data is missing")
#   if (is.null(temp)) stop("Temperature data is missing")
#   if (is.null(tree)) stop("Phylogenetic tree is missing")
#   
#   # Check data structure
#   if (!all(c("species", "x", "y") %in% names(occ_list))) {
#     stop("Occurrence data must contain 'species', 'x', and 'y' columns")
#   }
#   
#   # Check for missing values
#   if (any(is.na(occ_list))) {
#     warning("Missing values found in occurrence data")
#   }
#   
#   # Validate tree structure
#   if (!inherits(tree, "phylo")) {
#     stop("Tree must be a phylogenetic tree object")
#   }
# }

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

# Enhanced color palette using viridis
color_palette <- viridis::viridis(20, option = "plasma")
# Alternative color palettes
color_palette <- viridis::magma(20)
color_palette <- viridis::inferno(20)
color_palette <- colorRampPalette(c("darkblue", "blue", "green", "yellow", "red", "darkred"))(20)

# Generate a table dividing temperature data into bins (ranges) for each species
# This allows for a more structured comparison of environmental conditions
bin_tabl <- bin_table(Ms = m_list, occurrences = occ_list, species = "species", 
                      longitude = "x", latitude = "y", variable = temp, 
                      percentage_out = 5, n_bins = 20)  # Divides temperature values into 20 bins

# Ensure that the species names in the phylogenetic tree match the data table
tree$tip.label <- rownames(bin_tabl)

# Ladderize the tree to improve plotting layout
tree <- ladderize(tree)

# Prepare data for ancestral state reconstruction by combining tree and environmental data
tree_data <- treedata(tree, bin_tabl)

# Perform ancestral reconstruction using the Maximum Likelihood method
# This method estimates the most probable evolutionary changes
table_ml_rec <- bin_ml_rec(tree_data)

# Apply smoothing to Maximum Likelihood results to make trends clearer
s_ml_rec_table <- smooth_rec(table_ml_rec)

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

# Enhanced visualization with better formatting
png("ancestral_reconstruction_visualization.png", width = 1200, height = 800, res = 100)

# Set up the plot with better margins
par(mar = c(5, 4, 4, 8))

# Plot the tree with enhanced formatting
tree_plot <- plot.phylo(tree, label.offset = 0.04, type = "phylogram", 
                        direction = "rightwards", show.tip.label = TRUE, 
                        cex = 0.7, edge.width = 2)

# Get coordinates
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)

# Add pie charts at tips with enhanced formatting
for(i in 1:nrow(tip_pies)) {
  x <- lastPP$xx[i]
  y <- lastPP$yy[i]
  plotrix::floating.pie(x, y, tip_pies[i,], radius = 0.05, 
                       col = color_palette, border = "white", lwd = 0.5)
}

# Add pie charts at nodes with enhanced formatting
for(i in 1:nrow(node_pies)) {
  x <- lastPP$xx[length(tree$tip.label) + i]
  y <- lastPP$yy[length(tree$tip.label) + i]
  plotrix::floating.pie(x, y, node_pies[i,], radius = 0.05, 
                       col = color_palette, border = "white", lwd = 0.5)
}

# Enhanced legend with better formatting
legend("right", 
       legend = paste("Bin", 1:20),
       fill = color_palette,
       cex = 0.6,
       ncol = 2,
       title = "Temperature Bins",
       bty = "n",
       xpd = TRUE,
       inset = c(-0.2, 0))

# Add title and subtitle
title(main = "Ancestral Reconstruction of Temperature Niches",
      sub = "Maximum Likelihood Estimation with Enhanced Visualization",
      cex.main = 1.2,
      cex.sub = 0.8)

dev.off()

# Print confirmation messages
cat("\nEnhanced visualizations have been created:\n")
cat("1. ancestral_reconstruction_visualization.png\n")

# Calculate a continuous trait value for each tip and node (e.g., weighted mean of bins)
get_trait_value <- function(row) {
  bins <- as.numeric(row)
  bins[is.na(bins)] <- 0
  # Weighted mean: bin index * probability
  sum(bins * seq_along(bins)) / sum(bins)
}

# For tips
tip_trait <- apply(bin_tabl_matrix, 1, get_trait_value)
# For nodes (ancestral)
node_trait <- apply(s_ml_rec_matrix, 1, get_trait_value)

# Ensure tip_trait is ordered according to tree$tip.label
tip_trait_ordered <- tip_trait[tree$tip.label]

# Enhanced continuous trait visualization
png("ancestral_state_contmap.png", width = 1200, height = 800, res = 100)

# Reset graphical parameters for a clean plot
par(mfrow=c(1,1)) # Set to single plot layout
par(mar=c(5,4,4,2)+0.1) # Default margins plus a little padding

# Create and plot the continuous trait map with enhanced formatting
# Pass the ordered tip_trait to contMap and explicitly set the type and direction
contmap_obj <- phytools::contMap(tree, tip_trait_ordered, plot=FALSE, type="phylogram", direction="rightwards")
# Set the color palette for the contMap to match the example image (red-orange-yellow-green-blue gradient)
contmap_obj$cols <- colorRampPalette(c("red", "yellow", "green", "blue"))(1000) # Custom color ramp for red-to-blue

# Plot the contMap with enhanced visibility parameters and simplified label handling
plot(contmap_obj, 
     legend=0.7*max(nodeHeights(tree)), 
     fsize=0.8, # Slightly increased font size for better visibility of tip labels
     outline=TRUE, 
     lwd=2,
     type="phylogram",  
     direction="rightwards")

# Add title and subtitle
title(main = "Ancestral State Reconstruction: Temperature",
      sub = "Continuous Trait Evolution with Enhanced Visualization",
      cex.main = 1.2,
      cex.sub = 0.8)

dev.off()

cat("\nEnhanced visualizations have been created:\n")
cat("2. ancestral_state_contmap.png\n")

# Data format validation function (commented out)
# validate_data_format <- function(data, expected_format) {
#   if (!inherits(data, expected_format)) {
#     stop(sprintf("Data must be of class %s", expected_format))
#   }
# }

# Missing data handling function (commented out)
# handle_missing_data <- function(data) {
#   if (any(is.na(data))) {
#     warning("Missing values detected in data")
#     # Implement appropriate handling strategy
#     # data <- na.omit(data)  # or other appropriate method
#   }
#   return(data)
# }

# --- Diagnostic Code for contMap (Commented Out) ---
# This section is for troubleshooting if contMap visualization is still problematic.
# Uncomment and run this section to test basic contMap functionality.

# # Simulate a simple tree and trait data
# set.seed(123) # for reproducibility
# simple_tree <- rtree(5)
# simple_trait <- setNames(rnorm(Ntip(simple_tree)), simple_tree$tip.label)

# # Create a simple contMap object
# simple_contmap_obj <- phytools::contMap(simple_tree, simple_trait, plot=FALSE, type="phylogram", direction="rightwards")
# simple_contmap_obj$cols <- viridis::viridis(100) # Use a simple color palette

# # Plot the simple contMap
# png("simple_contmap_diagnostic.png", width = 800, height = 600, res = 100)
# par(mfrow=c(1,1))
# par(mar=c(5,4,4,2)+0.1)
# plot(simple_contmap_obj, 
#      legend=0.7*max(nodeHeights(simple_tree)), 
#      fsize=0.8, 
#      outline=FALSE,
#      lwd=3,
#      type="phylogram", 
#      direction="rightwards")
# title(main = "Simple contMap Diagnostic Plot", cex.main = 1.2)
# dev.off()
# cat("\nSimple contMap diagnostic plot created: simple_contmap_diagnostic.png\n")

# # Check tree properties (for your actual tree)
# cat("\n--- Original Tree Properties ---\n")
# cat("Is tree rooted? ", is.rooted(tree), "\n")
# cat("Is tree ultrametric? ", is.ultrametric(tree), "\n")
# cat("Has edge lengths? ", !is.null(tree$edge.length), "\n")
# cat("Number of tips: ", Ntip(tree), "\n")
# cat("Number of nodes: ", Nnode(tree), "\n")

# --- End Diagnostic Code ---
