# load necessary libraries
if (!requireNamespace("plotrix", quietly = TRUE)) {
  install.packages("plotrix")
}
library(nichevol)
library(terra)
library(ape)
library(geiger)
library(phytools)
library(plotrix)
library(viridis)
library(ggplot2)

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

# load data
data("occ_list", package = "nichevol")
data("tree", package = "nichevol")

# get niche model files
m_files <- list.files(system.file("extdata", package = "nichevol"), 
                      pattern="m\\d.gpkg", full.names=TRUE)
m_list <- lapply(m_files, terra::vect)

# load temperature raster
temp <- rast(system.file("extdata", "temp.tif", package = "nichevol"))

# set up color palette
color_palette <- viridis::viridis(20, option = "plasma")
# alternative palettes:
# color_palette <- viridis::magma(20)
# color_palette <- viridis::inferno(20)
# color_palette <- colorRampPalette(c("darkblue", "blue", "green", "yellow", "red", "darkred"))(20)

# create temperature bins for each species
bin_tabl <- bin_table(Ms = m_list, occurrences = occ_list, species = "species", 
                      longitude = "x", latitude = "y", variable = temp, 
                      percentage_out = 5, n_bins = 20)

# match tree tip labels to data
tree$tip.label <- rownames(bin_tabl)
tree <- ladderize(tree)

# combine tree and environmental data
tree_data <- treedata(tree, bin_tabl)

# ancestral reconstruction using maximum likelihood
table_ml_rec <- bin_ml_rec(tree_data)
s_ml_rec_table <- smooth_rec(table_ml_rec)

# function to convert bin probabilities to pie chart data
create_pie_data <- function(node_data) {
  node_data <- as.numeric(node_data)
  node_data[is.na(node_data)] <- 0
  probs <- node_data / sum(node_data)
  return(probs)
}

# prepare pie chart data
bin_tabl_matrix <- as.matrix(bin_tabl)
s_ml_rec_matrix <- as.matrix(s_ml_rec_table)

bin_tabl_matrix[is.na(bin_tabl_matrix)] <- 0
s_ml_rec_matrix[is.na(s_ml_rec_matrix)] <- 0

tip_pies <- t(apply(bin_tabl_matrix, 1, create_pie_data))
node_pies <- t(apply(s_ml_rec_matrix, 1, create_pie_data))

# create pie chart visualization
png("ancestral_reconstruction_visualization.png", width = 1200, height = 800, res = 100)

par(mar = c(5, 4, 4, 8))

tree_plot <- plot.phylo(tree, label.offset = 0.04, type = "phylogram", 
                        direction = "rightwards", show.tip.label = TRUE, 
                        cex = 0.7, edge.width = 2)

lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)

# add pie charts at tips
for(i in 1:nrow(tip_pies)) {
  x <- lastPP$xx[i]
  y <- lastPP$yy[i]
  plotrix::floating.pie(x, y, tip_pies[i,], radius = 0.05, 
                       col = color_palette, border = "white", lwd = 0.5)
}

# add pie charts at nodes
for(i in 1:nrow(node_pies)) {
  x <- lastPP$xx[length(tree$tip.label) + i]
  y <- lastPP$yy[length(tree$tip.label) + i]
  plotrix::floating.pie(x, y, node_pies[i,], radius = 0.05, 
                       col = color_palette, border = "white", lwd = 0.5)
}

# add legend
legend("right", 
       legend = paste("Bin", 1:20),
       fill = color_palette,
       cex = 0.6,
       ncol = 2,
       title = "Temperature Bins",
       bty = "n",
       xpd = TRUE,
       inset = c(-0.2, 0))

title(main = "Ancestral Reconstruction of Temperature Niches",
      sub = "Maximum Likelihood Estimation with Enhanced Visualization",
      cex.main = 1.2,
      cex.sub = 0.8)

dev.off()

cat("\nenhanced visualizations created:\n")
cat("1. ancestral_reconstruction_visualization.png\n")

# calculate continuous trait values (weighted mean of bins)
get_trait_value <- function(row) {
  bins <- as.numeric(row)
  bins[is.na(bins)] <- 0
  sum(bins * seq_along(bins)) / sum(bins)
}

tip_trait <- apply(bin_tabl_matrix, 1, get_trait_value)
node_trait <- apply(s_ml_rec_matrix, 1, get_trait_value)

# ensure tip_trait is ordered according to tree tip labels
tip_trait_ordered <- tip_trait[tree$tip.label]

# create continuous trait visualization
png("ancestral_state_contmap.png", width = 1200, height = 800, res = 100)

par(mfrow=c(1,1))
par(mar=c(8,4,4,2)+0.1)

# create and plot continuous trait map
contmap_obj <- phytools::contMap(tree, tip_trait_ordered, plot=FALSE, type="phylogram", direction="rightwards")
contmap_obj$cols <- colorRampPalette(c("red", "orange", "yellow", "green", "blue"))(1000)

plot(contmap_obj,
     legend=FALSE,
     fsize=0.8,
     outline=TRUE,
     lwd=4,
     type="phylogram",
     direction="rightwards")

title(main = "Ancestral State Reconstruction: Temperature",
      sub = "Continuous Trait Evolution with Enhanced Visualization",
      cex.main = 1.2,
      cex.sub = 0.8)

# add custom horizontal legend
usr <- par("usr")
x_range <- usr[2] - usr[1]
y_range <- usr[4] - usr[3]

min_trait_val <- min(tip_trait_ordered, na.rm = TRUE)
max_trait_val <- max(tip_trait_ordered, na.rm = TRUE)

legend_x_start <- usr[1] + x_range * 0.25
legend_x_end <- usr[1] + x_range * 0.75
legend_y_bottom_bar <- usr[3] - y_range * 0.15
legend_y_top_bar <- legend_y_bottom_bar + y_range * 0.02

legend_labels_values <- c(sprintf("%.1f", min_trait_val), sprintf("%.1f", max_trait_val))

plotrix::color.legend(xl = legend_x_start,
                      yb = legend_y_bottom_bar,
                      xr = legend_x_end,
                      yt = legend_y_top_bar,
                      legend = legend_labels_values,
                      cols = contmap_obj$cols,
                      gradient = "x",
                      align = "rb",
                      cex = 0.8,
                      xpd = TRUE)

text(x = mean(c(legend_x_start, legend_x_end)),
     y = legend_y_bottom_bar - (y_range * 0.04),
     labels = "trait value",
     cex = 0.8,
     xpd = TRUE)

text(x = mean(c(legend_x_start, legend_x_end)),
     y = legend_y_bottom_bar - (y_range * 0.07),
     labels = "length=0.7",
     cex = 0.8,
     xpd = TRUE)

dev.off()

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
