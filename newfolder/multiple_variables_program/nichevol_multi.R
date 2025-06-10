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

# load data
data("occ_list", package = "nichevol")
data("tree", package = "nichevol")

# get niche model files
m_files <- list.files(system.file("extdata", package = "nichevol"), 
                      pattern="m\\d.gpkg", full.names=TRUE)
m_list <- lapply(m_files, terra::vect)

# load temperature and precipitation rasters
temp <- rast(system.file("extdata", "temp.tif", package = "nichevol"))
# Use the new precipitation file path
datadir <- "newfolder/climate_data/precip.tif"
precip <- rast(datadir)

# set up color palette
color_palette <- viridis::viridis(20, option = "plasma")

# create multi-dimensional bins for temperature and precipitation
bin_tabl_temp <- bin_table(Ms = m_list, occurrences = occ_list, species = "species", 
                           longitude = "x", latitude = "y", variable = temp, 
                           percentage_out = 5, n_bins = 20)

bin_tabl_precip <- bin_table(Ms = m_list, occurrences = occ_list, species = "species", 
                             longitude = "x", latitude = "y", variable = precip, 
                             percentage_out = 5, n_bins = 20)

# match tree tip labels to data
tree$tip.label <- rownames(bin_tabl_temp)
tree <- ladderize(tree)

# combine tree and environmental data
tree_data_temp <- treedata(tree, bin_tabl_temp)
tree_data_precip <- treedata(tree, bin_tabl_precip)

# ancestral reconstruction using maximum likelihood
table_ml_rec_temp <- bin_ml_rec(tree_data_temp)
s_ml_rec_table_temp <- smooth_rec(table_ml_rec_temp)

table_ml_rec_precip <- bin_ml_rec(tree_data_precip)
s_ml_rec_table_precip <- smooth_rec(table_ml_rec_precip)

# function to convert bin probabilities to pie chart data
create_pie_data <- function(node_data) {
  node_data <- as.numeric(node_data)
  node_data[is.na(node_data)] <- 0
  probs <- node_data / sum(node_data)
  return(probs)
}

# prepare pie chart data
bin_tabl_matrix_temp <- as.matrix(bin_tabl_temp)
s_ml_rec_matrix_temp <- as.matrix(s_ml_rec_table_temp)

bin_tabl_matrix_precip <- as.matrix(bin_tabl_precip)
s_ml_rec_matrix_precip <- as.matrix(s_ml_rec_table_precip)

bin_tabl_matrix_temp[is.na(bin_tabl_matrix_temp)] <- 0
s_ml_rec_matrix_temp[is.na(s_ml_rec_matrix_temp)] <- 0

bin_tabl_matrix_precip[is.na(bin_tabl_matrix_precip)] <- 0
s_ml_rec_matrix_precip[is.na(s_ml_rec_matrix_precip)] <- 0

tip_pies_temp <- t(apply(bin_tabl_matrix_temp, 1, create_pie_data))
node_pies_temp <- t(apply(s_ml_rec_matrix_temp, 1, create_pie_data))

tip_pies_precip <- t(apply(bin_tabl_matrix_precip, 1, create_pie_data))
node_pies_precip <- t(apply(s_ml_rec_matrix_precip, 1, create_pie_data))

# create pie chart visualization for temperature
png("ancestral_reconstruction_visualization_temp.png", width = 1200, height = 800, res = 100)

par(mar = c(5, 4, 4, 8))

tree_plot <- plot.phylo(tree, label.offset = 0.04, type = "phylogram", 
                        direction = "rightwards", show.tip.label = TRUE, 
                        cex = 0.7, edge.width = 2)

lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)

# add pie charts at tips
for(i in 1:nrow(tip_pies_temp)) {
  x <- lastPP$xx[i]
  y <- lastPP$yy[i]
  plotrix::floating.pie(x, y, tip_pies_temp[i,], radius = 0.05, 
                       col = color_palette, border = "white", lwd = 0.5)
}

# add pie charts at nodes
for(i in 1:nrow(node_pies_temp)) {
  x <- lastPP$xx[length(tree$tip.label) + i]
  y <- lastPP$yy[length(tree$tip.label) + i]
  plotrix::floating.pie(x, y, node_pies_temp[i,], radius = 0.05, 
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

# create pie chart visualization for precipitation
png("ancestral_reconstruction_visualization_precip.png", width = 1200, height = 800, res = 100)

par(mar = c(5, 4, 4, 8))

tree_plot <- plot.phylo(tree, label.offset = 0.04, type = "phylogram", 
                        direction = "rightwards", show.tip.label = TRUE, 
                        cex = 0.7, edge.width = 2)

lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)

# add pie charts at tips
for(i in 1:nrow(tip_pies_precip)) {
  x <- lastPP$xx[i]
  y <- lastPP$yy[i]
  plotrix::floating.pie(x, y, tip_pies_precip[i,], radius = 0.05, 
                       col = color_palette, border = "white", lwd = 0.5)
}

# add pie charts at nodes
for(i in 1:nrow(node_pies_precip)) {
  x <- lastPP$xx[length(tree$tip.label) + i]
  y <- lastPP$yy[length(tree$tip.label) + i]
  plotrix::floating.pie(x, y, node_pies_precip[i,], radius = 0.05, 
                       col = color_palette, border = "white", lwd = 0.5)
}

# add legend
legend("right", 
       legend = paste("Bin", 1:20),
       fill = color_palette,
       cex = 0.6,
       ncol = 2,
       title = "Precipitation Bins",
       bty = "n",
       xpd = TRUE,
       inset = c(-0.2, 0))

title(main = "Ancestral Reconstruction of Precipitation Niches",
      sub = "Maximum Likelihood Estimation with Enhanced Visualization",
      cex.main = 1.2,
      cex.sub = 0.8)

dev.off()

cat("\nenhanced visualizations created:\n")
cat("1. ancestral_reconstruction_visualization_temp.png\n")
cat("2. ancestral_reconstruction_visualization_precip.png\n")

# calculate continuous trait values (weighted mean of bins)
get_trait_value <- function(row) {
  bins <- as.numeric(row)
  bins[is.na(bins)] <- 0
  sum(bins * seq_along(bins)) / sum(bins)
}

tip_trait_temp <- apply(bin_tabl_matrix_temp, 1, get_trait_value)
node_trait_temp <- apply(s_ml_rec_matrix_temp, 1, get_trait_value)

tip_trait_precip <- apply(bin_tabl_matrix_precip, 1, get_trait_value)
node_trait_precip <- apply(s_ml_rec_matrix_precip, 1, get_trait_value)

# ensure tip_trait is ordered according to tree tip labels
tip_trait_ordered_temp <- tip_trait_temp[tree$tip.label]
tip_trait_ordered_precip <- tip_trait_precip[tree$tip.label]

# create continuous trait visualization for temperature
png("ancestral_state_contmap_temp.png", width = 1200, height = 800, res = 100)

par(mfrow=c(1,1))
par(mar=c(8,4,4,2)+0.1)

# create and plot continuous trait map
contmap_obj_temp <- phytools::contMap(tree, tip_trait_ordered_temp, plot=FALSE, type="phylogram", direction="rightwards")
contmap_obj_temp$cols <- colorRampPalette(c("red", "orange", "yellow", "green", "blue"))(1000)

plot(contmap_obj_temp,
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

min_trait_val <- min(tip_trait_ordered_temp, na.rm = TRUE)
max_trait_val <- max(tip_trait_ordered_temp, na.rm = TRUE)

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
                      cols = contmap_obj_temp$cols,
                      gradient = "x",
                      align = "rb",
                      cex = 0.8,
                      xpd = TRUE,
                      rect.col = NA)

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

# create continuous trait visualization for precipitation
png("ancestral_state_contmap_precip.png", width = 1200, height = 800, res = 100)

par(mfrow=c(1,1))
par(mar=c(8,4,4,2)+0.1)

# create and plot continuous trait map
contmap_obj_precip <- phytools::contMap(tree, tip_trait_ordered_precip, plot=FALSE, type="phylogram", direction="rightwards")
contmap_obj_precip$cols <- colorRampPalette(c("red", "orange", "yellow", "green", "blue"))(1000)

plot(contmap_obj_precip,
     legend=FALSE,
     fsize=0.8,
     outline=TRUE,
     lwd=4,
     type="phylogram",
     direction="rightwards")

title(main = "Ancestral State Reconstruction: Precipitation",
      sub = "Continuous Trait Evolution with Enhanced Visualization",
      cex.main = 1.2,
      cex.sub = 0.8)

# add custom horizontal legend
usr <- par("usr")
x_range <- usr[2] - usr[1]
y_range <- usr[4] - usr[3]

min_trait_val <- min(tip_trait_ordered_precip, na.rm = TRUE)
max_trait_val <- max(tip_trait_ordered_precip, na.rm = TRUE)

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
                      cols = contmap_obj_precip$cols,
                      gradient = "x",
                      align = "rb",
                      cex = 0.8,
                      xpd = TRUE,
                      rect.col = NA)

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

cat("3. ancestral_state_contmap_temp.png\n")
cat("4. ancestral_state_contmap_precip.png\n") 
