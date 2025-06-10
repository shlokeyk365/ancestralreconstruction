# Niche Evolution Modeling: Ancestral State Reconstruction

This project uses R to reconstruct and visualize the evolution of ecological niches (specifically, temperature preferences) across a phylogenetic tree of species. The main script, `binnedancestralreconstruction.R`, leverages the `nichevol` package and other R tools to perform this analysis.

## What does the script do?
- Loads example species occurrence, environmental, and phylogenetic data from the `nichevol` package.
- Divides temperature data into bins for each species ("binned ancestral reconstruction").
- Reconstructs the likely temperature preferences of ancestral species using statistical methods.
- Visualizes the results on a phylogenetic tree, showing how temperature preferences may have evolved over time.

## How to run the script
1. **Install R** (if you haven't already): [Download R](https://cran.r-project.org/)
2. **Install required R packages** (the script will try to install any missing ones automatically):
   - `nichevol`, `terra`, `ape`, `geiger`, `phytools`, `plotrix`
3. **Open R or RStudio** and set your working directory to the folder containing `binnedancestralreconstruction.R`.
4. **Run the script**:
   - In RStudio: Open the script and click "Source".
   - In R: Use `source("binnedancestralreconstruction.R")`

The script will automatically load all necessary data from the `nichevol` package. You do not need to provide any additional data files.

## What outputs are generated?
- `ancestral_reconstruction_visualization.png`: A phylogenetic tree with pie charts at each tip and node, showing the probability distribution of temperature preferences for each species and ancestor.
- `ancestral_state_contmap.png`: A phylogenetic tree with a color gradient along the branches, showing the most likely temperature preference for each species and ancestor.

## How to interpret the plots

### 1. Pie Chart Tree (`ancestral_reconstruction_visualization.png`)
- **Tree structure**: Each branch represents evolutionary relationships between species.
- **Pie charts at tips**: Show the temperature preferences of living species. Each slice represents a temperature range (bin), and the size of the slice shows the probability that the species prefers that range.
- **Pie charts at nodes**: Show the reconstructed temperature preferences of ancestral species (ancestors). These are estimated using statistical models.
- **Color legend**: Indicates which colors correspond to which temperature bins (from cooler to warmer).
- **Interpretation**: You can see how temperature preferences may have shifted or diversified as species evolved.

### 2. Gradient Tree (`ancestral_state_contmap.png`)
- **Tree structure**: As above, shows evolutionary relationships.
- **Color gradient along branches**: Each color represents a different temperature preference, from cool (e.g., blue) to warm (e.g., red).
- **Legend**: Shows the mapping from color to temperature bin (or value).
- **Interpretation**: This plot makes it easy to see gradual changes in temperature preference along the tree, and to spot major shifts or trends in evolution.

