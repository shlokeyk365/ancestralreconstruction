# Project Progress Summary

## Overview
This document summarizes the progress and work completed on the entire project, including the main project and all its subprojects. The project focuses on ancestral state reconstruction for environmental variables using various methods, including a binned approach on a phylogenetic tree.

## Main Project: Ancestral Reconstruction
- **Objective**: Perform ancestral state reconstruction for environmental variables using multiple methods.
- **Progress**:
  - Successfully downloaded and prepared the precipitation data using WorldClim data.
  - Updated the script to use the new precipitation data file.
  - Generated visualizations for both temperature and precipitation ancestral reconstructions.
  - Created a comprehensive README file for the project.

## Subproject 1: Download Climate Data
- **Objective**: Download and prepare the precipitation data for use in the ancestral reconstruction.
- **Progress**:
  - Created a script (`download_climate_data.R`) to download WorldClim precipitation data.
  - Successfully generated the `precip.tif` file and a visualization (`precipitation_map.png`).

## Subproject 2: Update Ancestral Reconstruction Script
- **Objective**: Update the ancestral reconstruction script to use the new precipitation data.
- **Progress**:
  - Modified the script to use the new precipitation file path.
  - Fixed errors related to the `color.legend` function by adding the required `rect.col` argument.
  - Successfully ran the updated script, generating all intended outputs.

## Subproject 3: Create README Documentation
- **Objective**: Create comprehensive documentation for the project.
- **Progress**:
  - Created a README file for the multi-dimensional binned ancestral reconstruction program, detailing its purpose, usage, dependencies, and outputs.
  - Summarized the overall project progress in this document.

## Subproject 4: Birth-Death Models
- **Objective**: Implement and analyze birth-death models for the project.
- **Progress**:
  - Created a directory for birth-death models.
  - Implemented models and saved results in `environmental_birthdeath_results.RData`.

## Subproject 5: Phylogenetic Tree Analysis
- **Objective**: Analyze the phylogenetic tree for the project.
- **Progress**:
  - Utilized the `etheostoma_percina_chrono.tre` file for analysis.

## Subproject 6: Nichevol Package Integration
- **Objective**: Integrate and utilize the nichevol package for niche modeling and analysis.
- **Progress**:
  - Utilized the nichevol package for niche modeling and analysis.
  - Generated niche models and integrated them into the ancestral reconstruction process.

## Differences Between Methods
- **Nichevol**: Focuses on niche modeling and analysis, integrating environmental data to understand species' ecological niches and their evolution over time. The nichevol package uses a multivariate normal distribution in its niche modeling approach. Specifically, it uses this distribution to model the environmental niche space of species, where each species' niche is characterized by a mean vector (representing optimal conditions) and a covariance matrix (representing the breadth and shape of the niche). This multivariate normal approach allows the package to capture the complex relationships between multiple environmental variables and their influence on species distributions.

- **Bephyne**: Likely involves a different approach to niche modeling or phylogenetic analysis, possibly focusing on specific aspects of species evolution or environmental adaptation.
- **Birth-Death Models**: These models are used to simulate and analyze the evolutionary dynamics of species, focusing on speciation and extinction rates over time.

## Next Steps
- Review the generated visualizations and outputs.
- Address any warnings or errors in the scripts.
- Consider extending the project to include additional environmental variables or enhance the visualizations.
- Further analyze birth-death models and phylogenetic tree data.
- Explore additional features and capabilities of the nichevol package.

## Conclusion
The project has made significant progress in downloading climate data, updating the ancestral reconstruction script, and creating comprehensive documentation. The next steps involve reviewing the outputs and addressing any remaining issues. 
