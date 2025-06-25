
library(mvMORPH)
library(geiger)
library(plotrix)
library(raster)
#library(rgdal) -> not necessary for the simulation to run (rgdal doesn't exist anymore in CRAN)


traits2coefs <- function(traits, v=0.05){
  if (ncol(as.matrix(traits))==1){
    theta=traits[[1]]
    w=traits[[2]]
    height=traits[[3]]
  } else{
    theta=traits[,1]
    w=traits[,2]
    height=traits[,3]
  }
  
  H=log((1/height)-1)
  k=log(v/(1-v))
  b0<- ((theta^2)/(w^2))*k + H*(((theta^2)/(w^2))-1)
  b1<- (-2*(H+k)*theta)/w^2
  b2<- (1/w^2)*(H+k)
  return(data.frame(b0=b0,b1=b1, b2=b2))
}

traits2coefs_sp<-function(traits, v=0.05){
  theta=traits[1]
  w=traits[2]
  height=traits[3]
  H=log((1/height)-1)
  k=log(v/(1-v))
  b0<- ((theta^2)/(w^2))*k + H*(((theta^2)/(w^2))-1)
  b1<- (-2*(H+k)*theta)/w^2
  b2<- (1/w^2)*(H+k)
  return(data.frame(b0=b0,b1=b1, b2=b2))
}

coef2traits<-function(coefs, v=0.05){
  B0=coefs[[1]]
  B1=coefs[[2]]
  B2=coefs[[3]]
  
  theta=B1*-1/(2*B2)
  height=1 / (1 + exp( B1^2 / (4*B2) - B0) )
  width = -sqrt( B1^2 - 4*B2 * (B0 - log( (v)/(1-v) ) ) ) / (2*B2)
  
  return(c(theta, width, height))
}

forwardTransform1 <- function(x){
  newx <- x
  newx[2] <- log(x[2])
  Y=((x[3]-.05)/.95)
  newx[3] <- -1*log(Y/(1-Y))
  return(newx)
}

backTransform1 <- function(x){
  newx <- x
  newx[2] <- exp(x[2])
  newx[3] <- 0.05 + 0.95*(exp(-1*x[3])/(1+exp(-1*x[3])))
  return(newx)
}

sim_response_curve <- function(tree, A, corSig, R_sd){
  bA <- forwardTransform1(A)
  
  if(!all(eigen(corSig)$values >0)) stop("Correlation matrix not positive-definite")
  Sigma <- corpcor::rebuild.cov(corSig, R_sd^2)
  espace_traits <- mvMORPH::mvSIM(tree, nsim=1, model="BM1", param=list(theta=bA, sigma=Sigma, ntraits=3))
  traits <- t(apply(espace_traits, 1, backTransform1))
  colnames(traits) <- c("center", "width", "height")
  return(traits)
}






create_realistic_environment <- function(n_lat = 100, n_long = 200, seed = 123) {
  set.seed(seed)
  
  # Geographic coordinates
  lat <- seq(-60, 60, length.out = n_lat)
  long <- seq(-180, 180, length.out = n_long)
  coords <- expand.grid(lat = lat, long = long)
  
  # Environmental variables with realistic patterns
  coords$temperature <- 25 - 0.6 * abs(coords$lat) + 
    rnorm(nrow(coords), 0, 2) + 
    5 * sin(coords$long * pi / 180) # Temperature gradient with longitude effect
  
  coords$precipitation <- 1000 + 50 * coords$lat + 
    rnorm(nrow(coords), 0, 100) + 
    200 * cos(coords$long * pi / 90) # Precipitation pattern
  
  coords$elevation <- 500 + 20 * abs(coords$lat) + 
    rnorm(nrow(coords), 0, 200) + 
    100 * sin(coords$long * pi / 45) # Elevation gradient
  
  coords$soil_ph <- 6.5 + 0.02 * coords$lat + 
    rnorm(nrow(coords), 0, 0.3) + 
    0.5 * cos(coords$long * pi / 60) # Soil pH
  
  coords$vegetation_cover <- 0.8 - 0.005 * abs(coords$lat) + 
    rnorm(nrow(coords), 0, 0.1) + 
    0.2 * sin(coords$long * pi / 120) # Vegetation cover
  
  # Ensure realistic bounds
  coords$temperature <- pmax(pmin(coords$temperature, 40), -20)
  coords$precipitation <- pmax(coords$precipitation, 0)
  coords$elevation <- pmax(coords$elevation, 0)
  coords$soil_ph <- pmax(pmin(coords$soil_ph, 9), 4)
  coords$vegetation_cover <- pmax(pmin(coords$vegetation_cover, 1), 0)
  
  return(coords)
}

#Simulate realistic species traits with phylogenetic constraints

simulate_realistic_traits <- function(tree, n_env_vars = 5, 
                                      trait_correlations = NULL,
                                      evolutionary_rates = NULL,
                                      seed = 456) {
  set.seed(seed)
  
  # Default parameters if not provided
  if(is.null(trait_correlations)) {
    # Realistic correlation structure between center, width, and height
    trait_correlations <- matrix(c(1, -0.3, 0.1,    # center-width-height correlations
                                   -0.3, 1, -0.2,
                                   0.1, -0.2, 1), nrow = 3, byrow = TRUE)
  }
  
  if(is.null(evolutionary_rates)) {
    # Different rates for different environmental variables
    evolutionary_rates <- list(
      temperature = c(0.8, 0.15, 0.05),    # center, width, height rates
      precipitation = c(0.6, 0.12, 0.04),
      elevation = c(0.7, 0.14, 0.05),
      soil_ph = c(0.4, 0.08, 0.03),
      vegetation_cover = c(0.5, 0.10, 0.04)
    )
  }
  
  root_values <- list(
    temperature = c(15, 8, 0.85),      # moderate temperature preference
    precipitation = c(0, 2, 0.80),     # moderate precipitation
    elevation = c(0, 1.5, 0.75),       # low elevation preference
    soil_ph = c(0, 1, 0.70),           # neutral pH preference
    vegetation_cover = c(0, 1.2, 0.80) # moderate vegetation
  )
  
  all_traits <- list()
  
  for(i in 1:n_env_vars) {
    env_name <- names(root_values)[i]
    A <- root_values[[env_name]]
    R_sd <- evolutionary_rates[[env_name]]
    
    # Transform to evolutionary space
    bA <- forwardTransform1(A)
    
    # Simulate traits
    Sigma <- corpcor::rebuild.cov(trait_correlations, R_sd^2)
    espace_traits <- mvMORPH::mvSIM(tree, nsim = 1, model = "BM1", 
                                    param = list(theta = bA, sigma = Sigma, ntraits = 3))
    
    # Back-transform to niche space
    traits <- t(apply(espace_traits, 1, backTransform1))
    colnames(traits) <- paste(env_name, c("center", "width", "height"), sep = "_")
    
    all_traits[[env_name]] <- traits
  }
  
  # Combine all traits
  combined_traits <- do.call(cbind, all_traits)
  return(combined_traits)
}


calculate_occurrence_probabilities <- function(env_data, species_traits, 
                                               species_id = 1,
                                               interaction_strength = 0.1,
                                               noise_level = 0.05) {
  
  # Scale environmental variables
  env_scaled <- env_data
  env_vars <- c("temperature", "precipitation", "elevation", "soil_ph", "vegetation_cover")
  env_scaled[env_vars] <- scale(env_data[env_vars])
  
  # Extract traits for this species
  n_env <- length(env_vars)
  probabilities <- matrix(0, nrow = nrow(env_data), ncol = n_env)
  
  for(i in 1:n_env) {
    env_name <- env_vars[i]
    trait_cols <- paste(env_name, c("center", "width", "height"), sep = "_")
    traits <- species_traits[species_id, trait_cols]
    
    # Convert to logistic coefficients
    betas <- traits2coefs_sp(traits)
    
    # Calculate response for this environmental variable
    env_values <- env_scaled[[env_name]]
    response <- betas$b0 + betas$b1 * env_values + betas$b2 * env_values^2
    probabilities[, i] <- 1 / (1 + exp(-response))
  }
  
  # Combine responses across environmental variables
  # Option 1: Multiplicative (all conditions must be suitable)
  prob_mult <- apply(probabilities, 1, prod)
  
  # Option 2: Additive with interactions
  prob_add <- rowMeans(probabilities)
  
  # Option 3: Weighted combination (most realistic)
  weights <- c(0.3, 0.25, 0.2, 0.15, 0.1) # Temperature most important
  prob_weighted <- as.vector(probabilities %*% weights)
  
  # Add environmental interactions
  interaction_terms <- interaction_strength * 
    (probabilities[, 1] * probabilities[, 2] + 
       probabilities[, 2] * probabilities[, 3] + 
       probabilities[, 3] * probabilities[, 4])
  
  prob_final <- prob_weighted + interaction_terms
  
  # Add stochastic noise
  prob_final <- prob_final + rnorm(length(prob_final), 0, noise_level)
  
  # Ensure probabilities are between 0 and 1
  prob_final <- pmax(pmin(prob_final, 1), 0)
  
  return(prob_final)
}


simulate_species_occurrences <- function(env_data, species_traits, n_species = 20) {
  occurrence_probs <- matrix(0, nrow = nrow(env_data), ncol = n_species)
  
  for(i in 1:n_species) {
    occurrence_probs[, i] <- calculate_occurrence_probabilities(
      env_data, species_traits, species_id = i
    )
  }
  
  colnames(occurrence_probs) <- paste("species", 1:n_species, sep = "_")
  return(occurrence_probs)
}


simulate_community_assembly <- function(occurrence_probs,
                                        env_data,
                                        competition_strength = 0.3,
                                        dispersal_limitation = 0.2) {
  
  n_sites <- nrow(occurrence_probs)
  n_species <- ncol(occurrence_probs)
  
  # Initialize community matrix
  community <- matrix(0, nrow = n_sites, ncol = n_species)
  
  # Random colonization based on occurrence probabilities
  for(i in 1:n_sites) {
    for(j in 1:n_species) {
      if(runif(1) < occurrence_probs[i, j]) {
        community[i, j] <- 1
      }
    }
  }
  
  for(i in 1:n_sites) {
    present_species <- which(community[i, ] == 1)
    if(length(present_species) > 1) {
      # Calculate competitive effects, handle NAs
      competitive_pressure <- rowSums(occurrence_probs[i, present_species, drop = FALSE], na.rm = TRUE)
      
      # Remove some species based on competition
      for(idx in seq_along(present_species)) {
        j <- present_species[idx]
        cp <- competitive_pressure[idx]
        if(!is.na(cp) && runif(1) < competition_strength * cp) {
          community[i, j] <- 0
        }
      }
    }
  }
  
  # Apply dispersal limitation (reduce colonization in isolated areas)
  for(i in 1:n_sites) {
    for(j in 1:n_species) {
      if(community[i, j] == 0 && occurrence_probs[i, j] > 0.5) {
        # Check if nearby sites have this species
        nearby_sites <- which(abs(env_data$lat - env_data$lat[i]) < 5 & 
                                abs(env_data$long - env_data$long[i]) < 5)
        nearby_occurrence <- any(community[nearby_sites, j] == 1)
        
        if(!nearby_occurrence && runif(1) < dispersal_limitation) {
          community[i, j] <- 0
        }
      }
    }
  }
  
  return(community)
}

#' Load real environmental data from WorldClim

#' @param data_path Path to WorldClim data
#' @return Environmental data frame
load_real_environmental_data <- function(data_path = "Data/climate/wc2.1_10m/") {
  # List available climate variables
  climate_files <- list.files(data_path, pattern = "\\.tif$", full.names = TRUE)
  
  # Load key climate variables
  temp_mean <- raster(file.path(data_path, "wc2.1_10m_bio_1.tif")) # Annual mean temperature
  temp_range <- raster(file.path(data_path, "wc2.1_10m_bio_7.tif")) # Temperature annual range
  precip <- raster(file.path(data_path, "wc2.1_10m_bio_12.tif"))    # Annual precipitation
  
  # Create a sample region (e.g., North America)
  extent_na <- extent(-140, -50, 25, 60)
  temp_mean_na <- crop(temp_mean, extent_na)
  temp_range_na <- crop(temp_range, extent_na)
  precip_na <- crop(precip, extent_na)
  
  # Convert to data frame
  coords <- coordinates(temp_mean_na)
  env_real <- data.frame(
    long = coords[, 1],
    lat = coords[, 2],
    temperature = values(temp_mean_na),
    temp_range = values(temp_range_na),
    precipitation = values(precip_na)
  )
  
  # Remove NA values
  env_real <- na.omit(env_real)
  
  return(env_real)
}

# the main simulation function
run_realistic_simulation <- function(n_species = 20, 
                                     n_lat = 100, 
                                     n_long = 200,
                                     use_real_data = FALSE,
                                     seed = 123) {
  
  set.seed(seed)
  
  # 1. Generate phylogeny
  cat("Generating phylogeny...\n")
  tree <- sim.bdtree(b=0.1, d=0.05, stop="taxa", n=n_species, seed=seed)
  tree <- drop.extinct(tree)
  tree$tip.label <- paste("sp", 1:n_species, sep="")
  tree$edge.length <- tree$edge.length/max(branching.times(tree))*100
  
  # 2. Generate environmental data
  cat("Generating environmental data...\n")
  if(use_real_data) {
    tryCatch({
      env_data <- load_real_environmental_data()
      cat("Using real WorldClim data\n")
    }, error = function(e) {
      cat("Real data not available, using simulated data\n")
      env_data <- create_realistic_environment(n_lat, n_long, seed)
    })
  } else {
    env_data <- create_realistic_environment(n_lat, n_long, seed)
  }
  
  # 3. Simulate species traits
  cat("Simulating species traits...\n")
  species_traits <- simulate_realistic_traits(tree, seed = seed + 1)
  
  # 4. Calculate occurrence probabilities
  cat("Calculating occurrence probabilities...\n")
  occurrence_probs <- simulate_species_occurrences(env_data, species_traits, n_species)
  
  # 5. Simulate community assembly
  cat("Simulating community assembly...\n")
  community_matrix <- simulate_community_assembly(occurrence_probs, env_data)
  
  # 6. Calculate summary statistics
  cat("Calculating summary statistics...\n")
  species_richness <- rowSums(occurrence_probs > 0.5)
  community_richness <- rowSums(community_matrix)
  community_diversity <- apply(community_matrix, 1, function(x) {
    if(sum(x) == 0) return(0) else return(-sum(x/sum(x) * log(x/sum(x))))
  })
  
  # Return results
  results <- list(
    tree = tree,
    env_data = env_data,
    species_traits = species_traits,
    occurrence_probs = occurrence_probs,
    community_matrix = community_matrix,
    species_richness = species_richness,
    community_richness = community_richness,
    community_diversity = community_diversity
  )
  
  cat("Simulation complete!\n")
  return(results)
}


plot_environmental_gradients <- function(env_data) {
  par(mfrow = c(2, 3), mar = c(2, 2, 2, 1))
  for(var in c("temperature", "precipitation", "elevation", "soil_ph", "vegetation_cover")) {
    plot(env_data$long, env_data$lat, pch = 15, cex = 0.5,
         col = color.scale(env_data[[var]], color.spec = "rgb"),
         main = var, xlab = "", ylab = "")
  }
}


plot_species_distributions <- function(env_data, occurrence_probs, n_species = 20) {
  par(mfrow = c(4, 5), mar = c(1, 1, 1, 1))
  for(i in 1:n_species) {
    plot(env_data$long, env_data$lat, pch = 15, cex = 0.3,
         col = color.scale(occurrence_probs[, i], color.spec = "rgb"),
         main = paste("Species", i), xlab = "", ylab = "")
  }
}


plot_community_patterns <- function(env_data, community_richness, community_diversity) {
  par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
  
  # Richness map
  plot(env_data$long, env_data$lat, pch = 15, cex = 0.5,
       col = color.scale(community_richness, color.spec = "rgb"),
       main = "Community Richness", xlab = "Longitude", ylab = "Latitude")
  
  # Diversity map
  plot(env_data$long, env_data$lat, pch = 15, cex = 0.5,
       col = color.scale(community_diversity, color.spec = "rgb"),
       main = "Community Diversity", xlab = "Longitude", ylab = "Latitude")
  
  # Richness vs latitude
  plot(env_data$lat, community_richness, pch = 16, cex = 0.5,
       main = "Latitudinal Diversity Gradient",
       xlab = "Latitude", ylab = "Species Richness")
  abline(lm(community_richness ~ env_data$lat), col = "red", lwd = 2)
}
