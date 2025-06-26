# Enhanced Claudia Workflow - Realistic Species Distribution Simulation
# Follows the same process as Claudiatools.Rmd but with realistic improvements
# Author: Enhanced simulation framework
# Date: 2025-01-20

# Required libraries
library(mvMORPH)
library(geiger)
library(plotrix)

# =============================================================================
# ORIGINAL FUNCTIONS (from Claudia project) - UNCHANGED
# =============================================================================

# Function converts center, width and height to logistic regression coefficients
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

# Converts values for a vector from an individual species
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

# Converts logistic regression coefficients to center, width, and height
coef2traits<-function(coefs, v=0.05){
  B0=coefs[[1]]
  B1=coefs[[2]]
  B2=coefs[[3]]
  
  theta=B1*-1/(2*B2)
  height=1 / (1 + exp( B1^2 / (4*B2) - B0) )
  width = -sqrt( B1^2 - 4*B2 * (B0 - log( (v)/(1-v) ) ) ) / (2*B2)
  
  return(c(theta, width, height))
}

# Evolutionary space transformation
forwardTransform1 <- function(x){
  newx <- x
  newx[2] <- log(x[2])
  Y=((x[3]-.05)/.95)
  newx[3] <- -1*log(Y/(1-Y))
  return(newx)
}

# Back transform from BM/evolutionary space to niche space
backTransform1 <- function(x){
  newx <- x
  newx[2] <- exp(x[2])
  newx[3] <- 0.05 + 0.95*(exp(-1*x[3])/(1+exp(-1*x[3])))
  return(newx)
}

# Original simulation function
sim_response_curve <- function(tree, A, corSig, R_sd){
  bA <- forwardTransform1(A)
  
  if(!all(eigen(corSig)$values >0)) stop("Correlation matrix not positive-definite")
  Sigma <- corpcor::rebuild.cov(corSig, R_sd^2)
  espace_traits <- mvMORPH::mvSIM(tree, nsim=1, model="BM1", param=list(theta=bA, sigma=Sigma, ntraits=3))
  traits <- t(apply(espace_traits, 1, backTransform1))
  colnames(traits) <- c("center", "width", "height")
  return(traits)
}

# =============================================================================
# ENHANCED WORKFLOW - FOLLOWING CLAUDIA PROCESS BUT MORE REALISTIC
# =============================================================================

## Step 1: Simulate phylogeny (same as Claudia)
cat("Step 1: Simulating phylogeny...\n")
tree <- sim.bdtree(b=0.1, d=0.05, stop="taxa", n=100, seed=1)
tree <- drop.extinct(tree)
tree$tip.label <- paste("sp", 1:100, sep="")
tree$edge.length <- tree$edge.length/max(branching.times(tree))*100 #Rescale tree to 100 my
plot(tree, main="Phylogenetic Tree")

## Step 2: Enhanced trait simulation (improved from Claudia)
cat("Step 2: Simulating realistic traits...\n")

# ENHANCEMENT 1: More realistic correlation structure
# Original Claudia: Simple correlation matrix
# Enhanced: Realistic correlations based on ecological theory
corSig <- matrix(c(1, -0.4, 0.1,    # center-width-height correlations
                   -0.4, 1, -0.2,    # specialists vs generalists trade-off
                   0.1, -0.2, 1), nrow=3, ncol=3, byrow = TRUE)

# ENHANCEMENT 2: More realistic root values and evolutionary rates
# Set root trait centers to the mean of each environmental variable
A_temperature <- c(env_means["temperature"], 8, 0.85)      # center, width, height
A_precipitation <- c(env_means["precipitation"], 2, 0.80)
A_elevation <- c(env_means["elevation"], 1.5, 0.75)
A_soil_ph <- c(env_means["soil_ph"], 1, 0.70)
A_vegetation <- c(env_means["vegetation_cover"], 1.2, 0.80)

# ENHANCEMENT 3: Lower evolutionary rates for trait centers to reduce variance
R_sd_temperature <- c(0.3, 0.15, 0.05)    # center, width, height
R_sd_precipitation <- c(0.2, 0.12, 0.04)
R_sd_elevation <- c(0.2, 0.14, 0.05)
R_sd_soil_ph <- c(0.15, 0.08, 0.03)
R_sd_vegetation <- c(0.2, 0.10, 0.04)

# Simulate traits for each environmental variable
env_vars <- c("temperature", "precipitation", "elevation", "soil_ph", "vegetation_cover")
A_list <- list(A_temperature, A_precipitation, A_elevation, A_soil_ph, A_vegetation)
R_list <- list(R_sd_temperature, R_sd_precipitation, R_sd_elevation, R_sd_soil_ph, R_sd_vegetation)

Traits <- NULL
for(i in 1:length(env_vars)){
  A <- A_list[[i]]
  R_sd <- R_list[[i]]
  bA <- forwardTransform1(A)
  
  Sigma <- corpcor::rebuild.cov(corSig, R_sd^2)
  espace_traits <- mvMORPH::mvSIM(tree, nsim=1, model="BM1", 
                                  param=list(theta=bA, sigma=Sigma, ntraits=3))
  traits <- t(apply(espace_traits, 1, backTransform1))
  colnames(traits) <- paste(env_vars[i], c("center", "width", "height"), sep="_")
  Traits <- cbind(Traits, traits)
}

# ENHANCEMENT 4: Allow height to vary between environmental variables (more realistic)
# Original Claudia: Traits[,seq(3,ncol(Traits),3)] <- Traits[,3] #every height has to be the same
# Enhanced: Keep different heights for different environmental variables

head(Traits)

## Step 3: Enhanced environmental data generation (improved from Claudia)
cat("Step 3: Generating realistic environmental data...\n")

# ENHANCEMENT 5: More realistic environmental gradients
# Original Claudia: Simple linear gradients with random noise
# Enhanced: Realistic geographic patterns with spatial autocorrelation

lat <- seq(-60, 60, length.out=100)  # Increased resolution
long <- seq(-180, 180, length.out=200)  # Increased resolution
coords <- expand.grid(lat=lat, long=long)

# Realistic environmental patterns
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

# Visualize environmental gradients (same plot style as Claudia)
par(mfrow=c(2,3), mar=c(0,0,0,0))
plot(coords$long, coords$lat, pch=22, bg=color.scale(coords$temperature, color.spec="rgb"), 
     col=color.scale(coords$temperature, color.spec="rgb"), main="Temperature")
plot(coords$long, coords$lat, pch=22, bg=color.scale(coords$precipitation, color.spec="rgb"), 
     col=color.scale(coords$precipitation, color.spec="rgb"), main="Precipitation")
plot(coords$long, coords$lat, pch=22, bg=color.scale(coords$elevation, color.spec="rgb"), 
     col=color.scale(coords$elevation, color.spec="rgb"), main="Elevation")
plot(coords$long, coords$lat, pch=22, bg=color.scale(coords$soil_ph, color.spec="rgb"), 
     col=color.scale(coords$soil_ph, color.spec="rgb"), main="Soil pH")
plot(coords$long, coords$lat, pch=22, bg=color.scale(coords$vegetation_cover, color.spec="rgb"), 
     col=color.scale(coords$vegetation_cover, color.spec="rgb"), main="Vegetation Cover")

## Step 4: Scale environmental variables (same as Claudia)
cat("Step 4: Scaling environmental variables...\n")
coords_scaled <- coords
coords_scaled[,-c(1,2)] <- apply(coords_scaled[,-c(1,2)], 2, scale)

# ✅ 1. Store environment scaling factors
env_vars <- c("temperature", "precipitation", "elevation", "soil_ph", "vegetation_cover")
env_means <- apply(coords[ , env_vars], 2, mean)
env_sds <- apply(coords[ , env_vars], 2, sd)

# After Traits is created
# ✅ 2. Scale trait centers to match scaled environment
for (i in seq_along(env_vars)) {
  col_name <- paste0(env_vars[i], "_center")
  Traits[, col_name] <- (Traits[, col_name] - env_means[i]) / env_sds[i]
  # Clamp to -3:3
  Traits[, col_name] <- pmin(pmax(Traits[, col_name], -3), 3)
}

# Print diagnostic ranges for trait centers and scaled environment
cat("\nTrait center and environment diagnostics (after scaling):\n")
for (v in env_vars) {
  cat(v, "trait center (scaled) range:", range(Traits[, paste0(v, "_center")]), "\n")
  cat(v, "scaled env range:", range(coords_scaled[, v]), "\n")
}

# ✅ Clamp width values to avoid overly generalist/specialist species (min 1.0)
width_cols <- grep("_width$", colnames(Traits), value = TRUE)
Traits[, width_cols] <- pmin(pmax(Traits[, width_cols], 1.0), 2.5)

## Step 5: Enhanced occurrence probability calculation (improved from Claudia)
cat("Step 5: Calculating occurrence probabilities...\n")

# ENHANCEMENT 6: More sophisticated occurrence modeling
# Original Claudia: Simple additive combination
# Enhanced: Weighted combination with environmental interactions

# Function to calculate realistic occurrence probabilities
calculate_realistic_occurrences <- function(coords_scaled, Traits, species_id = 1) {
  K <- length(env_vars)
  traits_ind <- seq(1, ncol(Traits), 3)
  
  # Get traits for this species
  betas <- lapply(traits_ind, function(k) {
    traits <- Traits[species_id, k:(k+1)]
    height <- Traits[species_id, k+2]
    traits2coefs_sp(c(traits, height))
  })
  
  YY <- list()
  for(k in 1:K){
    D <- cbind(1, coords_scaled[,2+k], coords_scaled[,2+k]^2)
    Y1 <- as.vector(as.matrix(D) %*% t(as.matrix(betas[[k]])))
    YY[[k]] <- Y1
  }
  
  weights <- c(0.3, 0.25, 0.2, 0.15, 0.1)
  prob_weighted <- as.vector(do.call(cbind, YY) %*% weights)
  
  # Debug: Check logistic input range before applying sigmoid
  cat("Species", species_id, "- linear predictor range: ", range(prob_weighted), "\n")
  
  fixed_height <- 0.95
  # Use a gentler slope (-1) for the logistic transformation
  PP <- fixed_height * 1 / (1 + exp(-1 * prob_weighted))
  
  noise_level <- 0.01
  PP <- PP + rnorm(length(PP), 0, noise_level)
  
  PP <- pmax(pmin(PP, 1), 0)
  
  return(PP)
}

## Step 6: Plot occurrence probability maps (improved color scale)
cat("Step 6: Plotting occurrence probability maps...\n")

# Use heat.colors for better contrast
library(viridis)
par(mar=c(0,0,0,0), mfrow=c(3,3))
for(i in 1:9){
  PP <- calculate_realistic_occurrences(coords_scaled, Traits, species_id = i)
  range(PP)
  
  plot(coords$long, coords$lat, pch=22, bg=viridis::viridis(100)[as.numeric(cut(PP, breaks=100))], 
       col=viridis::viridis(100)[as.numeric(cut(PP, breaks=100))], main=paste("Species", i))
}

## Step 7: Enhanced analysis (additional to Claudia)
cat("Step 7: Analyzing patterns...\n")

# Calculate species richness patterns
species_richness <- matrix(0, nrow=nrow(coords), ncol=20)
for(i in 1:20){
  species_richness[,i] <- calculate_realistic_occurrences(coords_scaled, Traits, species_id = i)
}
total_richness <- rowSums(species_richness > 0.5)

# Plot richness patterns
par(mfrow=c(1,2), mar=c(4,4,2,1))
plot(coords$long, coords$lat, pch=22, bg=color.scale(total_richness, color.spec="rgb"), 
     col=color.scale(total_richness, color.spec="rgb"), 
     main="Species Richness", xlab="Longitude", ylab="Latitude")

# Latitudinal diversity gradient
plot(coords$lat, total_richness, pch=16, cex=0.5,
     main="Latitudinal Diversity Gradient",
     xlab="Latitude", ylab="Species Richness")
abline(lm(total_richness ~ coords$lat), col="red", lwd=2)

cat("Enhanced Claudia workflow complete!\n") 