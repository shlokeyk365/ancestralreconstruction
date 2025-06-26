# Load required libraries
library(geodata)
library(mvMORPH)
library(corpcor)
library(plotrix)

# 1. Load and prepare WorldClim data (bio_1 = temp, bio_12 = precip)
wclim <- worldclim_global(var = "bio", res = 10, path = "Data")
names(wclim) <- gsub("wc2.1_10m_", "", names(wclim))
wclim_df <- as.data.frame(wclim, xy = TRUE)
# Remove rows with NA values in bio_1 or bio_12
wclim_df <- wclim_df[!is.na(wclim_df$bio_1) & !is.na(wclim_df$bio_12), ]
# Scale variables
wclim_df$bio_1_scaled <- scale(wclim_df$bio_1)
wclim_df$bio_12_scaled <- scale(wclim_df$bio_12)

# 2. Simulate a phylogeny and species traits (as in Claudiatools.Rmd)
library(geiger)
tree <- sim.bdtree(b=0.1, d=0.05, stop="taxa", n=10, seed=1)
tree <- drop.extinct(tree)
tree$tip.label <- paste("sp", 1:10, sep="")

# Trait simulation helpers (from Claudiatools.Rmd)
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

# Simulate realistic traits for 2 variables (bio_1, bio_12)
simulate_realistic_traits <- function(tree, n_env_vars = 2) {
  trait_correlations <- matrix(c(1, -0.3, 0.1, -0.3, 1, -0.2, 0.1, -0.2, 1), nrow=3, byrow=TRUE)
  evolutionary_rates <- list(
    bio_1 = c(0.8, 0.15, 0.05),
    bio_12 = c(0.6, 0.12, 0.04)
  )
  root_values <- list(
    bio_1 = c(15, 8, 0.85),
    bio_12 = c(1000, 500, 0.8)
  )
  all_traits <- list()
  for(i in 1:n_env_vars) {
    env_name <- names(root_values)[i]
    A <- root_values[[env_name]]
    R_sd <- evolutionary_rates[[env_name]]
    bA <- forwardTransform1(A)
    Sigma <- corpcor::rebuild.cov(trait_correlations, R_sd^2)
    espace_traits <- mvMORPH::mvSIM(tree, nsim = 1, model = "BM1", param = list(theta = bA, sigma = Sigma, ntraits = 3))
    traits <- t(apply(espace_traits, 1, backTransform1))
    colnames(traits) <- paste(env_name, c("center", "width", "height"), sep = "_")
    all_traits[[env_name]] <- traits
  }
  combined_traits <- do.call(cbind, all_traits)
  return(combined_traits)
}
species_traits <- simulate_realistic_traits(tree)

# 3. Calculate occurrence probabilities for each species
calculate_occurrence_probabilities <- function(env_data, species_traits, species_id = 1) {
  env_vars <- c("bio_1_scaled", "bio_12_scaled")
  probabilities <- matrix(0, nrow = nrow(env_data), ncol = 2)
  for(i in 1:2) {
    env_name <- c("bio_1", "bio_12")[i]
    trait_cols <- paste(env_name, c("center", "width", "height"), sep = "_")
    traits <- species_traits[species_id, trait_cols]
    betas <- traits2coefs_sp(traits)
    env_values <- env_data[[env_vars[i]]]
    response <- betas$b0 + betas$b1 * env_values + betas$b2 * env_values^2
    probabilities[, i] <- 1 / (1 + exp(-response))
  }
  # Combine responses (additive, as in Claudia)
  prob_final <- rowMeans(probabilities)
  return(prob_final)
}

# Calculate for all species
occurrence_probs <- matrix(0, nrow = nrow(wclim_df), ncol = 10)
for(i in 1:10) {
  occurrence_probs[, i] <- calculate_occurrence_probabilities(wclim_df, species_traits, species_id = i)
}
colnames(occurrence_probs) <- paste0("species_", 1:10)

# 4. Visualize occurrence for a few example species (gridded, grayscale, no axes)
par(mfrow = c(2, 3), mar = c(0, 0, 0, 0))
for (i in 1:6) {
  plot(wclim_df$x, wclim_df$y, pch = 22, cex = 0.5,
       bg = gray.colors(100, start=0, end=1)[as.numeric(cut(occurrence_probs[, i], breaks=100))],
       col = gray.colors(100, start=0, end=1)[as.numeric(cut(occurrence_probs[, i], breaks=100))],
       axes = FALSE, xlab = "", ylab = "", main = "")
} 