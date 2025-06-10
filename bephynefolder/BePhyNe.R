# load required packages
library(devtools)
library(BePhyNE)
library(treeplyr)
library(readr)
library(tibble)
library(MCMCpack)
library(coda)
library(ape)
library(truncdist)
library(geiger)
library(phytools)
library(tidyr)
library(ratematrix)
library(mvMORPH)
library(dplyr)
library(mvtnorm)
library(MultiRNG)
library(Rphylopars)
library(Rcpp)
library(doParallel)
library(corpcor)
library(Matrix)
library(crayon)
library(shape)
library(scales)
library(robustbase)

# load example data
data(ENA_Pleth_PA)
data(ENA_Pleth_Tree)

# convert numeric columns
ENA_Pleth_PA[,c(1:2,4:6)] <- do.call(cbind, lapply(X = ENA_Pleth_PA[,c(1:2,4:6)], FUN=as.numeric))

tree <- ENA_Pleth_Tree

# prepare species data
species <- as.list(unique(ENA_Pleth_PA$species))
names(species) <- species

# format data for BePhyNE
data_final <- lapply(split(ENA_Pleth_PA[3:6], ENA_Pleth_PA$species), as.list)
for (sp in 1:length(data_final)) {
  data_final[[sp]]$species = data_final[[sp]]$species[[1]]
}

# match data order to tree tip order
data_final <- data_final[order(match(unlist(lapply(data_final, function(sp) sp$species)), tree$tip.label))]

# set up priors
center_fixed <- FALSE
v <- c(.05, .05)
k <- lapply(1:length(v), function(x) log(v[[x]]/(1-v[[x]])))

# breadth range parameters
max_bd <- 1.5
min_bd <- .1
bd <- ((log(max_bd)-log(min_bd))/4)^2

# number of environmental predictors
pred <- length(data_final$Stereochilus_marginatus[3:length(data_final[[1]])])

# prior distributions for mvBM rates
par.sd.lnorm.meanlog <- list(
  c(log(.2), log(bd)),
  c(log(.2), log(bd))
)

par.sd.lnorm.sdlog <- list(
  c(.1, .1),
  c(.1, .1)
)

par.sd.lnorm <- lapply(1:pred, function(x) cbind(par.sd.lnorm.meanlog[[x]], par.sd.lnorm.sdlog[[x]]))

# prior distributions for mvBM root means
par.mu.norm <- list(
  matrix(c(0, .5,
           log(0.3), .2),
         nrow = 2, ncol = 2, byrow = TRUE),
  matrix(c(0, .5,
           log(0.3), .2),
         nrow = 2, ncol = 2, byrow = TRUE)
)

# get GLM estimates for tolerance
GLM_only_ml <- suppressWarnings(MLglmStartpars(species_data = data_final, tree = tree, height = NULL))
heights_glm <- lapply(GLM_only_ml$start_pars_bt, function(pred) pred[,3])

# create prior functions
Prior_scale <- lapply(1:pred, function(x) makePrior_ENE(
  r = 2, 
  p = 1,
  den.mu = "norm",
  heights_by_sp = heights_glm[[x]],
  par.mu = par.mu.norm[[x]],
  den.sd = "lnorm",
  par.sd = par.sd.lnorm[[x]],
  plot = FALSE
))

# simulate starting parameters
repeat {
  startPars_scaled <- priorSim_pars(Prior_scale, tree, dist="norm", hard_coded_heights = heights_glm)
  if(length(findBadStart(res=startPars_scaled$sim_dat$sim_dat_bt, pa_data=data_final, plot=FALSE))==0) {
    break
  }
}

# tuning parameters
center_slide <- .18
center_mult <- .12
width_slide <- .15
width_mult <- .2
height_slide <- .5
height_mult <- .3

tuning <- list(
  niche_prop = lapply(1:pred, function(pred) list(
    slide = tibble(
      center = sample(center_slide, length(tree$tip.label), replace=TRUE),
      width = sample(width_slide, length(tree$tip.label), replace=TRUE),
      height = sample(height_slide, length(tree$tip.label), replace=TRUE)
    ),
    multi = tibble(
      center = sample(center_mult, length(tree$tip.label), replace=TRUE),
      width = sample(width_mult, length(tree$tip.label), replace=TRUE),
      height = sample(height_mult, length(tree$tip.label), replace=TRUE)
    )
  )),
  w_mu = lapply(1:pred, function(pred) list(
    slide = c(.6,.6),
    multi = c(.6,.6)
  )),
  w_sd = lapply(1:pred, function(pred) list(
    slide = c(.15,.12),
    multi = c(.15,.12)
  )),
  v_cor = lapply(1:pred, function(pred) 100)
)

# mcmc parameters
sets_full <- suppressWarnings(separate.data(data_final, ratio = .5))
sparse_sp <- FALSE
iterations <- 1000
trim_freq <- 10
burnin <- 0
chain_end <- (iterations-(burnin))/trim_freq
plot <- FALSE
plot_freq <- iterations/5

# move weights
moves_wieghts <- c(
  "height" = 2,
  "center" = 3,
  "width" = 3,
  "theta" = 1,
  "R_corr" = 1,
  "R_sd" = 1
)

move_prob <- c(
  "height" = moves_wieghts[[1]]/sum(moves_wieghts),
  "center" = moves_wieghts[[2]]/sum(moves_wieghts),
  "width" = moves_wieghts[[3]]/sum(moves_wieghts),
  "theta" = moves_wieghts[[4]]/sum(moves_wieghts),
  "R_corr" = moves_wieghts[[5]]/sum(moves_wieghts),
  "R_sd" = moves_wieghts[[6]]/sum(moves_wieghts)
)

# run mcmc
results <- metro_haste_full_MV(
  R_corr_start = startPars_scaled$R$R_cor,
  R_sd_start = startPars_scaled$R$R_sd,
  A_start = startPars_scaled$A$A_bt,
  Prior = Prior_scale,
  tree = tree,
  tibble_data = startPars_scaled$sim_dat$sim_td_bt,
  pa_data = sets_full$training,
  iterations = iterations,
  burnin = burnin,
  move_prob = move_prob,
  n = 2,
  print.i.freq = 1000,
  print.ac.freq = 100,
  printing = TRUE,
  trim = TRUE,
  trim_freq = trim_freq,
  H_fixed = FALSE,
  tuning = tuning,
  center_fixed = center_fixed,
  write_file = FALSE,
  IDlen = 5,
  dir = NA,
  outname = NA,
  prior_only = FALSE,
  glm_only = FALSE,
  plot = FALSE,
  plot_freq = iterations/5,
  plot_file = NA,
  True_pars = NA,
  k = k
)

# convert results to mcmc objects
mcmc_lt <- BePhyNE_out2coda_mcmc(results, tree)

# calculate hpd intervals
HPDs <- HPD_list(mcmc_lt, tree, prob=0.95)

# plot mcmc chains
pdf("mcmc_chains.pdf")
par(mfrow=c(2,2))

plot(mcmc_lt$A$pred_1$optimum, type="l", main="Optimum Chain (pred_1)", ylab="Optimum", xlab="Iteration")
plot(mcmc_lt$A$pred_1$breadth, type="l", main="Breadth Chain (pred_1)", ylab="Breadth", xlab="Iteration")

# plot pred_2 if available
if (!is.null(mcmc_lt$A$pred_2$optimum)) {
  plot(mcmc_lt$A$pred_2$optimum, type="l", main="Optimum Chain (pred_2)", ylab="Optimum", xlab="Iteration")
} else {
  plot.new(); title(main="No pred_2 optimum")
}
if (!is.null(mcmc_lt$A$pred_2$breadth)) {
  plot(mcmc_lt$A$pred_2$breadth, type="l", main="Breadth Chain (pred_2)", ylab="Breadth", xlab="Iteration")
} else {
  plot.new(); title(main="No pred_2 breadth")
}

dev.off()

# plot ancestral reconstruction
pdf("ancestral_reconstruction.pdf")
plot(tree, type="fan", show.tip.label=TRUE, tip.color="black", label.offset=0.1)
nodelabels(pie=HPDs$A$pred_1$optimum, piecol=c("red", "blue"), cex=0.5)
legend("bottomright", legend=c("State 1", "State 2"), fill=c("red", "blue"), bty="n")
dev.off()

# continuous trait ancestral state plot with color gradient
library(phytools)

# use tip means for contMap
tip_means <- apply(mcmc_lt$tip_curves$pred_1$optimum, 1, mean)
names(tip_means) <- tree$tip.label

# ensure tip_means matches tree tip labels
tip_means <- tip_means[tree$tip.label]

# check for missing matches
if (any(is.na(tip_means))) {
  cat("warning: some tip means are NA. check for mismatched tip labels.\n")
  print(tree$tip.label[is.na(tip_means)])
}

# plot if all tips match
if (!any(is.na(tip_means))) {
  pdf("ancestral_reconstruction_gradient.pdf")
  contmap_obj <- contMap(tree, tip_means, plot=FALSE)
  plot(contmap_obj, legend=0.7*max(nodeHeights(tree)), fsize=c(0.5,0.7), outline=FALSE, lwd=3)
  title("Ancestral State Reconstruction: Optimum (Color Gradient)")
  dev.off()
} else {
  cat("error: could not plot because of mismatched tip labels.\n")
}
