# libraries
library(phytools)
library(ape)
library(diversitree)
library(ggplot2)
library(tidyverse)

# simulate environmental data through time
simulate_environmental_data <- function(times, n_points = 100) {
  time_seq <- seq(min(times), max(times), length.out = n_points)
  temperature <- 15 + 5 * sin(time_seq/10) + rnorm(n_points, 0, 1)
  precipitation <- 1000 + 200 * cos(time_seq/8) + rnorm(n_points, 0, 50)
  
  data.frame(
    time = time_seq,
    temperature = temperature,
    precipitation = precipitation
  )
}

# calculate rates from environmental data
calculate_rates <- function(env_data, base_speciation = 0.1, base_extinction = 0.05) {
  temp_effect <- exp(-(env_data$temperature - 20)^2 / 50)
  precip_effect <- exp(-(env_data$precipitation - 1000)^2 / 50000)
  
  speciation_rate <- base_speciation * temp_effect * precip_effect
  extinction_rate <- base_extinction * (1 - temp_effect * precip_effect)
  
  data.frame(
    time = env_data$time,
    speciation_rate = speciation_rate,
    extinction_rate = extinction_rate
  )
}

# simulate tree using average rates
simulate_environmental_tree <- function(env_data, rates, n_tips = 100) {
  mean_speciation <- mean(rates$speciation_rate)
  mean_extinction <- mean(rates$extinction_rate)
  
  pbtree(b = mean_speciation, d = mean_extinction, n = n_tips, scale = max(env_data$time), method = "direct")
}

# main

# generate env data and rates
times <- seq(0, 100, length.out = 100)
env_data <- simulate_environmental_data(times)
rates <- calculate_rates(env_data)

# simulate tree
tree <- simulate_environmental_tree(env_data, rates)

# plots
par(mfrow = c(2, 2))

plot(env_data$time, env_data$temperature, type = "l", col = "red",
     xlab = "time", ylab = "temperature (c)",
     main = "temperature through time")

plot(env_data$time, env_data$precipitation, type = "l", col = "blue",
     xlab = "time", ylab = "precipitation (mm)",
     main = "precipitation through time")

plot(rates$time, rates$speciation_rate, type = "l", col = "green",
     xlab = "time", ylab = "rate",
     main = "speciation rate")

plot(rates$time, rates$extinction_rate, type = "l", col = "orange",
     xlab = "time", ylab = "rate",
     main = "extinction rate")

# tree plot
par(mfrow = c(1, 1))
plotTree(tree, ftype = "off", mar = c(4.1, 4.1, 2.1, 1.1))
title("simulated phylogenetic tree")

# ltt plot
ltt_data <- ltt(tree, plot = FALSE)
plot(ltt_data, log.lineages = TRUE, log = "y",
     col = "blue", lwd = 2, bty = "n",
     main = "lineage-through-time")

# save outputs
save(tree, env_data, rates, file = "environmental_birthdeath_results.RData")
