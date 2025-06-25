#load the simulation functions

source("realistic_species_simulation.R")

quick_results <- run_realistic_simulation(
  n_species = 10, 
  n_lat = 30, 
  n_long = 60, 
  use_real_data = FALSE,
  seed = 123
  
  
  
)

plot_environmental_gradients(quick_results$env_data)
plot_community_patterns(quick_results$env_data, 
                        quick_results$community_richness, 
                        quick_results$community_diversity)


#that was just a short simulation example, now lets try to run a full simulation

full_results <- run_realistic_simulation(
  n_species = 20,
  n_lat = 100,
  n_long = 200,
  use_real_data = FALSE,
  seed = 456
  
  
  
)

plot_environmental_gradients(full_results$env_data)
plot_species_distributions(full_results$env_data, full_results$occurrence_probs)
plot_community_patterns(full_results$env_data, 
                        full_results$community_richness, 
                        full_results$community_diversity)