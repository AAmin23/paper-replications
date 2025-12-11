# Scenario combinations
combinations <- expand.grid(
  scenario_hybrid = c("Hybrid0", "Hybrid1"), 
  scenario_ve = c("VE2", "VE0"), 
  scenario_rate = c("Rate0", "Rate1", "Rate2")) 

# Infection-conferred immunity is inferior to vaccine-conferred immunity ----
source("./2024_05_02-InfWorse-SimFxs.R") # functions
scenario_inf <- "InfWorse"

for (combo in 1:nrow(combinations)) {
  scenario_rate <- as.character(combinations$scenario_rate[combo])
  scenario_ve <- as.character(combinations$scenario_ve[combo])
  scenario_hybrid <- as.character(combinations$scenario_hybrid[combo])
  
  VEs <- if_else(scenario_ve == "VE2", VEs_20, VEs_0)
  inf_max <- if_else(scenario_hybrid == "Hybrid0", 0, 1)
  
  if(scenario_rate == "Rate0") {
    transmission_info <- rate0
  } else if(scenario_rate == "Rate1") {
    transmission_info <- rate1
  } else if(scenario_rate == "Rate2" ) {
    transmission_info <- rate2
  }
  
  VEp <- 1-((1-VEsp)/(1-VEs))
  VEh <- 1-((1-VEsph)/((1-VEs)*(1-VEp)))
  file_prefix <- paste(substr(pathogen,1,4), substr(scenario_inf,4,8), scenario_rate, scenario_ve, scenario_hybrid, sep="-")
  
  source("./2024_04_11-Simulation.R") # run with parameter combination

  gc()
}

# Infection-conferred immunity is equal to vaccine-conferred immunity ----
source("./2024_04_11-InfEqual-SimFxs.R") # functions
scenario_inf <- "InfEqual"

for (combo in 1:nrow(combinations)) {
  scenario_rate <- as.character(combinations$scenario_rate[combo])
  scenario_ve <- as.character(combinations$scenario_ve[combo])
  scenario_hybrid <- as.character(combinations$scenario_hybrid[combo])
  
  VEs <- if_else(scenario_ve == "VE2", VEs_20, VEs_0)
  inf_max <- if_else(scenario_hybrid == "Hybrid0", 0, 1)
  
  if(scenario_rate == "Rate0") {
    transmission_info <- rate0
  } else if(scenario_rate == "Rate1") {
    transmission_info <- rate1
  } else if(scenario_rate == "Rate2" ) {
    transmission_info <- rate2
  }
  
  VEp <- 1-((1-VEsp)/(1-VEs))
  VEh <- 1-((1-VEsph)/((1-VEs)*(1-VEp)))
  file_prefix <- paste(substr(pathogen,1,4), substr(scenario_inf,4,8), scenario_rate, scenario_ve, scenario_hybrid, sep="-")
  
  source("./2024_04_11-Simulation.R") # run with parameter combination
  
  gc()
}

# Examining bias with different VEs against diarrhea when immunity from vaccination and infection are comparable
  values <- c(seq(0.1, 0.3, by=0.1), seq(0.5, 0.9, by=0.1))
  
  for (combo in 1:length(values)) {
    
    scenario_rate <- "Rate1"
    scenario_ve <- "VE0"
    scenario_hybrid <- "Hybrid0"
    
    VEs <- VEs_0
    inf_max <- 0
    transmission_info <- rate1
    
    VEsp <- values[combo]
    VEsph <- VEsp # no attentuation of severity - but not of interest here
    
    VEp <- 1-((1-VEsp)/(1-VEs))
    VEh <- 1-((1-VEsph)/((1-VEs)*(1-VEp)))
    file_prefix <- paste0("Shig-VEsp", VEsp)
    
    source("./2024_04_11-Simulation.R") # run with parameter combination
    
    gc()
      
  }

