library(tidyverse)
library(epitools)

source("0-results-functions.R")

# all possible datasets ----
combinations <- expand.grid(
  pathogen=c("adeno", "astro", "campy", "crypto", "noro2", "rota", "sapo", "shig", "stetec", "tepec"), 
  episode=c("crude", "seq"), 
  attr=c("afe", "alg"))

# empty datasets for site-specific results ----
summary <- data.frame()
incidence_bysite <- data.frame()
firstage_bysite <- data.frame()

for (combo in 1:nrow(combinations)) {
  
  pathogen <- as.character(combinations$pathogen[combo])
  episode <- as.character(combinations$episode[combo])
  attr <- as.character(combinations$attr[combo])
  
  if (episode == "crude") { 
    
    full_df <- read_csv(here::here("Outputs", "Data", paste0(paste(pathogen, episode, attr, sep="-"), ".csv")), 
                        col_types = cols(pid=col_character(), country_id=col_character(), region=col_character(), agegroup3=col_character(), agegroup6=col_character(), first_sid=col_character(), prior_inf_cat=col_character(), prior_inf_bin=col_character(), .default = col_number()))
    
  } else {
      
    full_df <- read_csv(here::here("Outputs", "Data", paste0(paste(pathogen, episode, attr, sep="-"), ".csv")), 
                        col_types = cols(pid=col_character(), country_id=col_character(), region=col_character(), agegroup3=col_character(), agegroup6=col_character(), first_sid=col_character(), prior_inf_cat=col_character(), prior_inf_bin=col_character(), path_type=col_character(), .default = col_number()))
     
    }
  
  if (pathogen == "rota") { # only include data from non-introducing countries
    
    full_df <- full_df %>% 
      filter(!country_id %in% c("BR", "PE", "SA"))
    
  }
  
  if(episode == "crude") {
    
    inf_df <- full_df %>% filter(inf == 1)
    
  } else if(episode=="seq") {
    
    inf_df <- full_df %>% filter(!is.na(first_sid))
    
  }
  
  # summary stats: quantify unique infants and total person-time ----
  summary <- bind_rows(summary, 
                       data.frame(pathogen=pathogen, 
                                  episode=episode, 
                                  attr=attr, 
                                  n_infants=length(unique(full_df$pid))))
  
  
  incidence <- 
    bind_rows(# infection incidence
              rate_infection(full_df, inf, country_id) %>% mutate(outcome="infection"), 
              # attributable diarrhea incidence
              rate_diarrhea(full_df, max_sev, country_id) %>% mutate(outcome="diarrhea"), 
              # severe attributable diarrhea
              rate_severe(full_df, max_sev, country_id) %>% mutate(outcome="severe")) %>% 
    mutate(pathogen=pathogen, episode=episode, attr=attr)
  
  incidence_bysite <- bind_rows(incidence_bysite, incidence)
  
  # age at first infection ----
  firstage_bysite <- bind_rows(firstage_bysite, 
                               age_first(full_df, prior_inf_bin, agedays, country_id) %>% 
                                 mutate(pathogen=pathogen, episode=episode, attr=attr), 
                               age_first(full_df, prior_inf_bin, agedays) %>% 
                                 mutate(pathogen=pathogen, episode=episode, attr=attr, country_id="Overall"))
}

summary <- 
  full_join(x=summary, 
            y=incidence_bysite %>% group_by(pathogen, episode, attr, outcome) %>% summarize(pt=sum(total_pt), tot=sum(total_inf)), 
            by=c("pathogen", "episode", "attr"))

# save datasets ----
save(summary, file=here::here("Outputs", "Results", "3.0.2-summary.rdata"))
save(incidence_bysite, file=here::here("Outputs", "Results", "3.0.2-incidence_bysite.rdata"))
save(firstage_bysite, file=here::here("Outputs", "Results", "3.0.2-first_age.rdata"))
