library(tidyverse)
library(epitools)

source("0-results-functions.R")

# all possible datasets ----
combinations <- expand.grid(
  pathogen=c("adeno", "astro", "campy", "crypto", "noro2", "rota", "sapo", "shig", "stetec", "tepec"), 
  episode=c("crude", "seq"), 
  attr="alg")

# empty datasets for site-specific results ----
summary_sens <- data.frame()
incidence_bysite_sens <- data.frame()
firstage_bysite_sens <- data.frame()
incidence_byage_sens <- data.frame()
severity_byage_sens <- data.frame()

for (combo in 1:nrow(combinations)) {
  
  pathogen <- as.character(combinations$pathogen[combo])
  episode <- as.character(combinations$episode[combo])
  attr <- as.character(combinations$attr[combo])
  
  if (episode == "crude") { 
    
    full_df <- read_csv(here::here("Outputs", "Data", paste0(paste(pathogen, episode, attr, sep="-"), "-sens.csv")), 
                        col_types = cols(pid=col_character(), country_id=col_character(), region=col_character(), agegroup3=col_character(), agegroup6=col_character(), first_sid=col_character(), prior_inf_cat=col_character(), prior_inf_bin=col_character(), .default = col_number()))
    
  } else {
    
    full_df <- read_csv(here::here("Outputs", "Data", paste0(paste(pathogen, episode, attr, sep="-"), "-sens.csv")), 
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
  summary_sens <- bind_rows(summary_sens, 
                       data.frame(pathogen=pathogen, 
                                  episode=episode, 
                                  attr=attr, 
                                  n_infants=length(unique(full_df$pid))))
  
  incidence_bysite_sens <- bind_rows(incidence_bysite_sens, 
                                rate_infection(full_df, inf, country_id) %>% 
                                  mutate(pathogen=pathogen, episode=episode, attr=attr))
  
  # age at first infection ----
  firstage_bysite_sens <- bind_rows(firstage_bysite_sens, 
                               age_first(full_df, prior_inf_bin, agedays, country_id) %>% 
                                 mutate(pathogen=pathogen, episode=episode, attr=attr), 
                               age_first(full_df, prior_inf_bin, agedays) %>% 
                                 mutate(pathogen=pathogen, episode=episode, attr=attr, country_id="Overall"))
  
  # incidence ----
  incidence_byage_sens <- 
    bind_rows(incidence_byage_sens, 
              rate_infection(full_df, inf, agegroup3) %>% 
                mutate(pathogen=pathogen, episode=episode, attr=attr))
  
  # proportion of infections that are attributable  ----
  severity_byage_sens <- 
    bind_rows(severity_byage_sens,
              prop_symptoms(inf_df, max_sev, agegroup6) %>% 
                mutate(pathogen=pathogen, episode=episode, attr=attr, grouping="6mos", outcome="diarrhea") %>%
                rename(agegroup=agegroup6))
  
}

summary_sens <- 
  full_join(x=summary_sens, 
            y=incidence_bysite_sens %>% group_by(pathogen, episode, attr) %>% summarize(pt=sum(total_pt), tot=sum(total_inf)), 
            by=c("pathogen", "episode", "attr"))


# save datasets ----
save(summary_sens, file=here::here("Outputs", "Results", "3.5.1-summary.rdata"))
save(incidence_bysite_sens, file=here::here("Outputs", "Results", "3.5.2-incidence_bysite.rdata"))
save(firstage_bysite_sens, file=here::here("Outputs", "Results", "3.5.3-first_age.rdata"))

save(incidence_byage_sens, file=here::here("Outputs", "Results", "3.5.4-incidence_byage.rdata"))
save(severity_byage_sens, file=here::here("Outputs", "Results", "3.5.5-severity_byage.rdata"))