library(tidyverse)
library(epitools)

source("0-results-functions.R")

# all possible datasets ----
combinations <- expand.grid(
  pathogen=c("adeno", "astro", "campy", "crypto", "noro2", "rota", "sapo","shig", "stetec", "tepec"), 
  episode=c("crude", "seq"), 
  attr=c("afe", "alg"))

# initialize empty datasets ----
incidence_byage <- data.frame()
incidence_byboth <- data.frame()
severity_byage <- data.frame()
severity_byboth <- data.frame()
duration_byage <- data.frame()
duration_byboth <- data.frame()

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
  
  
  
  # incidence ----
  incidence_byage <- 
    bind_rows(incidence_byage, 
              rate_infection(full_df, inf, agegroup3) %>% 
                mutate(pathogen=pathogen, episode=episode, attr=attr))
  
  if(episode=="seq") {
    
    incidence_byboth <- 
      bind_rows(incidence_byboth,
                rate_infection(full_df, inf, agegroup3, prior_inf_bin) %>% 
                  mutate(pathogen=pathogen, episode=episode, attr=attr))
  }
  
  # proportion of infections that are attributable  ----
  severity_byage <- 
    bind_rows(severity_byage,
              prop_symptoms(inf_df, max_sev, agegroup6) %>% 
                mutate(pathogen=pathogen, episode=episode, attr=attr, grouping="6mos", outcome="diarrhea") %>%
                rename(agegroup=agegroup6))
  
  if(attr=="alg") {
    
    severity_byage <- 
      bind_rows(severity_byage, 
                prop_symptoms(inf_df, max_sev, agegroup3) %>% 
                  mutate(pathogen=pathogen, episode=episode, attr=attr, grouping="3mos", outcome="diarrhea") %>%
                  rename(agegroup=agegroup3))
    
  }
  
  
  if(episode=="seq") {
    
    severity_byboth <- 
      bind_rows(severity_byboth, 
                prop_symptoms(inf_df, max_sev, agegroup6, prior_inf_bin) %>% 
                  mutate(pathogen=pathogen, episode=episode, attr=attr, grouping="6mos", outcome="diarrhea") %>%
                  rename(agegroup=agegroup6))
    
    if(attr=="alg") {
      
      severity_byboth <- 
        bind_rows(severity_byboth, 
                  prop_symptoms(inf_df, max_sev, agegroup3, prior_inf_bin) %>% 
                    mutate(pathogen=pathogen, episode=episode, attr=attr, grouping="3mos", outcome="diarrhea") %>%
                    rename(agegroup=agegroup3))
    }
  }
  
  
  
  # proportion of infections that are severe and attributable --- 
  severity_byage <- 
    bind_rows(severity_byage, 
              prop_severe(inf_df, max_sev, agegroup6) %>% 
                mutate(pathogen=pathogen, episode=episode, attr=attr, grouping="6mos", outcome="severe") %>%
                rename(agegroup=agegroup6))
  
  if(attr=="alg") {
    
    severity_byage <- 
      bind_rows(severity_byage, 
                prop_severe(inf_df, max_sev, agegroup3) %>% 
                  mutate(pathogen=pathogen, episode=episode, attr=attr, grouping="3mos", outcome="severe") %>%
                  rename(agegroup=agegroup3))
    
  }
  
  if(episode=="seq") {
    
    severity_byboth <- 
      bind_rows(severity_byboth, 
                prop_severe(inf_df, max_sev, agegroup6, prior_inf_bin) %>% 
                  mutate(pathogen=pathogen, episode=episode, attr=attr, grouping="6mos", outcome="severe") %>%
                  rename(agegroup=agegroup6))
    
    if(attr=="alg") {
      
      severity_byboth <- 
        bind_rows(severity_byboth, 
                  prop_severe(inf_df, max_sev, agegroup3, prior_inf_bin) %>% 
                    mutate(pathogen=pathogen, episode=episode, attr=attr, grouping="3mos", outcome="severe") %>%
                    rename(agegroup=agegroup3))
    }
  }
  
  
  
  # duration of infection ---
  
  if(episode == "seq") { 
    
    duration_byage <- 
      bind_rows(duration_byage, 
                dur_infection(inf_df, first_age, last_age, agegroup6) %>% 
                  mutate(pathogen=pathogen, episode=episode, attr=attr, grouping="6mos", outcome="infection") %>%
                  rename(agegroup=agegroup6))
    
    duration_byboth <- 
      bind_rows(duration_byboth, 
                dur_infection(inf_df, first_age, last_age, agegroup6, prior_inf_bin) %>% 
                  mutate(pathogen=pathogen, episode=episode, attr=attr, grouping="6mos", outcome="infection") %>%
                  rename(agegroup=agegroup6))
    
    if(attr=="alg") {
      
      duration_byage <- 
        bind_rows(duration_byage, 
                  dur_infection(inf_df, first_age, last_age, agegroup3) %>% 
                    mutate(pathogen=pathogen, episode=episode, attr=attr, grouping="3mos", outcome="infection") %>%
                    rename(agegroup=agegroup3))
      
      duration_byboth <- 
        bind_rows(duration_byboth, 
                  dur_infection(inf_df, first_age, last_age, agegroup3, prior_inf_bin) %>% 
                    mutate(pathogen=pathogen, episode=episode, attr=attr, grouping="3mos", outcome="infection") %>%
                    rename(agegroup=agegroup3))
      
    }
  }
  
  # duration of symptoms when diarrhea is attributable ----
  duration_byage <- 
    bind_rows(duration_byage, 
              dur_diarrhea(inf_df, max_sev, first_agestart, last_ageend, agegroup6) %>% 
                mutate(pathogen=pathogen, episode=episode, attr=attr, grouping="6mos", outcome="diarrhea") %>%
                rename(agegroup=agegroup6))
  
  if(attr=="alg") {
    
    duration_byage <- 
      bind_rows(duration_byage, 
                dur_diarrhea(inf_df, max_sev, first_agestart, last_ageend, agegroup3) %>% 
                  mutate(pathogen=pathogen, episode=episode, attr=attr, grouping="3mos", outcome="diarrhea") %>%
                  rename(agegroup=agegroup3))
    
  }
  
  if(episode=="seq") {
    
    duration_byboth <- 
      bind_rows(duration_byboth, 
                dur_diarrhea(inf_df, max_sev, first_agestart, last_ageend, agegroup6, prior_inf_bin) %>% 
                  mutate(pathogen=pathogen, episode=episode, attr=attr, grouping="6mos", outcome="diarrhea") %>%
                  rename(agegroup=agegroup6))
    
    if(attr=="alg") {
      
      duration_byboth <- 
        bind_rows(duration_byboth, 
                  dur_diarrhea(inf_df, max_sev, first_agestart, last_ageend, agegroup3, prior_inf_bin) %>% 
                    mutate(pathogen=pathogen, episode=episode, attr=attr, grouping="3mos", outcome="diarrhea") %>%
                    rename(agegroup=agegroup3))
      
    }
    
  }
}

# save datasets ----
save(incidence_byage, file=here::here("Outputs", "Results", "3.0.2-incidence_byage.rdata"))
save(incidence_byboth, file=here::here("Outputs", "Results", "3.0.2-incidence_byboth.rdata"))
save(severity_byage, file=here::here("Outputs", "Results", "3.0.2-severity_byage.rdata"))
save(severity_byboth, file=here::here("Outputs", "Results", "3.0.2-severity_byboth.rdata"))
save(duration_byage, file=here::here("Outputs", "Results", "3.0.2-duration_byage.rdata"))
save(duration_byboth, file=here::here("Outputs", "Results", "3.0.2-duration_byboth.rdata"))
