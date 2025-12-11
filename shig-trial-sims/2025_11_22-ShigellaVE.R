library(tidyverse)
library(doParallel)
library(survival)
library(powerSurvEpi) # used for sample size calculations
library(data.table)

# create a reference to use for 
## 1 - generating file paths & names to read in data, and 
## 2 - generating file names to output data
combinations <- expand.grid(
  scenario_inf = c("InfWorse", "InfEqual"),
  scenario_rate = c("Rate0", "Rate1", "Rate2"), 
  scenario_ve = c("VE2", "VE0"), 
  scenario_hybrid = c("Hybrid0", "Hybrid1")) %>% 
  bind_rows(
    expand.grid(
      scenario_inf = c("Sens-VE"),
      scenario_ve = c("VEsp0.1", "VEsp0.2", "VEsp0.3", "VEsp0.4", "VEsp0.5", 
                      "VEsp0.6", "VEsp0.7", "VEsp0.8", "VEsp0.9")))

# common information across all data to be analyzed
pathogen <- "Shigella"
age_start <- 273
age_fu12 <- age_start+365
age_fu24 <- age_fu12+365
age_end <- 1095

# if running on a personal laptop; create cluster and determine number of available nodes
# if(detectCores() >= 7) {
#   nocores <- 5
#   cl <- makeCluster(nocores)
#   registerDoParallel(cl)
# } else {
#   errorCondition("Cannot process five at a time; adjust stop increment")
# }

# if running on cluster, make sure to request n+1 cores
# where n is a number that divides the total number of simulations evenly
nocores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))-1
cl <- makeCluster(nocores)
registerDoParallel(cl)

for (combo in 1:nrow(combinations)) {

  # determine parameter scenario ----
  scenario_inf <- as.character(combinations$scenario_inf[combo])
  scenario_rate <- as.character(combinations$scenario_rate[combo])
  scenario_ve <- as.character(combinations$scenario_ve[combo])
  scenario_hybrid <- as.character(combinations$scenario_hybrid[combo])
  
  # specify different naming conventions for each major scenario set (3 total)
  file_prefix <- case_when(scenario_inf == "Sens-VE" ~ paste("Shig", scenario_ve, sep="-"), 
                           scenario_inf == "InfWorse" ~ paste("Shig-Worse", scenario_rate, scenario_ve,
                                                              scenario_hybrid, sep="-"), 
                           TRUE ~ paste("Shig-Equal", scenario_rate, scenario_ve,
                                        scenario_hybrid, sep="-"))
  
  # specify the i_th run to begin and end data wrangling with 
  start <- 1 
  stop <- nocores
  
  while(start < 500) { # control loop so it stops once max number of simulations is hit
    
    write_me <- foreach(i=start:stop, .combine="rbind", .packages=c("tidyverse", "powerSurvEpi", "survival", "data.table", "here"), .inorder=F) %dopar% {
      
      results <- data.frame()
      
      if (scenario_inf == "Sens-VE") {
        
        infections <- fread(here::here("Output", pathogen, scenario_inf, 
                                    paste0(file_prefix, "-", formatC(i, width=3, format="d", flag="0"), "_O.csv")))
        
        # because we will use age as the underlying time scale for survival analysis, 
        # we only need the final observation (i.e., age when follow-up ends)
        # if there are changes that occur over calendar time, we need all time points
        end <- fread(here::here("Output", pathogen, scenario_inf, 
                                paste0(file_prefix, "-", formatC(i, width=3, format="d", flag="0"), "_B.csv"))) %>% 
          select(id, arm, event_id, t_event, age_event, outcome_inf) %>%
          filter(age_event == age_end) 
        
      } else {
        
        infections <- fread(here::here("Output", pathogen, scenario_inf, scenario_rate, scenario_ve, scenario_hybrid, 
                                    paste0(file_prefix, "-", formatC(i, width=3, format="d", flag="0"), "_O.csv")))
        
        end <- fread(here::here("Output", pathogen, scenario_inf, scenario_rate, scenario_ve, scenario_hybrid, 
                                paste0(file_prefix, "-", formatC(i, width=3, format="d", flag="0"), "_B.csv"))) %>% 
          select(id, arm, event_id, t_event, age_event, outcome_inf) %>%
          filter(age_event == age_end)
        
      }
      
      # create data for a perfect trial ----
      # this is the full set of individuals (n=10,000) with perfect data collection
      # i.e. all infections are known and accurately documented
      # this allows us to create follow-up intervals based on infections
      # everyone starts at the same age for the first interval, which ends at the first infection
      # since each infection confers temporary, total immunity for some x days 
      # the next interval starts x days after the recorded day of the first infection
      # this repeats until the censoring age
      
       perfect <- 
        bind_rows(infections, end) %>% 
        group_by(id) %>% 
        arrange(age_event, .by_group = T) %>% 
        mutate(
          t_start = case_when(
            row_number() == 1 ~ age_start, 
            lag(outcome_inf) == 1 ~ lag(age_event) + 2, # +2 reflects two days of temporary perfect immunity
            TRUE ~ lag(age_event)), 
          
          t_end = case_when(
            age_event == t_start ~ age_event+0.5, # offset so some follow-up time is recorded
            t_start >= age_end ~ -1, # handles when censoring age happens during the window of immunity
            TRUE ~ age_event), 
          
          # trial outcomes - event_id is severity when outcome_inf==1
          outcome_dia = if_else(outcome_inf == 1 & event_id > 1, 1, 0), 
          outcome_sev = if_else(outcome_inf == 1 & event_id == 3, 1, 0), 
          
          # adjust so current infection is not counted
          prior = cumsum(outcome_inf)-outcome_inf, 
          prior_cat = if_else(prior < 2, as.character(prior), "2+"), 
          sim=i)
      
      gc()
      
      # create large trial with realistic data collection ----
      # this is the full set of individuals (n=10,000) with symptom-based stool collection
      # i.e. only diarrheal episodes are known and accurately documented
      # here we create follow-up intervals based on diarrheal episodes
      # otherwise principles of follow-up intervals are the same as for the perfect trial
      symptom_based <- 
        perfect %>% 
        filter(outcome_dia == 1 | age_event == age_end) %>% 
        group_by(id) %>% 
        arrange(age_event, .by_group = T) %>% 
        mutate(
          t_start = case_when(
            row_number() == 1 ~ age_start,
            lag(outcome_dia) == 1 ~ lag(age_event) + 2,
            TRUE ~ lag(age_event)), 
          
          t_end = case_when(
            age_event == t_start ~ age_event+0.5,
            t_start >= age_end ~ -1, 
            TRUE ~ age_event),
          
          # adjust so current episode is not counted
          prior = cumsum(outcome_dia)-outcome_dia, 
          prior_cat = if_else(prior < 2, as.character(prior), "2+")) %>% 
        filter(t_end != -1)
      
      # now remove censoring age that happens during period of immunity
      # can't remove it beforehand since this could be due to an asymptomatic infection
      # which would not be observed in the realistic trial
      perfect <- perfect %>% filter(t_end != -1)
      gc()
      
      # VE for perfect trial (dat=100)----
      
      # now we calculate our VE estimates of interest
      # two outcomes of interest: 
      ## outcome=2: Shigella diarrhea
      ## outcome=3: severe Shigella diarrhea
      # three analytic approaches for each
      ## ana=2: analysis where hazard strata are constructed based on *infection* history during the trial
      ## ana=1: analysis where data after the first *infection* is discarded 
      ## (since asymptomatic infections also confer immunity)
      ## ana=0: crude analysis (no hazard strata and all infections are included)
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_dia)~factor(arm)+strata(prior_cat)+cluster(id), 
                     perfect, robust=T)
      
      results <- bind_rows(results, 
                           data.frame(outcome=2, ana=2, dat=100,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_dia)~factor(arm), 
                     perfect %>% filter(prior==0))
      
      results <- bind_rows(results, 
                           data.frame(outcome=2, ana=1, dat=100,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_dia)~factor(arm)+cluster(id), perfect, robust=T)
      
      results <- bind_rows(results, 
                           data.frame(outcome=2, ana=0, dat=100,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_sev)~factor(arm)+strata(prior_cat)+cluster(id), 
                     perfect, robust=T)
      
      results <- bind_rows(results, 
                           data.frame(outcome=3, ana=2, dat=100,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_sev)~factor(arm), 
                     perfect %>% filter(prior==0))
      
      results <- bind_rows(results, 
                           data.frame(outcome=3, ana=1, dat=100,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_sev)~factor(arm)+cluster(id), perfect, robust=T)
      
      results <- bind_rows(results, 
                           data.frame(outcome=3, ana=0, dat=100,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      # VE for large trial with symptom-based reporting (dat=80)----
      
      # same outcomes as the perfect trial, but different analytic approaches
      ## ana=2: analysis where hazard strata are constructed based on *diarrhea* history during the trial
      ## ana=1: analysis where data after the first *diarrheal episode* is discarded 
      ## ana=0: crude analysis (no hazard strata and all diarrheal episodes are included)
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_dia)~factor(arm)+strata(prior_cat)+cluster(id), 
                     symptom_based, robust=T)
      
      results <- bind_rows(results, 
                           data.frame(outcome=2, ana=2, dat=80,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_dia)~factor(arm), 
                     symptom_based %>% filter(prior==0))
      
      results <- bind_rows(results, 
                           data.frame(outcome=2, ana=1, dat=80,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_dia)~factor(arm)+cluster(id), symptom_based, robust=T)
      
      results <- bind_rows(results, 
                           data.frame(outcome=2, ana=0, dat=80,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_sev)~factor(arm)+strata(prior_cat)+cluster(id), 
                     symptom_based, robust=T)
      
      results <- bind_rows(results, 
                           data.frame(outcome=3, ana=2, dat=80,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_sev)~factor(arm), 
                     symptom_based %>% filter(prior==0))
      
      results <- bind_rows(results, 
                           data.frame(outcome=3, ana=1, dat=80,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_sev)~factor(arm)+cluster(id), symptom_based, robust=T)
      
      results <- bind_rows(results, 
                           data.frame(outcome=3, ana=0, dat=80,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      # VE for reasonably-sized trial with active detection of infections (dat=20)----
      # same analysis and outcomes as for the perfect trial
      # however, first we do a sample size calculation powered on the rarer outcome (severe Shigella diarrhea)
      # we use the last 400 infant IDs to represent a preliminary study
      # the package for sample size calculation returns an n *per arm*
      # then we filter the perfect trial data to include the first n times 2 infant IDs
      # different sample sizes based on what kind of censoring is done
      
      # no censoring models
      per_arm <- ssizeCT(formula = Surv(time=t_start, time2=t_end, event=outcome_sev) ~ arm+cluster(id), 
                         dat = perfect %>% 
                           filter(id >= 4601) %>% 
                           mutate(arm=factor(arm, levels=c(0,1), labels=c("C","E"))), 
                         power=0.8, k=1, RR=0.4, alpha=0.05)$ssize[1]
      
      smaller <- perfect %>% filter(id <= (per_arm*2))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_dia)~factor(arm)+strata(prior_cat)+cluster(id), 
                     smaller, robust=T)
      
      results <- bind_rows(results, 
                           data.frame(outcome=2, ana=2, dat=20, size=per_arm,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_dia)~factor(arm)+cluster(id), smaller, robust=T)
      
      results <- bind_rows(results, 
                           data.frame(outcome=2, ana=0, dat=20, size=per_arm,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_sev)~factor(arm)+strata(prior_cat)+cluster(id), 
                     smaller, robust=T)
      
      results <- bind_rows(results, 
                           data.frame(outcome=3, ana=2, dat=20, size=per_arm,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_sev)~factor(arm)+cluster(id), smaller, robust=T)
      
      results <- bind_rows(results, 
                           data.frame(outcome=3, ana=0, dat=20, size=per_arm,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      # single event model (censoring at first infection) 
      
      per_arm <- ssizeCT(formula = Surv(time=t_start, time2=t_end, event=outcome_sev) ~ arm, 
                         dat = perfect %>% 
                           filter(id >= 4601, prior == 0) %>% 
                           mutate(arm=factor(arm, levels=c(0,1), labels=c("C","E"))), 
                         power=0.8, k=1, RR=0.4, alpha=0.05)$ssize[1]
      
      smaller <- perfect %>% filter(id <= (per_arm*2), prior==0)
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_dia)~factor(arm), smaller)
      
      results <- bind_rows(results, 
                           data.frame(outcome=2, ana=1, dat=20, size=per_arm,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_sev)~factor(arm), smaller)
      
      results <- bind_rows(results, 
                           data.frame(outcome=3, ana=1, dat=20, size=per_arm,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      # VE for reasonably-sized trial with symptom-based reporting (dat=10)----
      # same analysis and outcomes as for the large trial with symptom-based reporting
      # sample size calculation is done the same way as for VE from a reasonably-sized trial 
      # with active detection of infections
      
      # recurrent event model (censoring at end of follow-up)
      per_arm <- ssizeCT(formula = Surv(time=t_start, time2=t_end, event=outcome_sev) ~ arm+cluster(id), 
                         dat = symptom_based %>% 
                           filter(id >= 4601) %>% 
                           mutate(arm=factor(arm, levels=c(0,1), labels=c("C","E"))), 
                         power=0.8, k=1, RR=0.4, alpha=0.05)$ssize[1]
      
      smaller <- symptom_based %>% filter(id <= (per_arm*2))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_dia)~factor(arm)+strata(prior_cat)+cluster(id), 
                     smaller, robust=T)
      
      results <- bind_rows(results, 
                           data.frame(outcome=2, ana=2, dat=10, size=per_arm,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_dia)~factor(arm)+cluster(id), smaller, robust=T)
      
      results <- bind_rows(results, 
                           data.frame(outcome=2, ana=0, dat=10, size=per_arm,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_sev)~factor(arm)+strata(prior_cat)+cluster(id), 
                     smaller, robust=T)
      
      results <- bind_rows(results, 
                           data.frame(outcome=3, ana=2, dat=10, size=per_arm,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_sev)~factor(arm)+cluster(id), smaller, robust=T)
      
      results <- bind_rows(results, 
                           data.frame(outcome=3, ana=0, dat=10, size=per_arm,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      # single event model (censoring at first infection) 
      per_arm <- ssizeCT(formula = Surv(time=t_start, time2=t_end, event=outcome_sev) ~ arm, 
                         dat = symptom_based %>% 
                           filter(id >= 4601, prior == 0) %>% 
                           mutate(arm=factor(arm, levels=c(0,1), labels=c("C","E"))), 
                         power=0.8, k=1, RR=0.4, alpha=0.05)$ssize[1]
      
      smaller <- symptom_based %>% filter(id <= (per_arm*2), prior==0)
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_dia)~factor(arm), smaller)
      
      results <- bind_rows(results, 
                           data.frame(outcome=2, ana=1, dat=10, size=per_arm,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_sev)~factor(arm), smaller)
      
      results <- bind_rows(results, 
                           data.frame(outcome=3, ana=1, dat=10, size=per_arm,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      results <- results %>% mutate(sim=i)
      
      results
    }
    
    if(start == 1) {
      write_csv(write_me, file = here::here("Processed", pathogen, paste(file_prefix, "VE-Results.csv", sep="-")), append=F)
    } else {
      write_csv(write_me, file = here::here("Processed", pathogen, paste(file_prefix, "VE-Results.csv", sep="-")), append=T)
    }
    
    gc()
    
    start <- start+nocores
    stop <- stop+nocores
    
  }
  
  gc()
}

# age chunk stratified ----
# doing only for the scenario in figure 1 (infection=vaccination, no ve vs infection, no hybrid immunity)
# also then only need severe outcome
strat_combo <- combinations %>% filter(scenario_inf=="InfEqual", scenario_ve=="VE0", scenario_hybrid=="Hybrid0")
for (combo in 1:nrow(strat_combo)) {
  
  # determine parameter scenario ----
  scenario_inf <- as.character(strat_combo$scenario_inf[combo])
  scenario_rate <- as.character(strat_combo$scenario_rate[combo])
  scenario_ve <- as.character(strat_combo$scenario_ve[combo])
  scenario_hybrid <- as.character(strat_combo$scenario_hybrid[combo])
  
  # specify different naming conventions for each major scenario set (3 total)
  file_prefix <- paste("Shig-Equal", scenario_rate, scenario_ve, scenario_hybrid, sep="-")
  
  # specify the i_th run to begin and end data wrangling with 
  start <- 1 
  stop <- nocores
  
  while(start < 500) { # control loop so it stops once max number of simulations is hit
    
    write_me <- foreach(i=start:stop, .combine="rbind", .packages=c("tidyverse", "powerSurvEpi", "survival", "data.table", "here"), .inorder=F) %dopar% {
      
      results <- data.frame()
      
      infections <- fread(here::here("Output", pathogen, scenario_inf, scenario_rate, scenario_ve, scenario_hybrid, 
                                     paste0(file_prefix, "-", formatC(i, width=3, format="d", flag="0"), "_O.csv"))) %>% 
        filter(age_event <= age_fu24)
      
      end <- fread(here::here("Output", pathogen, scenario_inf, scenario_rate, scenario_ve, scenario_hybrid, 
                              paste0(file_prefix, "-", formatC(i, width=3, format="d", flag="0"), "_B.csv"))) %>% 
        select(id, arm, event_id, t_event, age_event, outcome_inf) %>%
        filter(age_event == age_end) %>% 
        mutate(age_event = age_fu24)
      
      # create data for a perfect trial ----
      # this is the full set of individuals (n=10,000) with perfect data collection
      # i.e. all infections are known and accurately documented
      # this allows us to create follow-up intervals based on infections
      # everyone starts at the same age for the first interval, which ends at the first infection
      # since each infection confers temporary, total immunity for some x days 
      # the next interval starts x days after the recorded day of the first infection
      # this repeats until the censoring age
      
      perfect <- 
        bind_rows(infections, end) %>% 
        filter(age_event<=age_fu24) %>% 
        group_by(id) %>% 
        arrange(age_event, .by_group = T) %>% 
        mutate(
          t_start = case_when(
            row_number() == 1 ~ age_start, 
            lag(outcome_inf) == 1 ~ lag(age_event) + 2, # +2 reflects two days of temporary perfect immunity
            TRUE ~ lag(age_event)), 
          
          t_end = case_when(
            age_event == t_start ~ age_event+0.5, # offset so some follow-up time is recorded
            t_start >= age_fu24 ~ -1, # handles when censoring age happens during the window of immunity
            TRUE ~ age_event), 
          
          # trial outcomes - event_id is severity when outcome_inf==1
          outcome_dia = if_else(outcome_inf == 1 & event_id > 1, 1, 0), 
          outcome_sev = if_else(outcome_inf == 1 & event_id == 3, 1, 0), 
          
          # adjust so current infection is not counted
          prior = cumsum(outcome_inf)-outcome_inf, 
          prior_cat = if_else(prior < 2, as.character(prior), "2+"), 
          sim=i)
      
      gc()
      
      # create large trial with realistic data collection ----
      # this is the full set of individuals (n=10,000) with symptom-based stool collection
      # i.e. only diarrheal episodes are known and accurately documented
      # here we create follow-up intervals based on diarrheal episodes
      # otherwise principles of follow-up intervals are the same as for the perfect trial
      symptom_based <- 
        perfect %>% 
        filter(outcome_dia == 1 | age_event == age_fu24) %>% 
        group_by(id) %>% 
        arrange(age_event, .by_group = T) %>% 
        mutate(
          t_start = case_when(
            row_number() == 1 ~ age_start,
            lag(outcome_dia) == 1 ~ lag(age_event) + 2,
            TRUE ~ lag(age_event)), 
          
          t_end = case_when(
            age_event == t_start ~ age_event+0.5,
            t_start >= age_fu24 ~ -1, 
            TRUE ~ age_event),
          
          # adjust so current episode is not counted
          prior = cumsum(outcome_dia)-outcome_dia, 
          prior_cat = if_else(prior < 2, as.character(prior), "2+")) %>% 
        filter(t_end != -1)
      
      # now remove censoring age that happens during period of immunity
      # can't remove it beforehand since this could be due to an asymptomatic infection
      # which would not be observed in the realistic trial
      perfect <- perfect %>% filter(t_end != -1)
      gc()
      
      # creating separate datasets for first and second year of life
      perfect12 <- 
        perfect %>% 
        mutate(keep = case_when(
          age_event < age_fu12 ~ 1, 
          t_start <= age_fu12 & age_fu12 <= t_end ~ 1, 
          TRUE ~ 0)) %>% 
        filter(keep==1) %>% 
        mutate(t_end = if_else(t_end > age_fu12, age_fu12, t_end)) %>% 
        filter(t_start != t_end)
      
      perfect24 <- 
        perfect %>% 
        mutate(keep = case_when(
          age_event < age_fu12 ~ 0, 
          t_start <= age_fu12 & age_fu12 <= t_end ~ 1, 
          TRUE ~ 1)) %>% 
        filter(keep==1) %>% 
        mutate(t_start = if_else(t_start < age_fu12, age_fu12, t_start)) %>% 
        filter(t_start != t_end)
      
      symptom_based12 <- 
        symptom_based %>% 
        mutate(keep = case_when(
          age_event < age_fu12 ~ 1, 
          t_start <= age_fu12 & age_fu12 <= t_end ~ 1, 
          TRUE ~ 0)) %>% 
        filter(keep==1) %>% 
        mutate(t_end = if_else(t_end > age_fu12, age_fu12, t_end)) %>% 
        filter(t_start != t_end)
      
      symptom_based24 <- 
        symptom_based %>% 
        mutate(keep = case_when(
          age_event < age_fu12 ~ 0, 
          t_start <= age_fu12 & age_fu12 <= t_end ~ 1, 
          TRUE ~ 1)) %>% 
        filter(keep==1) %>% 
        mutate(t_start = if_else(t_start < age_fu12, age_fu12, t_start)) %>% 
        filter(t_start != t_end)
      
      rm(perfect, symptom_based)
      gc()
      
      # VE for perfect trial----
      #dat=120 for first year analyses and dat=240 for second year analyses
      
      # now we calculate our VE estimates of interest
      # two outcomes of interest: 
      ## outcome=2: Shigella diarrhea
      ## outcome=3: severe Shigella diarrhea
      # three analytic approaches for each
      ## ana=2: analysis where hazard strata are constructed based on *infection* history during the trial
      ## ana=1: analysis where data after the first *infection* is discarded 
      ## (since asymptomatic infections also confer immunity)
      ## ana=0: crude analysis (no hazard strata and all infections are included)
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_sev)~factor(arm)+strata(prior_cat)+cluster(id), 
                     perfect12, robust=T)
      
      results <- bind_rows(results, 
                           data.frame(outcome=3, ana=2, dat=120,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_sev)~factor(arm)+strata(prior_cat)+cluster(id), 
                     perfect24, robust=T)
      
      results <- bind_rows(results, 
                           data.frame(outcome=3, ana=2, dat=240,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_sev)~factor(arm), 
                     perfect12 %>% filter(prior==0))
      
      results <- bind_rows(results, 
                           data.frame(outcome=3, ana=1, dat=120,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_sev)~factor(arm), 
                     perfect24 %>% filter(prior==0))
      
      results <- bind_rows(results, 
                           data.frame(outcome=3, ana=1, dat=240,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_sev)~factor(arm)+cluster(id), perfect12, robust=T)
      
      results <- bind_rows(results, 
                           data.frame(outcome=3, ana=0, dat=120,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_sev)~factor(arm)+cluster(id), perfect24, robust=T)
      
      results <- bind_rows(results, 
                           data.frame(outcome=3, ana=0, dat=240,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      # VE for large trial with symptom-based reporting----
      # dat=12 for first year analyses and dat=24 for second year analyses
      
      # same outcomes as the perfect trial, but different analytic approaches
      ## ana=2: analysis where hazard strata are constructed based on *diarrhea* history during the trial
      ## ana=1: analysis where data after the first *diarrheal episode* is discarded 
      ## ana=0: crude analysis (no hazard strata and all diarrheal episodes are included)
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_sev)~factor(arm)+strata(prior_cat)+cluster(id), 
                     symptom_based12, robust=T)
      
      results <- bind_rows(results, 
                           data.frame(outcome=3, ana=2, dat=12,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_sev)~factor(arm)+strata(prior_cat)+cluster(id), 
                     symptom_based24, robust=T)
      
      results <- bind_rows(results, 
                           data.frame(outcome=3, ana=2, dat=24,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_sev)~factor(arm), 
                     symptom_based12 %>% filter(prior==0))
      
      results <- bind_rows(results, 
                           data.frame(outcome=3, ana=1, dat=12,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_sev)~factor(arm), 
                     symptom_based24 %>% filter(prior==0))
      
      results <- bind_rows(results, 
                           data.frame(outcome=3, ana=1, dat=24,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_sev)~factor(arm)+cluster(id), symptom_based12, robust=T)
      
      results <- bind_rows(results, 
                           data.frame(outcome=3, ana=0, dat=12,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      model <- coxph(Surv(time=t_start, time2=t_end, event=outcome_sev)~factor(arm)+cluster(id), symptom_based24, robust=T)
      
      results <- bind_rows(results, 
                           data.frame(outcome=3, ana=0, dat=24,
                                      beta=model$coefficients[1], 
                                      lower=confint(model)[1], upper=confint(model)[2]))
      
      results <- results %>% mutate(sim=i)
      
      results
    }
    
    write_csv(write_me, file = here::here("Processed", "Shigella", paste(file_prefix, "VE-Results.csv", sep="-")), append = T)
    
    gc()
    
    start <- start+nocores
    stop <- stop+nocores
    
  }


  gc()
}
