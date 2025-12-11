foreach(i=1:n_sim, .combine = "rbind", .packages = c("tidyverse")) %dopar% {
  
  # set reproducible seed
  seed <- init_seed+(i-1)
  set.seed(seed)

  # generate enrollment information and relevant ages & calendar times
  enrollment_long <- data.frame() # blank dataset
  enrollment_age <- age_info$age_start[1] # pull the age everyone starts at
  censoring_age <- age_info$age_start[nrow(age_info)] # pull the age everyone ends at
  
  for (n in 1:n_participants) {
  
    temp <- 
      age_info %>% # take the input age change data
      mutate(id=n, # assign to the nth participant
             arm = if_else(n %%2 == 0, 1, 0), # mechanism for assigning trial arm
             t_entry = floor(runif(1,1,days_recruitment)), # entry into simulation; all study doses + 14 days already elapsed
             t_start = if_else(age_group == 1, t_entry, # use entry time for first (aka starting) age group
                               t_entry+(age_start-enrollment_age)), # calculate calendar time at start of age group
             t_end = t_start+(age_end-age_start)) # calculate calendar time at end of age group
    
    enrollment_long <- bind_rows(enrollment_long, temp) # add to existing info
    }
  
  risk_sets_long <- data.frame()
  
  for(c in 1:nrow(transmission_info)) {
    
    c_start <- transmission_info$day_start[c]
    c_end <- transmission_info$day_end[c]
    c_rate <- transmission_info$rate[c]
    
    temp <- 
      enrollment_long %>% 
      # determine timing of age changes relative to calendar time
      mutate(type = case_when(
        t_start >= c_start & t_start <= c_end ~ 1, # time in the age group starts during the calendar period
        t_start < c_start & t_end >= c_start & t_end <= c_end ~ 5, # time in the age group starts before the calendar period but ends during the calendar period
        t_start < c_start & t_end > c_end ~ 10, # all time in the calendar period is for the same age group
        TRUE ~ 0), # no time in the age group occurs in the calendar period
        
        # define the calendar start time that applies to time in the age group
        start_me = case_when(
          type==1 ~ t_start, # time in age group starts during calendar period
          type==5 | type==10 ~ c_start, # time in age group started before the calendar period
          TRUE ~ NA), # no time in age group for this calendar period
        
        # define the calendar end time that applies to time in the age group
        stop_me = case_when(
          type==1 & t_end <= c_end ~ t_end, # time in age group ends during the calendar period
          type==1 & t_end > c_end ~ c_end, # time in age group ends after the calendar period
          type==5 ~ t_end, # time in age group ends during the calendar period
          type==10 ~ c_end, # time in age group ends after the calendar period
          TRUE ~ NA), # no time in age group for this calendar period
        
        cal_int=c, h_rate=c_rate) %>% # record the calendar period and population-level transmission rate
      filter(type!=0) # only keep times in age group that occur during the calendar period being evaluated
    
    risk_sets_long <- bind_rows(risk_sets_long, temp) # add to existing info
    
  }
  
  # order each individual's risk periods
  risk_sets_long <- 
    risk_sets_long %>% 
    group_by(id) %>% 
    arrange(age_start, .by_group = T) %>% 
    select(id, arm, age_group, cal_int, t_entry, start_me, stop_me, h_rate, rr_inf, rr_dia, rr_sev)
  
  
  write_csv(
    risk_sets_long %>% 
      mutate(event_id=row_number(), 
             outcome_inf=0, 
             age_event=enrollment_age + (start_me-t_entry)) %>% 
      select(id, arm, age_group, cal_int, outcome_inf, event_id, t_event=start_me, age_event), 
    paste("./Output", pathogen, scenario_inf, scenario_rate, scenario_ve, scenario_hybrid, paste0(file_prefix,"-", formatC(i, width=3, format="d", flag="0"), "_B.csv"), sep = "/"))
  
  # empty data frame to record infections
  observations <- data.frame(id=numeric(), arm=numeric(), event_id=numeric(), t_event=numeric(), age_event=numeric())
  
  # generate infections and severity
  
  for (n in 1:n_participants) {
    
    individual <- 
      risk_sets_long %>% 
      filter(id==n, !is.na(stop_me)) %>% 
      arrange(start_me) # pull the specific risk periods for participant n
    
    # pull time-invariant information
    arm <- individual$arm[1] # can be adapted to change over time by moving to p-loop below
    
    # initialize starting conditions; participants all start the same way
    last_age <- enrollment_age # starting age
    last_t <- individual$start_me[1] # starting calendar time
    begin <- last_t # start of at-risk time
    prior_i <- 0 # no infection history - an assumption - can be modified
    
    for (p in 1:nrow(individual)) {
      
      # pull time-variant information
      h_rate <- individual$h_rate[p] # hazard rate for risk period
      end <- individual$stop_me[p] # end of the homogeneous risk period
      rr_inf <- individual$rr_inf[p] # age-based modification of infection hazard
      rr_dia <- individual$rr_dia[p] # age-based modification of conditional diarrhea risk
      rr_sev <- individual$rr_sev[p] # age-based modification of conditional severe risk
      
      # calculate individual-specific parms - see script with functions for details
      rate <- modify_hr(VEs,arm, rr_inf, prior_i, inf_max)*h_rate # rate used to generate infection time
      rr_symptoms <- modify_rr_symptoms(VEp,arm, rr_dia, prior_i, inf_max) # modification to base probability of developing symptoms, CONDITIONAL on being infected
      rr_severe <- modify_rr_severe(VEh,arm, rr_sev, prior_i, inf_max) # modification to base probability of experiencing a severe infection, CONDITIONAL on being infected and developing symptoms
      
      # determine if infection time should be drawn - controls loop
      p_control <- if_else(rate <= 0, FALSE, TRUE) 
      
      while (p_control) {
        
        Ti <- infection_time(begin, rate) # draw infection time
        
        if (Ti <= end) { # occurs during the risk period under consideration
          
          severity <- infection_severity(p_symptoms*rr_symptoms, p_severe*rr_severe) # determine severity
          last_age <-last_age+(Ti-last_t) # calculate age at infection
          
          observations <- 
            observations %>% 
            add_row(id=n, arm=arm, event_id=severity, t_event=Ti, age_event=last_age) # record infection
          
          begin <- Ti+2 # assume total immunity from next infection for 2 days
          last_t <- Ti # update last t
          prior_i <- prior_i+1 # update infection history
          
          # update parms
          rate <- modify_hr(VEs,arm, rr_inf, prior_i, inf_max)*h_rate
          rr_symptoms <- modify_rr_symptoms(VEp,arm, rr_dia, prior_i, inf_max)
          rr_severe <- modify_rr_severe(VEh,arm, rr_sev, prior_i, inf_max)
          
          # determine if another infection time should be drawn
          p_control <- case_when(
            rate <= 0 ~ FALSE, # no longer at risk
            begin > end ~ FALSE, # resumption of risk occurs after the risk period ends
            TRUE ~ TRUE)
          
          # note that we leave <begin> as is so that, as the next risk period begins, the infection time draw is still appropriately conditioned
          
        } else {
          
          # update info for next iteration
          
          last_age <- last_age+(end-last_t)
          begin <- end+1 # since not infected during risk period, 
          last_t <- end
          
          p_control <- FALSE # end infection loop for this risk period
          
        }
        
        
      }
      
      
    }
    
    gc()
    
  }
  
  observations <- observations %>% mutate(outcome_inf=1)
    
  
  write_csv(observations, paste("./Output", pathogen, scenario_inf, scenario_rate, scenario_ve, scenario_hybrid, paste0(file_prefix,"-", formatC(i, width=3, format="d", flag="0"), "_O.csv"), sep = "/"))
  
  gc()

}
