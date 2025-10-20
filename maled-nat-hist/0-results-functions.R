# Functions used to calculate natural history parameters

# infection incidence rate ----

rate_infection <- function(source_data, outcome_var, ...) {
  
  if(episode == "crude") {
    
    count <- 
      source_data %>% 
      filter({{ outcome_var }} == 1) %>%
      mutate(duration = if_else(max_sev > 0, last_ageend-first_agestart+2, 6, 0)) %>% 
      group_by(pid, month_ss, ...) %>% 
      summarize(infections=sum({{ outcome_var }}), 
                minus=sum(duration))
    
    pt_months <- 
      source_data %>% 
      select(pid, month_ss, ...) %>% 
      distinct() %>% 
      group_by(pid, month_ss) %>% 
      mutate(n_obs = n(), obs=30.5*(1/n_obs))
    
    output_data <- 
      full_join(x=count, y=pt_months, by=join_by("pid", "month_ss", ...)) %>% 
      ungroup() %>% 
      mutate(infections = if_else(is.na(infections),0,infections), 
             minus = if_else(is.na(minus),0,minus)) %>%
      mutate(pt=obs-minus) %>%
      group_by(...) %>%
      summarize(total_inf = sum(infections), total_pt = sum(pt)/30.5) %>% 
      mutate(point = pois.exact(total_inf, total_pt)$rate*100, 
             lower = pois.exact(total_inf, total_pt)$lower*100, 
             upper = pois.exact(total_inf, total_pt)$upper*100) 
    
  } else if(episode == "seq" ) {
    
    count <-
      source_data %>%
      filter(!is.na(first_sid)) %>%
      mutate(duration = last_age-first_age) %>%
      group_by(pid, month_ss, ...) %>%
      summarize(infections=sum({{ outcome_var }}),
                days=sum(duration))

    pt_months <-
      source_data %>%
      select(pid, month_ss, ...) %>%
      distinct() %>%
      group_by(pid, month_ss) %>%
      mutate(n_obs = n(), weight=30.5*(1/n_obs))

    output_data <-
      full_join(x=count, y=pt_months, by=join_by("pid", "month_ss", ...)) %>%
      ungroup() %>%
      mutate(infections = if_else(is.na(infections),0,infections),
             days = if_else(is.na(days),0,days)) %>%
      group_by(pid, ...) %>%
      summarize(pid_inf=sum(infections), pid_obs=sum(weight), pid_minus=sum(days)) %>%
      mutate(pid_pt=pid_obs-pid_minus) %>%
      ungroup() %>%
      group_by(...) %>%
      summarize(total_inf = sum(pid_inf), total_pt = sum(pid_pt)/30.5) %>%
      mutate(point = pois.exact(total_inf, total_pt)$rate*100,
             lower = pois.exact(total_inf, total_pt)$lower*100,
             upper = pois.exact(total_inf, total_pt)$upper*100)
    
  }
  
  output_data
  
}

# incidence rate of pathogen-attributable diarrhea ----

rate_diarrhea <- function(source_data, outcome_var, ...) {
  
  if(episode == "crude") {
    
    count <- 
      source_data %>% 
      filter(inf == 1, {{ outcome_var }} > 0) %>%
      mutate(duration = last_ageend-first_agestart+2) %>% 
      group_by(pid, month_ss, ...) %>% 
      summarize(infections=sum({{ outcome_var }} > 0), 
                minus=sum(duration))
    
    pt_months <- 
      source_data %>% 
      select(pid, month_ss, ...) %>% 
      distinct() %>% 
      group_by(pid, month_ss) %>% 
      mutate(n_obs = n(), obs=30.5*(1/n_obs))
    
    output_data <- 
      full_join(x=count, y=pt_months, by=join_by("pid", "month_ss", ...)) %>% 
      ungroup() %>% 
      mutate(infections = if_else(is.na(infections),0,infections), 
             minus = if_else(is.na(minus),0,minus)) %>%
      mutate(pt=obs-minus) %>%
      group_by(...) %>%
      summarize(total_inf = sum(infections), total_pt = sum(pt)/30.5) %>% 
      mutate(point = pois.exact(total_inf, total_pt)$rate*100, 
             lower = pois.exact(total_inf, total_pt)$lower*100, 
             upper = pois.exact(total_inf, total_pt)$upper*100) 
    
  } else if(episode == "seq" ) {
    
    count <- 
      source_data %>% 
      filter(!is.na(first_sid), {{ outcome_var }} > 0) %>%
      mutate(duration = last_ageend-first_agestart+2) %>% 
      group_by(pid, month_ss, ...) %>% 
      summarize(infections=sum({{ outcome_var }} > 0), 
                days=sum(duration))
    
    pt_months <- 
      source_data %>% 
      select(pid, month_ss, ...) %>% 
      distinct() %>% 
      group_by(pid, month_ss) %>% 
      mutate(n_obs = n(), weight=30.5*(1/n_obs))
    
    output_data <- 
      full_join(x=count, y=pt_months, by=join_by("pid", "month_ss", ...)) %>% 
      ungroup() %>% 
      mutate(infections = if_else(is.na(infections),0,infections), 
             days = if_else(is.na(days),0,days)) %>%
      group_by(pid, ...) %>%
      summarize(pid_inf=sum(infections), pid_obs=sum(weight), pid_minus=sum(days)) %>% 
      mutate(pid_pt=pid_obs-pid_minus) %>% 
      ungroup() %>% 
      group_by(...) %>% 
      summarize(total_inf = sum(pid_inf), total_pt = sum(pid_pt)/30.5) %>% 
      mutate(point = pois.exact(total_inf, total_pt)$rate*100, 
             lower = pois.exact(total_inf, total_pt)$lower*100, 
             upper = pois.exact(total_inf, total_pt)$upper*100) 
    
  }
  
  output_data
  
}

# incidence rate of severe, pathogen-attributable diarrhea ----

rate_severe <- function(source_data, outcome_var, ...) {
  
  if(episode == "crude") {
    
    count <- 
      source_data %>% 
      filter(inf == 1, {{ outcome_var }} == 2) %>%
      mutate(duration = last_ageend-first_agestart+2) %>% 
      group_by(pid, month_ss, ...) %>% 
      summarize(infections=sum({{ outcome_var }} == 2), 
                minus=sum(duration))
    
    pt_months <- 
      source_data %>% 
      select(pid, month_ss, ...) %>% 
      distinct() %>% 
      group_by(pid, month_ss) %>% 
      mutate(n_obs = n(), obs=30.5*(1/n_obs))
    
    output_data <- 
      full_join(x=count, y=pt_months, by=join_by("pid", "month_ss", ...)) %>% 
      ungroup() %>% 
      mutate(infections = if_else(is.na(infections),0,infections), 
             minus = if_else(is.na(minus),0,minus)) %>%
      mutate(pt=obs-minus) %>%
      group_by(...) %>%
      summarize(total_inf = sum(infections), total_pt = sum(pt)/30.5) %>% 
      mutate(point = pois.exact(total_inf, total_pt)$rate*100, 
             lower = pois.exact(total_inf, total_pt)$lower*100, 
             upper = pois.exact(total_inf, total_pt)$upper*100) 
    
  } else if(episode == "seq" ) {
    
    count <- 
      source_data %>% 
      filter(!is.na(first_sid), {{ outcome_var }} == 2) %>%
      mutate(duration = last_ageend-first_agestart+2) %>% 
      group_by(pid, month_ss, ...) %>% 
      summarize(infections=sum({{ outcome_var }} == 2), 
                days=sum(duration))
    
    pt_months <- 
      source_data %>% 
      select(pid, month_ss, ...) %>% 
      distinct() %>% 
      group_by(pid, month_ss) %>% 
      mutate(n_obs = n(), weight=30.5*(1/n_obs))
    
    output_data <- 
      full_join(x=count, y=pt_months, by=join_by("pid", "month_ss", ...)) %>% 
      ungroup() %>% 
      mutate(infections = if_else(is.na(infections),0,infections), 
             days = if_else(is.na(days),0,days)) %>%
      group_by(pid, ...) %>%
      summarize(pid_inf=sum(infections), pid_obs=sum(weight), pid_minus=sum(days)) %>% 
      mutate(pid_pt=pid_obs-pid_minus) %>% 
      ungroup() %>% 
      group_by(...) %>% 
      summarize(total_inf = sum(pid_inf), total_pt = sum(pid_pt)/30.5) %>% 
      mutate(point = pois.exact(total_inf, total_pt)$rate*100, 
             lower = pois.exact(total_inf, total_pt)$lower*100, 
             upper = pois.exact(total_inf, total_pt)$upper*100) 
    
  }
  
  output_data
  
}

# age when first infected ----
age_first <- function(source_data, history_var, age_var, ...) {
  
  output_data <- 
    source_data %>% 
    filter({{ history_var }} == 0, inf==1) %>% 
    group_by(pid) %>%
    arrange(agedays, .by_group = T) %>% 
    filter(row_number()==1) %>%
    group_by(...) %>%
    summarize(point = quantile({{ age_var }}, 0.5, names=F), 
              lower = quantile({{ age_var }}, 0.25, names=F), 
              upper = quantile({{ age_var }}, 0.75, names=F))
  
  output_data
  
}

# proportion of infections that are symptomatic AND pathogen-attributable ----
prop_symptoms <- function(source_data, severity_var, ...) {
  
  output_data <- 
    source_data %>% 
    mutate(numerator = if_else({{ severity_var }} == 0, "no", "yes")) %>% 
    group_by(numerator, ...) %>% 
    summarize(n=n()) %>% 
    pivot_wider(names_from="numerator", values_from="n", values_fill = 0) %>% 
    mutate(denominator=yes+no, 
           point = binom.exact(yes, denominator)$prop, 
           lower = binom.exact(yes, denominator)$lower, 
           upper = binom.exact(yes, denominator)$upper)
  
  output_data
  
}

# proportion of infections that are severe, symptomatic, AND pathogen-attributable ----
prop_severe <- function(source_data, severity_var, ...) {
  
  output_data <- 
    source_data %>% 
    filter({{ severity_var }} > 0) %>%
    mutate(numerator = if_else({{ severity_var }} == 2, "yes", "no")) %>% 
    group_by(numerator, ...) %>% 
    summarize(n=n()) %>% 
    pivot_wider(names_from="numerator", values_from="n", values_fill = 0) %>% 
    mutate(denominator=yes+no,
           point = binom.exact(yes, denominator)$prop, 
           lower = binom.exact(yes, denominator)$lower, 
           upper = binom.exact(yes, denominator)$upper)
  
  output_data
  
}

# duration of infection ----
dur_infection <- function(source_data, begin, end, ...) {
  
  output_data <- 
    source_data %>% 
    mutate(duration = {{ end }} - {{ begin }}) %>% 
    group_by(...) %>% 
    summarize(point = quantile(duration, 0.5, names=F), 
              lower = quantile(duration, 0.25, names=F), 
              upper = quantile(duration, 0.75, names=F))
  
  output_data
}

# duration of attributable diarrhea ----
dur_diarrhea <- function(source_data, severity_var, begin, end, ...) {
  
  output_data <- 
    source_data %>% 
    filter({{ severity_var }} > 0) %>% 
    mutate(duration = {{ end }} - {{ begin }}) %>% 
    group_by(...) %>% 
    summarize(point = quantile(duration, 0.5, names=F), 
              lower = quantile(duration, 0.25, names=F), 
              upper = quantile(duration, 0.75, names=F))
  
  output_data
}