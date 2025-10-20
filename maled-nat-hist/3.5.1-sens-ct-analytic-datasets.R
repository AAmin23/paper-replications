library(tidyverse)
tac <- read_csv(here::here("Outputs", "Data", "3.5-path-attr-data.csv"))

# custom function for unique infections; crude approach ----
create_crude <- function(dataset, var_infection, var_attribution) {
  
  temp <- 
    dataset %>% 
    select(pid, sid, obs_n, agedays, agestart, ageend, month_ss, stooltype, d1_location, score_binary, {{ var_infection }}, {{ var_attribution }}) %>% 
    group_by(pid) %>% 
    arrange(agedays, .by_group = T) %>%
    mutate(group = case_when(
      d1_location == "Lag" ~ lag(obs_n), 
      lag(d1_location) == "Lead" ~ lag(obs_n), 
      d1_location == "Both" ~ lag(obs_n), 
      lag(d1_location) == "Both" ~ lag(obs_n,2), 
      TRUE ~ obs_n)) 
  
  one_pos <- 
    temp %>% 
    filter({{ var_infection }} == 1) %>%
    group_by(pid, group) %>% 
    filter(n()==1) %>% 
    mutate(
      first_sid=sid,
      
      # for attributed diarrhea only - duration of illness
      first_agestart=if_else(stooltype=="D1" & {{ var_attribution }}==1, agestart, NA_real_),
      last_ageend=if_else(stooltype=="D1" & {{ var_attribution }}==1, ageend, NA_real_),
      max_sev=if_else(stooltype=="D1" & {{ var_attribution }}==1, score_binary+1, 0), 
      # for all episodes - only consider diarrheal ep info if attributable
      first_age=pmin(agedays, first_agestart, na.rm=T), 
      last_age=pmax(agedays, last_ageend, na.rm=T),
      first_month=month_ss)
  
  two_pos <- 
    temp %>% 
    filter({{ var_infection }}==1) %>% 
    group_by(pid, group) %>% 
    filter(n()==2) %>% 
    mutate(
      keep = case_when(
        stooltype=="D1" & {{ var_attribution }}==1 ~ 1,
        stooltype=="M1" & any({{ var_attribution }}==0) ~ 1, 
        TRUE~0), # want diarrheal specimen retained if attributable
      
      # want earliest and latest possible ages
      # for attributed diarrhea only - duration of illness
      first_agestart=if_else(any(stooltype=="D1") & any({{ var_attribution }}==1), min(agestart, na.rm=T), NA_real_),
      last_ageend=if_else(any(stooltype=="D1") & any({{ var_attribution }}==1), max(ageend, na.rm=T), NA_real_),
      max_sev=if_else(any(stooltype=="D1") & any({{ var_attribution }}==1), max(score_binary,na.rm=T)+1,0,0),
      
      # for all episodes - only consider diarrheal ep info if attributable and set to one day duration otherwise to flag sub-clin infections
      age1=min(agedays), age2=max(agedays), 
      first_age=pmin(age1,first_agestart,na.rm=T),
      last_age=if_else({{ var_attribution }}==1, pmax(age2,last_ageend), first_age, first_age),
      first_month=min(month_ss), 
      first_sid=if_else(obs_n==min(obs_n),sid,NA_character_)) %>% 
    fill(first_sid, .direction="updown") %>% 
    filter(keep==1) %>% 
    select(-keep, -age1, -age2)
  
  # the rare ones that have pos D1-M1-D1 within a short period - treat as separate episodes of diarrhea and ignore monthly stool sample
  mult_pos <- 
    temp %>% 
    filter({{ var_infection }}==1) %>% 
    group_by(pid, group) %>% 
    filter(n()>2, stooltype=="D1") %>% 
    ungroup() %>% 
    mutate(
      first_sid=sid,
      
      # for attributed diarrhea only - duration of illness
      first_agestart=if_else(stooltype=="D1" & {{ var_attribution }}==1, agestart, NA_real_),
      last_ageend=if_else(stooltype=="D1" & {{ var_attribution }}==1, ageend, NA_real_),
      max_sev=if_else(stooltype=="D1" & {{ var_attribution }}==1, score_binary+1, 0), 
      # for all episodes - only consider diarrheal ep info if attributable
      first_age=pmin(agedays, first_agestart, na.rm=T), 
      last_age=pmax(agedays, last_ageend, na.rm=T),
      first_month=month_ss)
  
  all_pos <- 
    bind_rows(one_pos, two_pos, mult_pos) %>% 
    ungroup() %>% 
    group_by(pid) %>% 
    arrange(agedays, .by_group=T) %>% 
    select(-sid, -obs_n, -agestart, -ageend, -stooltype, -d1_location, -score_binary, -group) %>% 
    rename(inf={{ var_infection }})
  
  obs <- 
    dataset %>% 
    filter(!is.na({{ var_infection }})) %>% 
    select(pid, country_id, region, month_ss, agegroup3, agegroup6) %>% 
    distinct()
  
  final <- 
    left_join(x=obs, y=all_pos, by=c("pid", "month_ss")) %>% 
    group_by(pid) %>% 
    arrange(month_ss, agedays, .by_group=T) %>% 
    mutate(inf = if_else(is.na(inf), 0, inf), 
           prior_inf=cumsum(inf)-inf, 
           prior_inf_cat = if_else(prior_inf > 1, "2+", as.character(prior_inf)), 
           prior_inf_bin = if_else(prior_inf > 0, "1+", "0"))
  
  final
}

# custom function for unique infections - sequential approach ----

create_seq <- function(dataset, var_infection, var_attribution) {
  
  dia_flag <- 
    dataset %>% 
    filter(!is.na({{ var_infection }})) %>% 
    mutate(dia_cause = case_when(
      stooltype=="D1" & {{ var_attribution }} == 1 ~ "Attributable", 
      stooltype=="D1" & {{ var_attribution }} == 0 & {{ var_infection }} == 1 ~ "Unrelated detection", 
      stooltype=="D1" & {{ var_infection }}==0 ~ "No detection",
      TRUE ~ "No diarrhea")) %>% 
    group_by(pid) %>% 
    arrange(agedays, .by_group=T) %>% 
    mutate(grp=cumsum(dia_cause=="Attributable")) %>% 
    group_by(pid, grp) %>% 
    mutate(dia_start=if_else(any(dia_cause=="Attributable"), min(agestart, na.rm=T),NA_real_), 
           days_post_dia=agedays-dia_start) %>% 
    ungroup()
  
  # create episode groupings 
  det_grp <- 
    dia_flag %>% 
    group_by(pid) %>% 
    arrange(agedays, .by_group=T) %>%
    mutate(
      inc = if_else({{ var_infection }} == 1 & dia_cause != "Attributable" & 
                      lead(dia_cause) == "Attributable" & lead(agedays)-agedays <= 5, 1, 0),
      path_det_brk = case_when(month_ss-lag(month_ss) > 1 ~ 1, 
                               inc == 1 ~ 1, 
                               dia_cause=="Attributable" & 
                                 (lag(inc) == 0 | is.na(lag(inc))) ~ 1,
                               lag({{ var_infection }}) != {{ var_infection }} ~ 1, 
                               TRUE ~ 0), 
      path_det_grp = cumsum(path_det_brk)) %>% 
    group_by(pid, path_det_grp) %>% 
    mutate(path_det_type = case_when(
      dia_cause=="Attributable" ~ "Diarrhea", 
      inc==1 ~ "Incubation", 
      any(dia_cause=="Attributable") & dia_cause != "Attributable" ~ "Shedding",
      all({{ var_infection }}==0) ~ "Not detected", 
      all({{ var_infection }}==1 & dia_cause != "Attributable") ~ "Sub-clinical")) %>% 
    ungroup() 
  
  # time on study captured by each observation
  time_cap <- 
    det_grp %>% 
    select(pid, sid, country_id, region, 
           month_ss, agedays, agegroup3, agegroup6, agestart, ageend, 
           path_det_type, path_det_grp, score_binary) %>% 
    group_by(pid) %>%
    arrange(agedays, .by_group=T) %>%
    mutate(midday_pre = case_when(row_number()==1 & month_ss == 1 ~ 0, 
                                  row_number()==1 & month_ss > 1 ~ agedays-15, 
                                  month_ss-lag(month_ss) <= 1 ~ agedays-((agedays-lag(agedays))/2),
                                  TRUE ~ agedays-15), 
           midday_post = case_when(lead(month_ss)-month_ss <= 1 ~ agedays+((lead(agedays)-agedays)/2), 
                                   TRUE ~ agedays+15), 
           first_agestart = if_else(path_det_type=="Diarrhea", agestart, NA_real_), 
           last_ageend = if_else(path_det_type=="Diarrhea", ageend, NA_real_)) %>% 
    ungroup()
  
  # estimating start and end of infections
  episode_dur <- 
    time_cap %>% 
    group_by(pid,path_det_grp) %>%
    filter(all(path_det_type != "Not detected")) %>%
    mutate(first_age = min(midday_pre),
           last_age = max(midday_post),
           inf_type = if_else(all(path_det_type == "Sub-clinical"),
                              "Sub-clinical infection","Diarrhea (& shedding)"), 
           inf=1, 
           first_sid=case_when(
             path_det_type=="Diarrhea" ~ sid,
             row_number()==1 & !any(path_det_type=="Diarrhea") ~ sid, 
             TRUE ~ NA_character_), 
           max_sev=case_when(
             path_det_type=="Diarrhea" ~ score_binary+1, 
             inf_type == "Sub-clinical infection" ~ 0)) %>% 
    select(pid, sid, inf, first_sid, max_sev, first_age, last_age, path_det_grp, path_det_type) %>%
    ungroup()
  
  merged <- 
    left_join(x=time_cap, y=episode_dur) %>% 
    group_by(pid) %>% 
    arrange(agedays, .by_group=T) %>%
    mutate(inf = if_else(is.na(inf), 0, inf), 
           first = if_else(!is.na(first_sid), 1, 0),
           prior_inf = cumsum(first)-inf, 
           prior_inf_cat = if_else(prior_inf > 1, "2+", as.character(prior_inf)), 
           prior_inf_bin = if_else(prior_inf > 0, "1+", "0")) %>% 
    select(pid, sid, country_id, region, month_ss, agedays, agegroup3, agegroup6, inf, first_sid, path_det_type, path_det_grp, max_sev, first_age, last_age, first_agestart, last_ageend, midday_pre, midday_post, prior_inf, prior_inf_cat, prior_inf_bin)
  
  merged
}  

# create adenovirus analytic datasets ----
adeno_crude_alg <- create_crude(tac, inf_adeno, attr_alg_adeno)
write_csv(adeno_crude_alg, file=here::here("Outputs", "Data", "adeno-crude-alg-sens.csv"))

adeno_seq_alg <- create_seq(tac, inf_adeno, attr_alg_adeno)
write_csv(adeno_seq_alg, file=here::here("Outputs", "Data", "adeno-seq-alg-sens.csv"))

# create astrovirus analytic datasets ----
astro_crude_alg <- create_crude(tac, inf_astro, attr_alg_astro)
write_csv(astro_crude_alg, file=here::here("Outputs", "Data", "astro-crude-alg-sens.csv"))

astro_seq_alg <- create_seq(tac, inf_astro, attr_alg_astro)
write_csv(astro_seq_alg, file=here::here("Outputs", "Data", "astro-seq-alg-sens.csv"))

# create campylobacter analytic datasets ---- 
campy_crude_alg <- create_crude(tac, inf_campy, attr_alg_campy)
write_csv(campy_crude_alg, file=here::here("Outputs", "Data", "campy-crude-alg-sens.csv"))

campy_seq_alg <- create_seq(tac, inf_campy, attr_alg_campy)
write_csv(campy_seq_alg, file=here::here("Outputs", "Data", "campy-seq-alg-sens.csv"))

# create cryptosporidium analytic datasets ----
crypto_crude_alg <- create_crude(tac, inf_crypto, attr_alg_crypto)
write_csv(crypto_crude_alg, file=here::here("Outputs", "Data", "crypto-crude-alg-sens.csv"))

crypto_seq_alg <- create_seq(tac, inf_crypto, attr_alg_crypto)
write_csv(crypto_seq_alg, file=here::here("Outputs", "Data", "crypto-seq-alg-sens.csv"))

# create norovirus G.II analytic datasets ----
noro2_crude_alg <- create_crude(tac, inf_noro2, attr_alg_noro2)
write_csv(noro2_crude_alg, file=here::here("Outputs", "Data", "noro2-crude-alg-sens.csv"))

noro2_seq_alg <- create_seq(tac, inf_noro2, attr_alg_noro2)
write_csv(noro2_seq_alg, file=here::here("Outputs", "Data", "noro2-seq-alg-sens.csv"))

# create rotavirus analytic datasets ----
rota_crude_alg <- create_crude(tac, inf_rota, attr_alg_rota)
write_csv(rota_crude_alg, file=here::here("Outputs", "Data", "rota-crude-alg-sens.csv"))

rota_seq_alg <- create_seq(tac, inf_rota, attr_alg_rota)
write_csv(rota_seq_alg, file=here::here("Outputs", "Data", "rota-seq-alg-sens.csv"))

# create sapovirus analytic datasets ----
sapo_crude_alg <- create_crude(tac, inf_sapo, attr_alg_sapo)
write_csv(sapo_crude_alg, file=here::here("Outputs", "Data", "sapo-crude-alg-sens.csv"))

sapo_seq_alg <- create_seq(tac, inf_sapo, attr_alg_sapo)
write_csv(sapo_seq_alg, file=here::here("Outputs", "Data", "sapo-seq-alg-sens.csv"))

# create Shigella analytic datasets ----
shig_crude_alg <- create_crude(tac, inf_shig, attr_alg_shig)
write_csv(shig_crude_alg, file=here::here("Outputs", "Data", "shig-crude-alg-sens.csv"))

shig_seq_alg <- create_seq(tac, inf_shig, attr_alg_shig)
write_csv(shig_seq_alg, file=here::here("Outputs", "Data", "shig-seq-alg-sens.csv"))

# create ST-ETEC analytic datasets ----
stetec_crude_alg <- create_crude(tac, inf_stetec, attr_alg_stetec)
write_csv(stetec_crude_alg, file=here::here("Outputs", "Data", "stetec-crude-alg-sens.csv"))

stetec_seq_alg <- create_seq(tac, inf_stetec, attr_alg_stetec)
write_csv(stetec_seq_alg, file=here::here("Outputs", "Data", "stetec-seq-alg-sens.csv"))

# create typical EPEC analytic datasets ----
tepec_crude_alg <- create_crude(tac, inf_tepec, attr_alg_tepec)
write_csv(tepec_crude_alg, file=here::here("Outputs", "Data", "tepec-crude-alg-sens.csv"))

tepec_seq_alg <- create_seq(tac, inf_tepec, attr_alg_tepec)
write_csv(tepec_seq_alg, file=here::here("Outputs", "Data", "tepec-seq-alg-sens.csv"))
