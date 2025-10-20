library(tidyverse)
library(lubridate)

# custom functions and operators ----
falseifNA <- function(x) {ifelse(is.na(x), FALSE, x)}
ifelse2 <- function(x, a, b) {ifelse(falseifNA(x), a, b)}
specify_decimal <- function(x, k) {trimws(format(round(x, k), nsmall=k))}
`%notin%` <- Negate(`%in%`)

# import the required datasets ----

# results from tac arrays and original attribution
tac <- read_csv(here::here("Inputs", "Tac_Sep2018.csv")) %>% 
  select(country_id, pid, sid, date, month_ss, agedays, stooltype, stooltype7, score, postrota,
         adeno_ct=adenovirus_40_41, adeno_afe=adenovirus_40_41_afe,
         aero_ct=aeromonas, aero_afe=aeromonas_afe,
         ancyclo_ct=ancyclostoma, ancyclo_afe=ancyclostoma_afe,
         astro_ct=astrovirus, astro_afe=astrovirus_afe, 
         campy_ct=campylobacter_jejuni_coli, campy_afe=campylobacter_jejuni_coli_afe,
         crypto_ct=cryptosporidium, crypto_afe=cryptosporidium_afe,
         cyclo_ct=cyclospora, cyclo_afe=cyclospora_afe,
         ebien_ct=e_bieneusi, ebien_afe=e_bieneusi_afe, 
         ehist_ct=e_histolytica, ehist_afe=e_histolytica_afe, 
         eint_ct=e_intestinalis, eint_afe=e_intestinalis_afe, 
         giar_ct=giardia, giar_afe=giardia_afe, 
         hpyl_ct=h_pylori, hpyl_afe=h_pylori_afe, 
         iso_ct=isospora, iso_afe=isospora_afe, 
         ples_ct=plesiomonas, ples_afe=plesiomonas_afe, 
         rota_ct=rotavirus, rota_afe=rotavirus_afe, 
         salm_ct=salmonella, salm_afe=salmonella_afe, 
         shig_ct=shigella_eiec, shig_afe=shigella_eiec_afe, 
         sapo_ct=sapovirus, sapo_afe=sapovirus_afe, 
         strong_ct=strongyloides, strong_afe=strongyloides_afe, 
         trich_ct=trichuris, trich_afe=trichuris_afe, 
         chol_ct=v_cholerae, chol_afe=v_cholerae_afe, 
         tepec_ct=tEPEC, tepec_afe=tEPEC_afe, 
         stec_ct=STEC, stec_afe=STEC_afe, 
         stetec_ct=ST_ETEC, stetec_afe=ST_ETEC_afe, 
         ltetec_ct=LT_ETEC, ltetec_afe=LT_ETEC_afe, 
         noro1_ct=norovirus_gi, noro1_afe=norovirus_gi_afe, 
         noro2_ct=norovirus_gii, noro2_afe=norovirus_gii_afe) %>% 
  mutate(date = mdy(date)) 

# start and end of diarrheal episodes + fix two errors
age <- read_csv(here::here("Inputs", "microtacages.csv")) %>% 
  mutate(agestart = ifelse(sid == "BG1002105",88,agestart),
         ageend = ifelse(sid == "BG1002105",89,ageend)) %>%
  mutate(agestart = ifelse(sid == "PK1005658",129,agestart),
         ageend = ifelse(sid == "PK1005658",131,ageend))

# create key vars in the tac dataset ----
specimens <- 
  inner_join(x=age, y=tac, by=c("pid", "sid", "agedays")) %>% 
  group_by(pid) %>% 
  arrange(date, agedays, .by_group = T) %>% 
  mutate(
    obs_n = row_number(), 
    days_post = if_else(obs_n > 1 & stooltype == "M1" & lag(stooltype) == "D1", agedays-lag(ageend), NA_real_), 
    days_pre = if_else(obs_n<max(obs_n) & stooltype == "M1" & lead(stooltype) == "D1", lead(agestart)-agedays, NA_real_),
    d1_location = case_when(
      days_post <= 7 & days_pre <=7  & days_post==days_pre ~ "Both",
      days_post <= 7 & days_pre <=7  & days_post<days_pre ~ "Lag",
      days_post <= 7 & days_pre <=7  & days_post>days_pre ~ "Lead",
      days_post <= 7 ~ "Lag", 
      days_pre <= 7 ~ "Lead", 
      is.na(stooltype7) ~ "Missing",
      TRUE ~ "N/A"),
      # obs_n == 1 & lead(days_since) <= 7 & lead(stooltype7) == "D1" ~ "Lead", 
      # obs_n == max(obs_n) & days_since <= 7 & lag(stooltype7) == "D1" ~ "Lag", 
      # lead(days_since) <= 7 & lead(stooltype7) == "D1" & days_since <= 7 & lag(stooltype7) == "D1" & lead(days_since) == days_since ~ "Both",
      # lead(days_since) <= 7 & lead(stooltype7) == "D1" & days_since <= 7 & lag(stooltype7) == "D1" & lead(days_since) > days_since ~ "Lag",
      # lead(days_since) <= 7 & lead(stooltype7) == "D1" & days_since <= 7 & lag(stooltype7) == "D1" & lead(days_since) < days_since ~ "Lead",
      # lead(days_since) <= 7 & lead(stooltype7) == "D1" ~ "Lead", 
      # days_since <= 7 & lag(stooltype7) == "D1" ~ "Lag", 
      # TRUE ~ "Missing"), 
    
    month_ss = if_else(month_ss == 0, 1, month_ss), 
    agegroup3 = cut(month_ss, breaks=seq(0,24,by=3), 
                    labels=c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")), 
    agegroup6 = cut(month_ss, breaks=seq(0,24,by=6), 
                    labels=c("1-6", "7-12", "13-18", "19-24")), 
    
    region = case_when(
      country_id %in% c("BG", "NP", "PK", "IN") ~ "SEA", 
      country_id %in% c("TZ", "SA") ~ "AFR", 
      TRUE ~ "SAM"), 
    
    rvvIntro = if_else(country_id %in% c("BR", "PE", "SA"), 1, 0),
    rota_ct = if_else(postrota > 0, NA_real_, rota_ct), 
    rota_afe = if_else(postrota > 0, NA_real_, rota_afe),
    
    score_binary = if_else(score > 6, 1, 0),
    
    across(.cols = ends_with("_ct"), 
           .fns = function(x) {if_else(x <= 30, 1, 0)}, # binary indicator for infection
           .names = "inf_{.col}"), 
    
    # replace certain thresholds based on Liz's protection from infection paper
    inf_adeno = if_else(adeno_ct <= 30.424, 1, 0),
    inf_noro2 = if_else(noro2_ct <= 30.357, 1, 0),
    inf_rota = if_else(rota_ct <= 32.638, 1, 0),
    inf_shig = if_else(shig_ct <= 30.507, 1, 0)) %>% 
  rename_with(.fn = function(x) {substr(x,1,nchar(x)-3)}, .cols=starts_with("inf_"))

# characterize all diarrheal episodes ----
diarrhea <- specimens %>% 
  select(pid, sid, agedays, agestart, ageend, country_id, month_ss, stooltype) %>% 
  arrange(pid, agedays) %>% 
  group_by(pid) %>% 
  mutate(grp=cumsum(stooltype=="D1")) %>% 
  group_by(pid, grp) %>% 
  mutate(diarrhea_start = ifelse2(any(stooltype == "D1"), min(agestart,na.rm = T),NA), 
         diarrhea_end = ifelse2(any(stooltype == "D1"), min(ageend,na.rm = T),NA),
         diarrhea_days = diarrhea_end - diarrhea_start +1, 
         days_post=agedays - diarrhea_start) %>%
  ungroup()

# create long datasets for ct and afe values for pathogens ----
ct <- diarrhea %>%
  select(sid,pid,agedays,agestart,stooltype) %>%
  group_by(pid) %>%
  mutate(lag.sid = lag(sid),
         lag2.sid = lag(sid, n=2)) %>%
  left_join(specimens %>% select(sid, ends_with("_ct")), by=c("pid","sid")) %>%
  pivot_longer(cols = ends_with("_ct"),
               names_to = "pathogen",
               values_to = "ct") %>%
  ungroup() %>%
  mutate(pathogen = str_sub(pathogen, end=-4))

afe <- diarrhea %>% 
  filter(stooltype=="D1") %>%
  select(sid,pid,agedays) %>% 
  left_join(specimens,by=c("pid","sid","agedays")) %>% 
  select(pid,sid,agedays,ends_with("_afe")) %>% 
  pivot_longer(cols=ends_with("_afe"), names_to="pathogen", values_to="afe") %>% 
  group_by(sid) %>% 
  filter(!is.na(afe)) %>%
  distinct() %>% 
  select(sid,pathogen,afe) %>% 
  mutate(pathogen=str_sub(pathogen,end=-5))

ct.afe <- ct %>% left_join(afe, by = c("sid","pathogen"))

# create variables for lagged specimen info ----

lag1 <- ct.afe %>% select(lag.sid=sid,lag.agedays=agedays,pathogen,lag.ct=ct)

lag2 <- ct.afe %>% select(lag2.sid=sid,lag2.agedays=agedays,pathogen,lag2.ct=ct) 

# pathogen attribution for diarrheal episodes ONLY ----
attribution <- ct.afe %>%
  left_join(lag1,by = c("lag.sid","pathogen")) %>%
  left_join(lag2,by = c("lag2.sid","pathogen")) %>%
  mutate(lag.days = agestart - lag.agedays,
         lag2.days = agestart - lag2.agedays) %>%
  select(pid,lag2.sid,lag.sid,sid,stooltype,
         pathogen,
         lag2.ct,lag.ct,ct,afe,
         lag2.days,lag.days,agestart,agedays) %>%
  ungroup() %>%
  filter(stooltype == "D1") %>%
  mutate(ctlag.diff = lag.ct-ct,
         ctlag2.diff = lag2.ct-ct) %>%
  group_by(sid) %>%
  # if the first specimen collected is diarrhea
  mutate(alg.step = ifelse(is.na(lag.sid) &
                               afe > 0.5, "Cause - 1st specimen", NA)) %>%
  mutate(alg.step = ifelse(is.na(lag.sid) &
                               all(is.na(alg.step)) &
                               afe == max(afe) &
                               afe > 0, "Cause - 1st specimen", alg.step)) %>%
  # all other specimens
  mutate(alg.step = ifelse(!is.na(lag.sid) &
                               all(is.na(alg.step)) &
                               afe > 0.5 &
                               ((lag.days <= 60 &
                                   ctlag.diff >= 3.322) |
                                  ((lag.days <= 3 | is.na(lag.ct)) &
                                     lag2.days <= 60 &
                                     ctlag2.diff >= 3.322)), "Cause - Step 1", alg.step)) %>%
  mutate(alg.step = ifelse(!is.na(lag.sid) &
                               all(is.na(alg.step)) &
                               afe > 0 &
                               ((lag.days > 3 &
                                   lag.days <= 60 &
                                   ctlag.diff >= 3.322) |
                                  ((is.na(lag.ct) | lag.days <= 3) &
                                     lag2.days <= 60 &
                                     ctlag2.diff >= 3.322)), "Cause - Step 2", alg.step)) %>%
  mutate(alg.step = ifelse(!is.na(lag.sid) &
                               all(is.na(alg.step)) &
                               afe > 0.5, "Cause - Step 3", alg.step)) %>%
  mutate(alg_cause = ifelse(is.na(alg.step),0,1),
         afe_cause = ifelse2(afe > 0.5, 1, 0)) %>%
  mutate(alg_cause = ifelse(is.na(ct),NA,alg_cause),
         afe_cause = ifelse(is.na(afe),NA,afe_cause))

# creating primary dataset for subsequent analyses ----

# whether each pathogen was causal based on new algorithm - wide
alg_attr <- attribution %>%
  select(pid,sid,stooltype,pathogen,alg_cause) %>%
  arrange(pathogen) %>%
  pivot_wider(names_from = "pathogen",values_from = "alg_cause", names_prefix = "attr_alg_") %>%
  rowwise() %>%
  mutate(alg_tot = sum(across(starts_with("attr_alg_")), na.rm = T)) %>%
  mutate(alg_multi = ifelse(alg_tot > 1,1,0))

# whether each pathogen was causal based on AFe alone - wide
afe_attr <- attribution %>%
  select(pid,sid,stooltype,pathogen,afe_cause) %>%
  arrange(pathogen) %>%
  pivot_wider(names_from = "pathogen",values_from = "afe_cause", names_prefix = "attr_afe_") %>%
  rowwise() %>%
  mutate(afe_tot = sum(across(starts_with("attr_afe_")), na.rm = T)) %>%
  mutate(afe_multi = ifelse(afe_tot > 1,1,0))

# adding this info back into main dataset
primary <- 
  specimens %>% 
  left_join(y=afe_attr, by=c("pid", "sid", "stooltype")) %>% 
  left_join(y=alg_attr, by=c("pid", "sid", "stooltype"))

write_csv(primary, file=here::here("Outputs", "Data", "3.0-path-attr-data.csv"))
