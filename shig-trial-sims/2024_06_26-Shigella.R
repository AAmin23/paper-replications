library(tidyverse)
library(ggplot2)
library(gam)
library(survival)
library(doParallel)
library(foreach)

nocores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
#nocores <- detectCores()-1
cl <- makeCluster(nocores-1)
registerDoParallel(cl)

# overall conditions
n_sim <- 500
init_seed <- 121121
pathogen <- "Shigella"

# trial conditions
days_recruitment <- 365
n_participants <- 10000

# vaccine parameters
VEs_0 <- 0
VEs_20 <- 0.2
VEsp <- 0.40
VEsph <- 0.60

# age trends as pulled from paper; use consecutive detections and algorithm attribution
# infection rates using 3-month bands
inf <- data.frame(
  age = c("10-12", "13-15", "16-18", "19-21", "22-24"), 
  rate = c(6.03, 7.49, 9.62, 11.19, 10.30)) %>% 
  mutate(rate = rate/100/30.5, 
         rr_inf = rate/(11.19/100/30.5)) # use 19-21 as referent

# probability of diarrhea, given infection; 3-month bands for 10-12 and 6-month bands otherwise
dia <- data.frame(
  age = c("10-12", "13-15", "16-18", "19-21", "22-24"), 
  probability = c(0.36, 0.39, 0.39, 0.35, 0.35)) %>% 
  mutate(rr_dia = probability/0.35) # use 19-21 as referent

# probability of severe diarrhea, given diarrhea; 3-month bands for 10-12 and 6-month bands otherwise
sev <- data.frame(
  age = c("10-12", "13-15", "16-18", "19-21", "22-24"), 
  probability = c(0.125, 0.155, 0.155, 0.084, 0.084)) %>% 
  mutate(rr_sev = probability/0.084) # use 19-21 as referent

trends <- full_join(x=inf, y=dia, by="age") %>% 
  full_join(y=sev, by="age") %>% 
  select(age, starts_with("rr")) %>% 
  ## 3rd yr of life is hypothetical
  add_row(age=c("25-30", "31-36"), rr_inf=c(0.9,0.8), rr_dia=c(0.9,0.8), rr_sev=c(0.9, 0.8))

## last age group is always the censoring age - so no entries for end & RRs
age_info <- data.frame(
  age_group = c(1:8), 
  age = c("10-12", "13-15", "16-18", "19-21", "22-24", "25-30", "31-36", ">36"),
  age_start = c(91*3,366,366+91,366+(91*2),366+(91*3),(365*2)+1,(365*2)+1+(91*2),365*3), # last age should be the last age (in days) *included* in follow-up so censoring occurs appropriately
  age_end = c(365,365+91,365+(91*2),365+(91*3),365*2,(365*2)+(91*2),365*3,NA)) %>% 
  full_join(y=trends, by="age")

# natural history for 19-21 reference
p_symptoms <- 0.35
p_severe <- 0.084

# transmission scenarios - rate is for the 19-21 referent
# numeric rate is based on specific sites; requires conversion from x/100 infant-months to x/day
rate0 <- data.frame(calendar_int=1, day_start=1, day_end=1461, rate=2.18/100/30.5) # low rate; based on Brazil
rate1 <- data.frame(calendar_int=1, day_start=1, day_end=1461, rate=median(c(2.18/100/30.5, 11.57/100/30.5))) # medium rate; midpoint between the low and high
rate2 <- data.frame(calendar_int=1, day_start=1, day_end=1461, rate=11.57/100/30.5) # high rate; based on Bangladesh

source("./2024_05_02-Scenarios.R") # master script for all scenarios regardless of pathogen

