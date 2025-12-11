library(tidyverse)
library(ggpubr)

ci_bound <- function(p, se, bound) {
  
  if(bound=="ll") {
    calc <- round(p-(1.96*se), digits=1)
  }
  
  if(bound=="ul") {
    calc <- round(p+(1.96*se), digits=1)
  }
  
  return(calc)
}

files <- list.files("../../Processed/Shigella/", full.names = T) %>% as.data.frame()
main24 <- files %>% 
  filter(grepl("Equal", .)==T | grepl("Worse", .) == T) %>% 
  mutate(inf_imm = if_else(grepl("Equal", .) == T, "Equal", "Worse"), 
         rate = case_when(grepl("Rate0", .) == T ~ "Low", 
                          grepl("Rate1", .) == T ~ "Medium", 
                          TRUE ~ "High"), 
         hybrid = if_else(grepl("Hybrid0", .) == T, "No", "Yes"), 
         ves = if_else(grepl("VE0", .) == T, "No", "Yes"))

main24_df <- data.frame()

for (i in 1:nrow(main24)) {
  
  temp <- read_csv(main24$.[i]) %>% 
    mutate(ve = (1-exp(beta))*100, ll = (1-exp(upper))*100, ul = (1-exp(lower))*100)
  
  fmt_temp <- 
    temp %>% 
    mutate(approach = case_when(
      dat %in% c(240, 120, 100, 20) & ana == 1 ~ 1, 
      dat %in% c(240, 120, 100, 20) & ana == 2 ~ 3,
      dat %in% c(24, 12, 80, 10) & ana == 1 ~ 2, 
      dat %in% c(24, 12, 80, 10) & ana == 2 ~ 4, 
      dat %in% c(240, 120, 100, 20) & ana == 0 ~ 5,
      dat %in% c(24, 12, 80, 10) & ana == 0 ~ 6), 
      diff_ve = if_else(outcome==2, ve-40, ve-60), 
      diff_sq = diff_ve^2,
      contains = case_when(
        outcome == 2 & ll <= 40 & 40 <= ul ~ 1, 
        outcome == 3 & ll <= 60 & 60 <= ul ~ 1, 
        TRUE ~ 0), 
      outside = if_else(contains == 0, 1, 0), 
      reject = if_else(ll <= 0 & 0 <= ul, 1, 0))
  
  summary_temp <- 
  fmt_temp %>% 
    mutate(truth=if_else(outcome==2, 40, 60)) %>%
    group_by(outcome, approach, ana, dat) %>% 
    summarise(
      n=n(),
      mean_ve = mean(ve),
      min_ve = min(ve), 
      max_ve = max(ve),
      mean_size=mean(size, na.rm=T), 
      min_size=min(size, na.rm=T), 
      max_size=max(size, na.rm=T),
      bias_p = mean(ve)-mean(truth),
      bias_mcse = sqrt(sum((ve-mean(ve))^2)/(500*499)),
      empse_p = sqrt(sum((ve-mean(ve))^2)/(499)), 
      mse_p = mean(diff_sq), 
      mse_mcse = sqrt(sum((diff_sq-mse_p)^2)/(500*499)), 
      coverage_p = mean(contains), 
      rejrate_p = mean(reject)) %>% 
    mutate(
      empse_mcse = empse_p/sqrt(2*499),
      coverage_mcse = sqrt((coverage_p*(1-coverage_p))/500), 
      rejrate_mcse = sqrt((rejrate_p*(1-rejrate_p))/500), 
      
      fmt_approach = factor(approach, levels=c(1:6), 
                            labels=c("Single Outcome,\nActive Surveillance", "Single Outcome,\nSymptom-Based Reporting", 
                                     "Stratified Recurrent Outcome,\nActive Surveillance", "Stratified Recurrent Outcome,\nSymptom-Based Reporting", 
                                     "Crude Recurrent Outcome,\nActive Surveillance", "Crude Recurrent Outcome,\nSymptom-Based Reporting"), ordered=T), 
      fmt_panel = factor(dat, levels=c(100, 80, 120, 12, 240, 24), 
                         labels=c("All Follow-Up", "All Follow-Up", 
                                  "First Year of Follow-Up", "First Year of Follow-Up", 
                                  "Second Year of Follow-Up", "Second Year of Follow-Up"), ordered=T),
      inf_imm = main24$inf_imm[i], 
      rate = main24$rate[i], 
      hybrid = main24$hybrid[i], 
      ves = main24$ves[i], 
      scenario = i, 
      bias_ll = bias_p-(1.96*bias_mcse), bias_ul = bias_p+(1.96*bias_mcse), 
      mse_ll = mse_p-(1.96*mse_mcse), mse_ul = mse_p+(1.96*mse_mcse), 
      coverage_ll = coverage_p-(1.96*coverage_mcse), coverage_ul = coverage_p+(1.96*coverage_mcse), 
      rejrate_ll = rejrate_p-(1.96*rejrate_mcse), rejrate_ul = rejrate_p+(1.96*rejrate_mcse))
  
  main24_df <- bind_rows(main24_df, summary_temp)
  
}

main24_df_long <- 
  main24_df %>% 
  mutate(ves = factor(ves, levels=c("No", "Yes"), ordered=T), 
         inf_imm = factor(inf_imm, levels=c("Equal", "Worse"), ordered=T), 
         hybrid = factor(hybrid, levels=c("No", "Yes"), ordered=T), 
         rate = factor(rate, levels=c("Low", "Medium", "High"), ordered=T), 
         fmt_approach = factor(approach, levels=c(1:6), 
                               labels=c("Single Outcome,\nActive Surveillance", "Single Outcome,\nSymptom-Based Reporting", 
                                        "Stratified Recurrent Outcome,\nActive Surveillance", "Stratified Recurrent Outcome,\nSymptom-Based Reporting", 
                                        "Crude Recurrent Outcome,\nActive Surveillance", "Crude Recurrent Outcome,\nSymptom-Based Reporting"), ordered=T)) %>% 
  pivot_longer(cols=c(ends_with("_p"), ends_with("_mcse")), names_to=c("metric", "type"), names_sep="_") %>% 
  pivot_wider(names_from="type", values_from="value") %>% 
  mutate(fmt_metric = factor(metric, levels=c("bias", "empse", "mse", "coverage", "rejrate"), 
                             labels=c("Bias", "Empirical Squared Error", "Mean Squared Error", "Coverage", "False Negative")))

fmt_main24_df <- 
  main24_df %>% 
  mutate(ves = factor(ves, levels=c("No", "Yes"), ordered=T), 
         inf_imm = factor(inf_imm, levels=c("Equal", "Worse"), ordered=T), 
         hybrid = factor(hybrid, levels=c("No", "Yes"), ordered=T), 
         rate = factor(rate, levels=c("Low", "Medium", "High"), ordered=T), 
         fmt_size = paste0(round(mean_size), "<br>(", 
                           round(min_size), ", ", round(max_size), ")"),
         fmt_ve = paste0(round(mean_ve, digits=1), "%<br>(", 
                         round(min_ve, digits=1), "%, ", round(max_ve, digits=1), "%)"), 
         fmt_bias = paste0(round(bias_p, digits=1), "<br>(", 
                           round(bias_ll, digits=1), ", ", round(bias_ul, digits=1), ")"), 
         fmt_mse = paste0(round(mse_p, digits=1), "<br>(", 
                           round(mse_ll, digits=1), ", ", round(mse_ul, digits=1), ")"), 
         fmt_cov = paste0(round(coverage_p*100, digits=1), "%<br>(", 
                           round(coverage_ll*100, digits=1), "%, ", round(coverage_ul*100, digits=1), "%)"), 
         fmt_rej = paste0(round(rejrate_p*100, digits=1), "%<br>(", 
                           round(rejrate_ll*100, digits=1), "%, ", round(rejrate_ul*100, digits=1), "%)")) %>% 
  select(outcome, dat, fmt_approach, ves, inf_imm, hybrid, rate, fmt_size, fmt_ve, fmt_bias, fmt_mse, fmt_cov, fmt_rej)

# large trials, severe VE ----

# figures for bias from highly powered trials, VE against severe Shigella diarrhea
main24_df_long %>% 
  filter(outcome == 3, dat %in% c(12, 24, 80, 100, 120, 240), metric=="bias", hybrid=="No", inf_imm=="Equal", ves=="No") %>% 
  ggplot(aes(x=fmt_approach, y=p, ymin=p-(1.96*mcse), ymax=p+(1.96*mcse), color=fmt_approach, shape=rate)) + 
  geom_point(size=2.25, position = position_dodge(width=0.3)) + 
  geom_errorbar(width=0.2, position = position_dodge(width=0.3)) + 
  scale_color_manual(values=c("#419BF5","#08519c", "#FB6725", "#a63603", "#777777", "#252525"), 
                     name=element_blank(), labels=element_blank(), guide="none") + 
  scale_shape_discrete(name="Burden at Site") + scale_y_continuous(name="Bias (Estimate - True VE)") + 
  facet_wrap(.~fmt_panel, ncol=1) + 
  # facet_wrap(.~fmt_panel, ncol=1, scales="free_y") # use to zoom y-axes on upper and middle panels
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.box = "vertical",
        axis.text.x=element_text(angle=45, hjust=1)) + 
  labs(x="")

ggsave("../../Figures/Shigella/Figure1.tiff", device="tiff", width=4.5, height=9, units="in")

main24_df_long %>% 
  filter(outcome == 3, dat %in% c(80, 100), metric=="bias", rate=="Medium") %>% 
  mutate(imm_scen = case_when(
    hybrid=="No" & ves == "No" ~ 1, 
    hybrid=="No" & ves == "Yes" ~ 2, 
    hybrid=="Yes" & ves == "No" ~ 3, 
    hybrid=="Yes" & ves == "Yes" ~ 4), 
    imm_scen = factor(imm_scen, levels=c(1:4), labels=c("No hybrid immunity,\n no vaccine effect on infection hazard", 
                                                        "No hybrid immunity,\n vaccine reduces infection hazard", 
                                                        "Hybrid immunity possible,\n no vaccine effect on infection hazard", 
                                                        "Hybrid immunity possible,\n vaccine reduces infection hazard"), ordered=T), 
    inf_imm = factor(inf_imm, levels = c("Worse", "Equal"), labels=c("Infection Immunity < Vaccine", 
                                                                     "Infection Immunity = Vaccine"), ordered=T)) %>%
  ggplot(aes(x=fmt_approach, y=p, ymin=p-(1.96*mcse), ymax=p+(1.96*mcse), color=fmt_approach, shape=imm_scen)) + 
  geom_point(size=1.75, position = position_dodge(width=0.35)) + 
  geom_errorbar(width=0.1, position = position_dodge(width=0.35)) + 
  scale_color_manual(values=c("#419BF5","#08519c", "#FB6725", "#a63603", "#777777", "#252525"), 
                     name=element_blank(), labels=element_blank(), guide="none") + 
  scale_shape_manual(name="", values=c(0,15,1,16)) + 
  scale_y_continuous(name="Bias (Estimate - True VE)") + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.box = "vertical",
        axis.text.x=element_text(angle=45, hjust=1)) + 
  labs(x="") + 
  facet_wrap(~inf_imm) + 
  guides(shape=guide_legend(nrow=2,byrow=TRUE))

ggsave("../../Figures/Shigella/Figure3.tiff", device="tiff", width=6.5, height=5, units="in")

# tables for highly powered trials, VE against severe Shigella diarrhea

fmt_main24_df %>% 
  filter(fmt_approach == "Single Outcome,\nActive Surveillance", outcome == 3, dat==100) %>% 
  arrange(ves, inf_imm, hybrid, rate) %>% 
  sjPlot::tab_df()

fmt_main24_df %>% 
  filter(fmt_approach == "Single Outcome,\nSymptom-Based Reporting", outcome == 3, dat==80) %>% 
  arrange(ves, inf_imm, hybrid, rate) %>% 
  sjPlot::tab_df()

fmt_main24_df %>% 
  filter(fmt_approach == "Stratified Recurrent Outcome,\nActive Surveillance", outcome == 3, dat==100) %>% 
  arrange(ves, inf_imm, hybrid, rate) %>% 
  sjPlot::tab_df()

fmt_main24_df %>% 
  filter(fmt_approach == "Stratified Recurrent Outcome,\nSymptom-Based Reporting", outcome == 3, dat==80) %>% 
  arrange(ves, inf_imm, hybrid, rate) %>% 
  sjPlot::tab_df()

fmt_main24_df %>% 
  filter(fmt_approach == "Crude Recurrent Outcome,\nActive Surveillance", outcome == 3, dat==100) %>% 
  arrange(ves, inf_imm, hybrid, rate) %>% 
  sjPlot::tab_df()

fmt_main24_df %>% 
  filter(fmt_approach == "Crude Recurrent Outcome,\nSymptom-Based Reporting", outcome == 3, dat==80) %>% 
  arrange(ves, inf_imm, hybrid, rate) %>% 
  sjPlot::tab_df()

fmt_main24_df %>% 
  filter(dat %in% c(12, 24, 120, 240), outcome == 3) %>% 
  arrange(fmt_approach, rate, dat) %>% 
  select(fmt_approach, rate, dat, starts_with("fmt_")) %>% 
  sjPlot::tab_df()

# realistic trials, severe VE ----

# figures for realistically sized trials in example scenarios
main24_df_long %>% 
  filter(outcome == 3, dat <=20, metric %in% c("bias", "coverage", "rejrate"), hybrid=="No", inf_imm=="Equal", ves=="No") %>% 
  ggplot(aes(x=fmt_approach, y=p, ymin=p-(1.96*mcse), ymax=p+(1.96*mcse), color=fmt_approach, shape=rate)) + 
  geom_point(size=2, position = position_dodge(width=0.3)) + 
  geom_errorbar(width=0.1, position = position_dodge(width=0.3)) + 
  scale_color_manual(values=c("#419BF5","#08519c", "#FB6725", "#a63603", "#777777", "#252525"), 
                     name=element_blank(), labels=element_blank(), guide="none") + 
  scale_shape_discrete(name="Burden at Site") + scale_y_continuous(name="Bias (True VE - Estimate)") + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.box = "vertical",
        axis.text.x=element_text(angle=45, hjust=1)) + 
  labs(x="") + 
  facet_wrap(.~metric, scales="free_y")

f2df <- 
bind_rows(
  read_csv("../../Processed/Shigella/Shig-Equal-Rate0-VE0-Hybrid0-VE-Results.csv") %>% mutate(rate="Low"),
  read_csv("../../Processed/Shigella/Shig-Equal-Rate1-VE0-Hybrid0-VE-Results.csv") %>% mutate(rate="Medium"),
  read_csv("../../Processed/Shigella/Shig-Equal-Rate2-VE0-Hybrid0-VE-Results.csv") %>% mutate(rate="High")) %>% 
  filter(outcome==3, dat %in% c(10,20,80,100)) %>% 
  mutate(ve = (1-exp(beta))*100, ll = (1-exp(upper))*100, ul = (1-exp(lower))*100, 
         approach = case_when(
           (dat == 100 | dat == 20) & ana == 1 ~ 1, 
           (dat == 100 | dat == 20) & ana == 2 ~ 3,
           (dat == 80 | dat == 10) & ana == 1 ~ 2, 
           (dat == 80 | dat == 10) & ana == 2 ~ 4, 
           (dat == 100 | dat == 20) & ana == 0 ~ 5,
           (dat == 80 | dat == 10) & ana == 0 ~ 6), 
         diff_ve = ve-60, 
         contains = if_else(ll <= 60 & 60 <= ul, 1, 0), 
         outside = if_else(contains == 0, 1, 0), 
         reject = if_else(ll <= 0 & 0 <= ul, 1, 0), 
         conc = case_when(
           reject== 0 ~ "Vaccine has some effect", 
           reject == 1 ~ "Vaccine has no effect"))

f2df %>% 
  mutate(rate = factor(rate, levels=c("Low", "Medium", "High"), ordered=T), 
         approach = factor(approach, levels=c(1:6), 
                               labels=c("Single Outcome,\nActive Surveillance", "Single Outcome,\nSymptom-Based\nReporting", 
                                        "Stratified Recurrent\nOutcome,\nActive Surveillance", "Stratified Recurrent\nOutcome, Symptom-\nBased Reporting", 
                                        "Crude Recurrent\nOutcome, Active\nSurveillance", "Crude Recurrent\nOutcome, Symptom-\nBased Reporting"), ordered=T)) %>% 
  group_by(rate, approach) %>% 
  arrange(diff_ve, reject, contains, .by_group = T) %>% 
  mutate(posvar=row_number(), conc = factor(conc, levels=c("Vaccine has some effect", "Vaccine has no effect"), ordered=T)) %>% 
  ggplot(aes(x=ve, xmin=ll, xmax=ul, y=posvar, color=conc, fill=conc)) + 
  geom_vline(xintercept=0, color="darkgrey") + 
  geom_errorbarh() + 
  facet_grid(rows = vars(rate), cols=vars(approach), switch="y") + 
  geom_vline(xintercept=60, color="black") + 
  scale_y_continuous(name="") + 
  scale_x_continuous(name="Estimated Vaccine Efficacy Against Severe Shigella Diarrhea") + 
  scale_color_brewer(name="Interpretation", type="qual") + scale_fill_discrete(name="Interpretation") + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(), 
        strip.text.x=element_text(size=9)) + 
  coord_cartesian(xlim=c(-25,100))

ggsave("../../Figures/Shigella/Figure2.tiff", device="tiff", width=8.9, height=5, units="in")

# tables for realistically sized trials, VE against severe Shigella diarrhea

fmt_main24_df %>% 
  filter(fmt_approach == "Single Outcome,\nActive Surveillance", outcome == 3, dat==20) %>% 
  arrange(ves, inf_imm, hybrid, rate) %>% 
  sjPlot::tab_df()

fmt_main24_df %>% 
  filter(fmt_approach == "Single Outcome,\nSymptom-Based Reporting", outcome == 3, dat==10) %>% 
  arrange(ves, inf_imm, hybrid, rate) %>% 
  sjPlot::tab_df()

fmt_main24_df %>% 
  filter(fmt_approach == "Stratified Recurrent Outcome,\nActive Surveillance", outcome == 3, dat==20) %>% 
  arrange(ves, inf_imm, hybrid, rate) %>% 
  sjPlot::tab_df()

fmt_main24_df %>% 
  filter(fmt_approach == "Stratified Recurrent Outcome,\nSymptom-Based Reporting", outcome == 3, dat==10) %>% 
  arrange(ves, inf_imm, hybrid, rate) %>% 
  sjPlot::tab_df()

fmt_main24_df %>% 
  filter(fmt_approach == "Crude Recurrent Outcome,\nActive Surveillance", outcome == 3, dat==20) %>% 
  arrange(ves, inf_imm, hybrid, rate) %>% 
  sjPlot::tab_df()

fmt_main24_df %>% 
  filter(fmt_approach == "Crude Recurrent Outcome,\nSymptom-Based Reporting", outcome == 3, dat==10) %>% 
  arrange(ves, inf_imm, hybrid, rate) %>% 
  sjPlot::tab_df()

# changes to VE, fixed scenario ----
f4df <- bind_rows(
  read_csv(files$.[13]) %>% mutate(truth=10),
  read_csv(files$.[14]) %>% mutate(truth=20),
  read_csv(files$.[15]) %>% mutate(truth=30), 
  read_csv(files$.[16]) %>% mutate(truth=40),
  read_csv(files$.[17]) %>% mutate(truth=50), 
  read_csv(files$.[18]) %>% mutate(truth=60),
  read_csv(files$.[19]) %>% mutate(truth=70),
  read_csv(files$.[20]) %>% mutate(truth=80),
  read_csv(files$.[21]) %>% mutate(truth=90)) %>% 
  mutate(ve = (1-exp(beta))*100, ll = (1-exp(upper))*100, ul = (1-exp(lower))*100) 

fmt_f4df <- 
  f4df %>% 
  mutate(approach = case_when(
    (dat == 100 | dat == 20) & ana == 1 ~ 1, 
    (dat == 100 | dat == 20) & ana == 2 ~ 3,
    (dat == 80 | dat == 10) & ana == 1 ~ 2, 
    (dat == 80 | dat == 10) & ana == 2 ~ 4, 
    (dat == 100 | dat == 20) & ana == 0 ~ 5,
    (dat == 80 | dat == 10) & ana == 0 ~ 6), 
    diff_ve = if_else(outcome==2, ve-40, ve-60), 
    diff_sq = diff_ve^2,
    contains = case_when(
      outcome == 2 & ll <= 40 & 40 <= ul ~ 1, 
      outcome == 3 & ll <= 60 & 60 <= ul ~ 1, 
      TRUE ~ 0), 
    outside = if_else(contains == 0, 1, 0), 
    reject = if_else(ll <= 0 & 0 <= ul, 1, 0)) %>% 
  group_by(truth, outcome, approach, ana, dat) %>% 
  summarise(
    n=n(),
    bias_p = mean(ve)-mean(truth),
    bias_mcse = sqrt(sum((ve-mean(ve))^2)/(500*499)),
    empse_p = sqrt(sum((ve-mean(ve))^2)/(499)), 
    mse_p = mean(diff_sq), 
    mse_mcse = sqrt(sum((diff_sq-mse_p)^2)/(500*499)), 
    coverage_p = mean(contains), 
    rejrate_p = mean(reject)) %>% 
  mutate(
    empse_mcse = empse_p/sqrt(2*499),
    coverage_mcse = sqrt((coverage_p*(1-coverage_p))/500), 
    rejrate_mcse = sqrt((rejrate_p*(1-rejrate_p))/500), 
    fmt_approach = factor(approach, levels=c(1:6), 
                          labels=c("Single Outcome,\nActive Surveillance", "Single Outcome,\nSymptom-Based Reporting", 
                                   "Stratified Recurrent Outcome,\nActive Surveillance", "Stratified Recurrent Outcome,\nSymptom-Based Reporting", 
                                   "Crude Recurrent Outcome,\nActive Surveillance", "Crude Recurrent Outcome,\nSymptom-Based Reporting"), ordered=T)) %>% 
  pivot_longer(cols=c(ends_with("_p"), ends_with("_mcse")), names_to=c("metric", "type"), names_sep="_") %>% 
  pivot_wider(names_from="type", values_from="value") %>% 
  mutate(fmt_metric = factor(metric, levels=c("bias", "empse", "mse", "coverage", "rejrate"), 
                             labels=c("Bias", "Empirical Squared Error", "Mean Squared Error", "Coverage", "False Negative")))

fmt_f4df %>% 
  filter(outcome == 2, dat >20, metric=="bias") %>% 
  ggplot(aes(x=truth, y=p, ymin=p-(1.96*mcse), ymax=p+(1.96*mcse), color=fmt_approach, shape=fmt_approach)) + 
  geom_point(size=1, position = position_dodge(width=0.35)) + geom_line() + 
  geom_errorbar(width=0.1, position = position_dodge(width=0.35)) + 
  scale_shape_discrete(name="") + 
  scale_color_manual(values=c("#419BF5","#08519c", "#FB6725", "#a63603", "#777777", "#252525"), name="") + 
  theme_bw() + ylab("Mean Bias and 95% Confidence Intervals") + 
  scale_x_continuous(name="True Efficacy Against Shigella Diarrhea", breaks=c(seq(0,100,by=10))) + 
  theme(legend.position = "bottom", 
        legend.box = "vertical",
        axis.text.x=element_text(angle=45, hjust=1)) + labs(x="") + 
  guides(shape=guide_legend(nrow=2,byrow=TRUE))

ggsave("../Figures/Shigella/Figure4.tiff", device="tiff", width=6.5, height=5, units="in")

# tables for any severity VE ----

# highly powered trials
fmt_main24_df %>% 
  filter(fmt_approach == "Single Outcome,\nActive Surveillance", outcome == 2, dat==100) %>% 
  arrange(ves, inf_imm, hybrid, rate) %>% 
  sjPlot::tab_df()

fmt_main24_df %>% 
  filter(fmt_approach == "Single Outcome,\nSymptom-Based Reporting", outcome == 2, dat==80) %>% 
  arrange(ves, inf_imm, hybrid, rate) %>% 
  sjPlot::tab_df()

fmt_main24_df %>% 
  filter(fmt_approach == "Stratified Recurrent Outcome,\nActive Surveillance", outcome == 2, dat==100) %>% 
  arrange(ves, inf_imm, hybrid, rate) %>% 
  sjPlot::tab_df()

fmt_main24_df %>% 
  filter(fmt_approach == "Stratified Recurrent Outcome,\nSymptom-Based Reporting", outcome == 2, dat==80) %>% 
  arrange(ves, inf_imm, hybrid, rate) %>% 
  sjPlot::tab_df()

fmt_main24_df %>% 
  filter(fmt_approach == "Crude Recurrent Outcome,\nActive Surveillance", outcome == 2, dat==100) %>% 
  arrange(ves, inf_imm, hybrid, rate) %>% 
  sjPlot::tab_df()

fmt_main24_df %>% 
  filter(fmt_approach == "Crude Recurrent Outcome,\nSymptom-Based Reporting", outcome == 2, dat==80) %>% 
  arrange(ves, inf_imm, hybrid, rate) %>% 
  sjPlot::tab_df()

# realistically sized trials
fmt_main24_df %>% 
  filter(fmt_approach == "Single Outcome,\nActive Surveillance", outcome == 2, dat==20) %>% 
  arrange(ves, inf_imm, hybrid, rate) %>% 
  sjPlot::tab_df()

fmt_main24_df %>% 
  filter(fmt_approach == "Single Outcome,\nSymptom-Based Reporting", outcome == 2, dat==10) %>% 
  arrange(ves, inf_imm, hybrid, rate) %>% 
  sjPlot::tab_df()

fmt_main24_df %>% 
  filter(fmt_approach == "Stratified Recurrent Outcome,\nActive Surveillance", outcome == 2, dat==20) %>% 
  arrange(ves, inf_imm, hybrid, rate) %>% 
  sjPlot::tab_df()

fmt_main24_df %>% 
  filter(fmt_approach == "Stratified Recurrent Outcome,\nSymptom-Based Reporting", outcome == 2, dat==10) %>% 
  arrange(ves, inf_imm, hybrid, rate) %>% 
  sjPlot::tab_df()

fmt_main24_df %>% 
  filter(fmt_approach == "Crude Recurrent Outcome,\nActive Surveillance", outcome == 2, dat==20) %>% 
  arrange(ves, inf_imm, hybrid, rate) %>% 
  sjPlot::tab_df()

fmt_main24_df %>% 
  filter(fmt_approach == "Crude Recurrent Outcome,\nSymptom-Based Reporting", outcome == 2, dat==10) %>% 
  arrange(ves, inf_imm, hybrid, rate) %>% 
  sjPlot::tab_df()
