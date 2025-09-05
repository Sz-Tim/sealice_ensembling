# Project: Sealice IP Ensemble
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Results



# setup -------------------------------------------------------------------
library(tidyverse); library(glue)
library(sf)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
library(biotrackR) # devtools::install_github("Sz-Tim/biotrackR")
library(rstan)
library(yardstick)
library(terra)
library(recipes)
library(ggpubr)
library(cowplot)
library(scico)
library(ggdist)
library(ggnewscale)
theme_set(theme_bw() + theme(panel.grid=element_blank()))
source("code/00_fn.R")


palette_cols <- colorspace::diverge_hcl(8, palette="Vik")
modType_cols <- c("Ens['Blend']"=palette_cols[length(palette_cols)-1],
                  "Ens['Avg']"=palette_cols[length(palette_cols)-2],
                  "Constituent"=palette_cols[3])
lab_expressions <- c(expression(Ens['Blend']),
                     expression(Ens['Avg']),
                     expression(Constituent))
modType2_cols <- c("Ens['Blend']"=palette_cols[length(palette_cols)-1],
                   "Ens['Avg']"=palette_cols[length(palette_cols)-2],
                   "Opt['Param']"=palette_cols[2],
                   "Other"=palette_cols[3],
                   "Constituent"=palette_cols[3])
lab2_expressions <- c(expression(Ens['Blend']),
                      expression(Ens['Avg']),
                      expression(Opt['Param']),
                      expression(Other))
modType3_cols <- c("Ens['Fcst']"="black",
                   "Ens['Blend']"=palette_cols[length(palette_cols)-1],
                   "Ens['Avg']"=palette_cols[length(palette_cols)-2],
                   "Opt['Param']"=palette_cols[2],
                   "Other"=palette_cols[3],
                   "Constituent"=palette_cols[3],
                   "2D"=palette_cols[3],
                   "3D"=palette_cols[3])


# cmr <- readRDS("../00_misc/cmr_cmaps.RDS")

# Full dataset
ensFull_df <- read_csv("out/valid_df_2021-2024.csv") |>
  mutate(lice_g05=factor(licePerFish_rtrt^4 > 0.5))

site_i <- read_csv("data/farm_sites.csv") 
sim_i <- read_csv("out/sim_2021-2024/sim_i.csv") |>
  mutate(sim=paste0("sim_", i),
         lab_short=if_else(fixDepth, "2D", "3D")) |>
  group_by(lab_short) |>
  mutate(lab=paste0("'", lab_short, ".", row_number(), "'")) |>
  ungroup() |>
  select(sim, lab_short, lab) |>
  bind_rows(
    tibble(sim=c("predFcst", "predBlend", "sim_avgAll", 
                 "null0", "nullTime", "nullFarm"),
           lab_short=c("Ens['Fcst']", "Ens['Blend']", "Ens['Avg']", 
                       "Null[0]", "Null['time']", "Null['farm']"),
           lab=c("Ens['Fcst']", "Ens['Blend']", "Ens['Avg']",
                 "Null[0]", "Null['time']", "Null['farm']"))
  ) |>
  mutate(lab=factor(lab, 
                    levels=c("Ens['Fcst']", "Ens['Blend']", "Ens['Avg']", 
                             paste0("'3D.", 1:16, "'"), paste0("'2D.", 1:4, "'"),
                             "Null[0]", "Null['time']", "Null['farm']")),
         lab_short=factor(lab_short, 
                          levels=c("Ens['Fcst']", "Ens['Blend']", "Ens['Avg']", 
                                   "3D", "2D", 
                                   "Null[0]", "Null['time']", "Null['farm']")))



# simulation settings -----------------------------------------------------

var_pretty <- list("D_h"="Diffusion: horizontal", 
                   "D_hVert"="Diffusion: vertical",
                   "eggTemp_fn"="Gravid egg production", 
                   "mortSal_fn"="Larval mortality rate", 
                   "fixDepth"="Movement dimensionality",
                   "swimUpSpeedMean"="Upward swimming speed", 
                   "swimDownSpeedMean"="Downward swimming speed",
                   "salinityThreshMin"="Salinity: Lower threshold", 
                   "salinityThreshMax"="Salinity: Upper threshold", 
                   "lightThreshNauplius"="Light threshold: Nauplius", 
                   "lightThreshCopepodid"="Light threshold: Copepodid",
                   "viableDegreeDays"="Development: Nauplius") |>
  as_tibble() |>
  pivot_longer(everything(), names_to="var", values_to="var_pretty") |>
  mutate(var_order=factor(var, 
                          levels=c("fixDepth", "D_h", "D_hVert", 
                                   "eggTemp_fn", "mortSal_fn", "viableDegreeDays",
                                   "salinityThreshMin", "salinityThreshMax",
                                   "swimDownSpeedMean", "swimUpSpeedMean",
                                   "lightThreshNauplius", "lightThreshCopepodid"))) |>
  arrange(var_order) |>
  mutate(units=list("Spatial dynamics", expression(m^2 %.% s^-1), expression(m^2 %.% s^-1),
                 "Function", "Function", expression(degree*C %.% d),
                 "psu", "psu", 
                 expression(cm %.% s^-1), expression(cm %.% s^-1),
                 expression(mu*'mol' %.%~'m'^-2 %.% s^-1), expression(mu*'mol' %.%~'m'^-2 %.% s^-1)),
         varNumeric=c(0, 1, 1, 
                      0, 0, 1, 
                      1, 1, 
                      1, 1, 
                      1, 1),
         var_cor=c(0, 0, 0, 
                   0, 0, 0, 
                   1, 1, 
                   1, 1, 
                   1, 1),
         var_pairs=c(1, 2, 3,
                     4, 5, 6,
                     7, 7, 
                     8, 8, 
                     9, 9))

read_csv("out/sim_2021-2024/sim_i.csv") |>
  mutate(sim=paste0("sim_", i),
         lab_short=if_else(fixDepth, "2D", "3D")) |>
  group_by(lab_short) |>
  mutate(Simulation=paste0(lab_short, ".", row_number())) |>
  ungroup() |>
  arrange(desc(fixDepth), i) |>
  mutate(across(where(is.numeric), ~signif(.x, 3))) |>
  mutate(salinityThresh=paste0(salinityThreshMin, "-", salinityThreshMax),
         swimUpSpeedMean=100*swimUpSpeedMean,
         swimDownSpeedMean=100*swimDownSpeedMean,
         eggTemp_fn=if_else(eggTemp_fn=="constant", "28.2", "f(Temp.)"),
         mortSal_fn=if_else(mortSal_fn=="constant", "0.01", "f(Sal.)")) |>
  mutate(across(where(is.numeric), ~format(.x, scientific=T, digits=3))) |> 
  select(Simulation, D_h, D_hVert, eggTemp_fn, mortSal_fn, viableDegreeDays,
         swimUpSpeedMean, swimDownSpeedMean,
         salinityThresh, lightThreshNauplius, lightThreshCopepodid) |>
  write_csv("ms/table_1_new.csv")

n_total <- 50000
n_sim3D <- 0.8*n_total
n_sim2D <- 0.2*n_total

light_mx <- MASS::mvrnorm(n_sim3D, c(0,0), matrix(c(1, 0.8, 0.8, 1), nrow=2))
swim_mx <- MASS::mvrnorm(n_sim3D, c(0,0), matrix(c(1, 0.8, 0.8, 1), nrow=2))
par_dist_df <- bind_rows(
  tibble(fixDepth="false",
         D_hVert=exp(runif(n_sim3D, log(1e-5), log(1e-1))),
         mortSal_fn=sample(c("constant", "logistic"), n_sim3D, replace=T),
         eggTemp_fn=sample(c("constant", "logistic"), n_sim3D, replace=T),
         salinityThreshMin=runif(n_sim3D, 20, 28),
         salinityThreshMax=pmin(salinityThreshMin + runif(n_sim3D, 0.1, 6), 32),
         lightThreshCopepodid=qunif(pnorm(light_mx[,1]), (2e-6)^0.5, (2e-4)^0.5)^2,
         lightThreshNauplius=qunif(pnorm(light_mx[,2]), (0.05)^0.5, (0.5)^0.5)^2,
         swimUpSpeedMean=-(qunif(pnorm(swim_mx[,1]), (1e-4)^0.5, (2e-2)^0.5))^2,
         swimDownSpeedMean=(qunif(pnorm(swim_mx[,2]), (1e-4)^0.5, (2e-2)^0.5))^2),
  tibble(mortSal_fn=sample(c("constant", "logistic"), n_sim2D, replace=T),
         eggTemp_fn=sample(c("constant", "logistic"), n_sim2D, replace=T),
         fixDepth="true")
) |>
  # mutate(across(where(is.numeric), ~if_else(is.na(.x), 0, .x))) |>
  mutate(D_h=exp(runif(n_sim2D+n_sim3D, log(1e-3), log(1e1))),
         viableDegreeDays=runif(n_sim2D+n_sim3D, 30, 50),
         swimUpSpeedMean=100*swimUpSpeedMean,
         swimDownSpeedMean=100*swimDownSpeedMean) |> 
  mutate(fixDepth=if_else(fixDepth=="true", "2D", "3D"))

sim_used <- read_csv("out/sim_2021-2024/sim_i.csv") |>
  mutate(sim=paste0("sim_", i),
         lab_short=if_else(fixDepth, "2D", "3D")) |>
  group_by(lab_short) |>
  mutate(Simulation=paste0(lab_short, ".", row_number())) |>
  ungroup() |>
  arrange(desc(fixDepth), i) |>
  mutate(across(where(is.numeric), ~signif(.x, 3))) |>
  mutate(salinityThresh=paste0(salinityThreshMin, "-", salinityThreshMax),
         swimUpSpeedMean=100*swimUpSpeedMean,
         swimDownSpeedMean=100*swimDownSpeedMean) |>
  mutate(across(any_of(c("swimDownSpeedMean", "swimUpSpeedMean",
                         "lightThreshNauplius", "lightThreshCopepodid",
                         "salinityThreshMin", "salinityThreshMax",
                         "D_hVert")),
                ~if_else(fixDepth, NA, .x)),
         fixDepth=if_else(fixDepth, "2D", "3D"))

p_ls <- vector("list", nrow(var_pretty))
for(i in seq_along(p_ls)) {
  if(var_pretty$varNumeric[i]) {
    if(grepl("Diffusion", var_pretty$var_pretty[i])) {
      x_scale <- scale_x_log10()
    } else {
      x_scale <- scale_x_continuous()
    }
    p_ls[[i]] <- ggplot(par_dist_df, aes(.data[[var_pretty$var[i]]])) +
      geom_histogram(aes(y=after_stat(count)/sum(after_stat(count))), 
                     colour="grey30", fill=modType3_cols[["Constituent"]], bins=9) +
      geom_rug(data=sim_used, aes(.data[[var_pretty$var[i]]]), sides="b", 
               colour=modType3_cols[["Constituent"]]) +
      geom_rug(data=sim_used |> filter(Simulation=="3D.7"), 
               sides="b", colour=modType3_cols[["Opt['Param']"]], linewidth=1) +
      x_scale +
      labs(subtitle=var_pretty$var_pretty[i],
           x=var_pretty$units[[i]],
           y="Probability") + 
      ylim(0, 0.21) 
  } else {
    p_ls[[i]] <- ggplot(par_dist_df, aes(.data[[var_pretty$var[i]]])) +
      geom_bar(aes(y=after_stat(count)/sum(after_stat(count))), 
               colour="grey30", fill=modType3_cols[["Constituent"]]) +
      geom_point(data=sim_used |> group_by(.data[[var_pretty$var[i]]]) |> 
                   summarise(N=n()) |> ungroup() |> mutate(p=N/sum(N)),
                 aes(y=p), shape=3) +
      geom_rug(data=sim_used |> filter(Simulation=="3D.7"), 
               sides="b", colour=modType3_cols[["Opt['Param']"]], linewidth=1) +
      labs(subtitle=var_pretty$var_pretty[i],
           x=var_pretty$units[[i]],
           y="Probability") + 
      ylim(0, 1)
  }
}
p <- cowplot::plot_grid(plotlist=p_ls, nrow=4, align="hv", axis="tblr")
ggsave("figs/pub_new/parameter_distributions.png", p, width=10, height=14)

p_cor_sal <- par_dist_df |>
  select(starts_with("salinity")) |>
  ggplot(aes(salinityThreshMin, salinityThreshMax)) + 
  geom_point(shape=1, alpha=0.05, size=0.5) +
  xlim(20, 28) + ylim(20, 32) +
  labs(x="Lower salinity threshold (100% sink)",
       y="Upper salinity threshold (0% sink)")
p_cor_swim <- par_dist_df |>
  select(starts_with("swim")) |>
  ggplot(aes(swimUpSpeedMean, swimDownSpeedMean)) + 
  geom_point(shape=1, alpha=0.05, size=0.5) +
  labs(x="Upward swim speed (copepodid)",
       y="Downward swim speed (copepodid)")
p_cor_light <- par_dist_df |>
  select(starts_with("light")) |>
  ggplot(aes(lightThreshNauplius, lightThreshCopepodid)) + 
  geom_point(shape=1, alpha=0.05, size=0.5) +
  labs(x="Light threshold (nauplius)",
       y="Light threshold (copepodid)")

p <- cowplot::plot_grid(p_cor_sal, p_cor_swim, p_cor_light, nrow=3, align="hv", axis="tblr")
ggsave("figs/pub_new/parameter_correlations.png", p, width=3.75, height=10)



# ensemble results --------------------------------------------------------

ensCV_df <- ensFull_df |> 
  select(rowNum, sepaSite, CV_k, year, date, licePerFish_rtrt, lice_g05) |>
  left_join(read_csv("out/candidates/CV_candidate_predictions.csv")) |>
  left_join(read_csv("out/ensembles/CV_avg_predictions.csv") |>
              select(rowNum, IP_sim_avgAll)) |>
  left_join(read_csv("out/ensembles/CV_ensBlend_predictions.csv") |>
              select(rowNum, IP_D4_n20) |> rename(IP_predBlend=IP_D4_n20)) |>
  left_join(read_csv("out/ensembles/CV_ensFc-1_rmse.csv") |>
              select(rowNum, .pred) |> rename(IP_predFwk_RMSE=.pred)) |>
  left_join(read_csv("out/ensembles/CV_ensFc-1_rsq.csv") |>
              select(rowNum, .pred) |> rename(IP_predFwk_rsq=.pred)) |>
  left_join(read_csv("out/ensembles/CV_ensFc-1_roc_auc.csv") |>
              select(rowNum, .pred_TRUE) |> rename(IP_predFwk_ROC=.pred_TRUE)) |>
  left_join(read_csv("out/ensembles/CV_ensFc-1_average_precision.csv") |>
              select(rowNum, .pred_TRUE) |> rename(IP_predFwk_PR=.pred_TRUE)) |>
  left_join(read_csv("out/ensembles/CV_ensFc-CV-farms_rmse.csv") |>
              select(rowNum, .pred) |> rename(IP_predFf_RMSE=.pred)) |>
  left_join(read_csv("out/ensembles/CV_ensFc-CV-farms_rsq.csv") |>
              select(rowNum, .pred) |> rename(IP_predFf_rsq=.pred)) |>
  left_join(read_csv("out/ensembles/CV_ensFc-CV-farms_roc_auc.csv") |>
              select(rowNum, .pred_TRUE) |> rename(IP_predFf_ROC=.pred_TRUE)) |>
  left_join(read_csv("out/ensembles/CV_ensFc-CV-farms_average_precision.csv") |>
              select(rowNum, .pred_TRUE) |> rename(IP_predFf_PR=.pred_TRUE))

folds <- unique(ensFull_df$CV_k)
ensNull_0 <- ensNull_time <- ensNull_farm <- vector("list", length(folds))
for(k in seq_along(folds)) {
  ensNull_0[[k]] <- ensFull_df |>
    filter(CV_k != folds[k]) |>
    summarise(IP_null0=mean(licePerFish_rtrt)) |>
    ungroup() |>
    mutate(CV_k=folds[k])
  ensNull_time[[k]] <- ensFull_df |>
    filter(CV_k != folds[k]) |>
    mutate(week=floor(week(date)/2)) |>
    group_by(week) |>
    summarise(IP_nullTime=mean(licePerFish_rtrt)) |>
    ungroup() |>
    mutate(CV_k=folds[k])
  ensNull_farm[[k]] <- ensFull_df |>
    filter(CV_k != folds[k]) |>
    group_by(sepaSite) |>
    summarise(IP_nullFarm=mean(licePerFish_rtrt)) |>
    ungroup() |>
    full_join(site_i |> select(sepaSite), by=join_by(sepaSite)) |>
    mutate(IP_nullFarm=if_else(is.na(IP_nullFarm), mean(IP_nullFarm, na.rm=T), IP_nullFarm),
           CV_k=folds[k])
}
ensCV_df <- ensCV_df |>
  left_join(reduce(ensNull_0, bind_rows), by=join_by(CV_k)) |>
  mutate(week=floor(week(date)/2)) |>
  left_join(reduce(ensNull_time, bind_rows), by=join_by(CV_k, week)) |>
  select(-week) |>
  left_join(reduce(ensNull_farm, bind_rows), by=join_by(CV_k, sepaSite))


write_csv(ensCV_df, "out/ensemble_CV.csv")



# performance plot --------------------------------------------------------

ensCV_df <- read_csv("out/ensemble_CV.csv") |>
  filter(date >= "2021-05-01") |>
  mutate(lice_g05=factor(lice_g05)) 

# Mean within site
metrics_by_farm <- ensCV_df |>
  pivot_longer(starts_with("IP_"), names_to="sim") |>
  mutate(sim=str_remove(sim, "IP_")) |>
  group_by(sepaSite, sim) |>
  summarise(rmse=rmse_vec(value, truth=licePerFish_rtrt),
            rho=cor(value, licePerFish_rtrt, method="spearman", use="pairwise"),
            ROC_AUC=roc_auc_vec(value, truth=lice_g05, event_level="second"),
            N=n(),
            prop_g05=mean(lice_g05=="TRUE"),
            prop_0=mean(licePerFish_rtrt==0)) |>
  ungroup()
ensF_metrics_by_farm <- metrics_by_farm |>
  filter(grepl("predF", sim)) |>
  group_by(sepaSite) |>
  summarise(rmse=sum(rmse * (sim=="predFf_RMSE")),
            rho=sum(rho * (sim=="predFf_rsq")),
            ROC_AUC=sum(ROC_AUC * (sim=="predFf_ROC")),
            sim="predFcst", 
            across(any_of(c("N", "prop_g05", "prop_0", "minPRAUC")), first))
metrics_by_farm <- metrics_by_farm |>
  filter(!grepl("predF", sim)) |>
  bind_rows(ensF_metrics_by_farm) |>
  mutate(rho=if_else(is.na(rho), 0, rho),
         ROC_AUC=if_else(is.na(ROC_AUC), 0.5, ROC_AUC))

# Mean among site
metrics_by_week <- ensCV_df |>
    pivot_longer(starts_with("IP_"), names_to="sim") |>
    mutate(sim=str_remove(sim, "IP_")) |>
    group_by(date, sim) |>
    summarise(rmse=rmse_vec(value, truth=licePerFish_rtrt),
              rho=cor(value, licePerFish_rtrt, method="spearman", use="pairwise"),
              ROC_AUC=roc_auc_vec(value, truth=lice_g05, event_level="second"),
              N=n(),
              prop_g05=mean(lice_g05=="TRUE"),
              prop_0=mean(licePerFish_rtrt==0)) |>
    ungroup()
ensF_metrics_by_week <- metrics_by_week |>
  filter(grepl("predF", sim)) |>
  group_by(date) |>
  summarise(rmse=sum(rmse * (sim=="predFwk_RMSE")),
            rho=sum(rho * (sim=="predFwk_rsq")),
            ROC_AUC=sum(ROC_AUC * (sim=="predFwk_ROC")),
            sim="predFcst", 
            across(any_of(c("N", "prop_g05", "prop_0", "minPRAUC")), first))
metrics_by_week <- metrics_by_week |>
  filter(!grepl("predF", sim)) |>
  bind_rows(ensF_metrics_by_week) |>
  mutate(rho=if_else(is.na(rho), 0, rho),
         ROC_AUC=if_else(is.na(ROC_AUC), 0.5, ROC_AUC))

# Medians
metrics_by_farm_md <- metrics_by_farm |>
  filter(N >= 10) |>
  group_by(sim) |>
  summarise(rmse=median(rmse, na.rm=T),
            rho=median(rho, na.rm=T),
            ROC_AUC=median(ROC_AUC, na.rm=T),
            N=mean(N, na.rm=T),
            prop_g05=mean(prop_g05),
            prop_0=mean(prop_0, na.rm=T)) |>
  ungroup()
metrics_by_week_md <- metrics_by_week |>
  filter(N >= 10) |>
  group_by(sim) |>
  summarise(rmse=median(rmse, na.rm=T),
            rho=median(rho, na.rm=T),
            ROC_AUC=median(ROC_AUC, na.rm=T),
            N=mean(N, na.rm=T),
            prop_g05=mean(prop_g05),
            prop_0=mean(prop_0, na.rm=T)) |>
  ungroup()
metrics_by_farm_mn <- metrics_by_farm |>
  filter(N >= 10) |>
  group_by(sim) |>
  summarise(rmse=mean(rmse, na.rm=T),
            rho=mean(rho, na.rm=T),
            ROC_AUC=mean(ROC_AUC, na.rm=T),
            N=mean(N, na.rm=T),
            prop_g05=mean(prop_g05),
            prop_0=mean(prop_0, na.rm=T)) |>
  ungroup()
metrics_by_week_mn <- metrics_by_week |>
  filter(N >= 10) |>
  group_by(sim) |>
  summarise(rmse=mean(rmse, na.rm=T),
            rho=mean(rho, na.rm=T),
            ROC_AUC=mean(ROC_AUC, na.rm=T),
            N=mean(N, na.rm=T),
            prop_g05=mean(prop_g05),
            prop_0=mean(prop_0, na.rm=T)) |>
  ungroup()

highlight_sims <- c("sim_03", "sim_07", #paste0("sLonLatD3_n", c(5, 10, 20)), 
                    "predBlend", "predFcst",
                    "D3_n20", "D4_n20", "D8_n20",
                    paste0("predFwk_", c("RMSE", "rsq", "ROC")),
                    paste0("predFf_", c("RMSE", "rsq", "ROC")))
map(list(metrics_by_farm_md, metrics_by_week_md),
    ~plot_metric_ordered(.x, rmse, highlight_sims)) |>
  plot_grid(plotlist=_, ncol=2)
map(list(metrics_by_farm_mn, metrics_by_week_mn),
    ~plot_metric_ordered(.x, rmse, highlight_sims)) |>
  plot_grid(plotlist=_, ncol=2)

map(list(metrics_by_farm_md, metrics_by_week_md),
    ~plot_metric_ordered(.x, rho, highlight_sims)) |>
  plot_grid(plotlist=_, ncol=2)
map(list(metrics_by_farm_mn, metrics_by_week_mn),
    ~plot_metric_ordered(.x, rho, highlight_sims)) |>
  plot_grid(plotlist=_, ncol=2)

map(list(metrics_by_farm_md, metrics_by_week_md),
    ~plot_metric_ordered(.x, ROC_AUC, highlight_sims)) |>
  plot_grid(plotlist=_, ncol=2)
map(list(metrics_by_farm_mn, metrics_by_week_mn),
    ~plot_metric_ordered(.x, ROC_AUC, highlight_sims)) |>
  plot_grid(plotlist=_, ncol=2)


metric_ranks <- bind_rows(
  metrics_by_farm |>
    filter(sim %in% c("predFcst", "predBlend", "sim_avgAll",
                      paste0("sim_0", 1:9), paste0("sim_", 10:20),
                      paste0("D", 3:8, "_n20"))) |>
    select(sepaSite, sim, N, rmse, rho, ROC_AUC) |>
    pivot_longer(4:6, names_to="metric", values_to="value") |>
    mutate(value_lowGood=if_else(metric=="rmse", value, -value),
           type="byFarm") |>
    drop_na() |>
    group_by(sepaSite, metric) |>
    mutate(rank=min_rank(value_lowGood)) |>
    ungroup(),
  metrics_by_week |>
    filter(sim %in% c("predFcst", "predBlend", "sim_avgAll",
                      paste0("sim_0", 1:9), paste0("sim_", 10:20),
                      paste0("D", 3:8, "_n20"))) |>
    select(date, sim, N, rmse, rho, ROC_AUC) |>
    pivot_longer(4:6, names_to="metric", values_to="value") |>
    mutate(value_lowGood=if_else(metric=="rmse", value, -value),
           type="byWeek") |>
    drop_na() |>
    group_by(date, metric) |>
    mutate(rank=min_rank(value_lowGood)) |>
    ungroup()
  ) |>
  left_join(sim_i)


metric_ranks |>
  filter(N >= 10) |>
  group_by(metric, sim, lab, lab_short, type) |>
  summarise(mn=mean(rank, na.rm=T)) |>
  ungroup() |>
  filter(grepl("D|Avg", lab)) |>
  group_by(type, metric) |>
  arrange(lab) |>
  summarise(pBetterThan=sum(first(mn) < mn)/20)

metric_ranks |>
  filter(N >= 10) |>
  group_by(metric, sim, lab, lab_short, type) |>
  summarise(mn=mean(rank, na.rm=T)) |>
  ungroup() |>
  filter(grepl("D|Blend", lab)) |>
  group_by(type, metric) |>
  arrange(lab) |>
  summarise(pBetterThan=sum(first(mn) < mn)/20)

metric_ranks |>
  filter(N >= 10) |>
  filter(grepl("D", lab)) |>
  group_by(metric, type, date, sepaSite) |>
  mutate(rank=min_rank(rank)) |>
  ungroup() |>
  filter(rank==1) |>
  count(sim) |>
  mutate(p=n/sum(n)*100)

metric_ranks |>
  filter(N >= 10) |>
  group_by(sim, lab, lab_short, type) |>
  summarise(mn=mean(rank, na.rm=T)) |>
  ungroup() |>
  mutate(type=factor(type,
                     levels=c("global", "byFarm", "byWeek"),
                     labels=c("Global", "By farm (median)", "By week (median)"))) |>
  ggplot(aes(mn, sim)) +
  geom_point() +
  geom_rug(sides="b") +
  facet_grid(type~.) +
  theme_bw() +
  theme(panel.grid.major=element_line(linewidth=0.5, colour="grey80"),
        panel.grid.minor=element_line(linewidth=0.1, colour="grey95"))

metric_ranks |>
  filter(N >= 10) |>
  group_by(sim, lab, lab_short, type) |>
  summarise(mn=mean(rank, na.rm=T)) |>
  group_by(sim) |>
  summarise(mn=mean(mn)) |>
  arrange(mn)

metric_ranks |>
  ggplot(aes(rank, sim)) +
  geom_boxplot()
metric_ranks |>
  ggplot(aes(rank, sim)) + 
  stat_halfeye()

all_metrics_df <- bind_rows(
  metrics_by_farm_mn |> mutate(type="byFarm"),
  metrics_by_week_mn |> mutate(type="byWeek")
) |>
  filter(sim != "null0") |>
  pivot_longer(any_of(c("rmse", "rho", "ROC_AUC")), names_to="metric") |>
  mutate(metric=factor(metric, levels=c("ROC_AUC", "rho", "rmse"),
                       labels=c("'AUC'['ROC']", "rho", "RMSE"))) |>
  left_join(sim_i) |>
  arrange(lab) |>
  mutate(type=factor(type, 
                     levels=c("byFarm", "byWeek"),
                     labels=c("By farm", "By week"))) |>
  drop_na() |>
  arrange(desc(lab)) 

all_metrics_labs <- all_metrics_df |>
  filter(metric=="RMSE",
         type=="By farm",
         grepl("Ens", lab_short)) |>
  arrange(lab) |>
  mutate(label=c("Ens['Fcst']", "Ens['Blend']", "Ens['Avg']")) |>
  bind_rows(tibble(sim=c("3D.1", "2D.1"),
                   N=1, prop_g05=1, prop_0=1,
                   type="By farm",
                   metric="RMSE",
                   lab=c("3D.1", "2D.1"),
                   lab_short=c("3D", "2D"),
                   label=c("'3D'", "'2D'"))) |>
  mutate(value=seq(0.975, 0.84, length.out=n()))
  # bind_rows(., 
  #           . |> 
  #             filter(grepl("Mean", lab_short)) |>
  #             mutate(lab=c("3D.1", "2D.1"),
  #                    lab_short=c("3D", "2D"),
  #                    label=c(NA, NA)))


ms_rmse <- all_metrics_df |> filter(metric=="RMSE") |>
  filter(!grepl("null", sim)) |>
  metric_plot_base(theme="ms", modType3_cols) + 
  scale_y_continuous("RMSE", limits=c(0.25, 0.35), oob=scales::oob_keep, expand=c(0,0),
                     breaks=seq(0, 1, by=0.05), minor_breaks=seq(0, 1, by=0.01)) +
  scale_x_discrete(limits=c("By farm", "By week"), 
                   labels=c("Mean\nwithin\nfarm", "Mean\nwithin\nweek"))
ms_r <- all_metrics_df |> filter(metric=="rho") |>
  filter(!grepl("null", sim)) |>
  metric_plot_base(theme="ms", modType3_cols) + 
  scale_y_continuous(expression('Spearmans'~~rho), limits=c(0, 1), oob=scales::oob_keep, expand=c(0,0),
                     breaks=seq(0, 1, by=0.25), minor_breaks=seq(0, 1, by=0.05)) +
  scale_x_discrete(limits=c("By farm", "By week"), 
                   labels=c("Mean\nwithin\nfarm", "Mean\nwithin\nweek"))
ms_ROC <- all_metrics_df |> filter(metric=="'AUC'['ROC']") |>
  filter(!grepl("null", sim)) |>
  metric_plot_base(theme="ms", modType3_cols) + 
  scale_y_continuous(expression('AUC'['ROC']), limits=c(0.5, 1), oob=scales::oob_keep, expand=c(0,0),
                     breaks=seq(0.5, 1, by=0.1), minor_breaks=seq(0.5, 1, by=0.02)) +
  scale_x_discrete(limits=c("By farm", "By week"), 
                   labels=c("Mean\nwithin\nfarm", "Mean\nwithin\nweek"))

ms_legend <- all_metrics_labs |>
  filter(!grepl("null", sim)) |>
  mutate(label=factor(label, levels=unique(label)),
         lab_short=factor(lab_short, levels=unique(lab_short))) |>
  ggplot() +
  geom_text(aes(type, value, label=label, colour=lab_short),
            hjust=0, nudge_x=-0.15, vjust=0.5, size=2.5, parse=T) +
  geom_point(position=position_nudge(x=-0.35), stroke=0.7,
             aes(type, value, colour=lab_short, shape=lab_short, size=lab_short)) +
  scale_colour_manual(values=modType3_cols) +
  # scale_colour_manual(values=c("black", "#b2182b", "#d6604d",
  #                              scico(2, begin=0.2, end=0.7, palette="broc", direction=1))) +
  scale_shape_manual(values=c(1, 1, 1, 4, 3)) +
  scale_size_manual(values=c(rep(2.5, 3), rep(1, 2))) +
  scale_alpha_manual(values=c(1, 1, 1, 0.5, 0.5)) +
  ylim(0.575, 1.175) +
  theme(legend.position="none",
        plot.margin=margin(t=0, b=0, l=0, r=0),
        panel.border=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())
p <- plot_grid(ms_rmse, ms_r, ms_ROC, ms_legend, 
               align="h", axis="tb", nrow=1, rel_widths=c(1.12, 1.12, 1.12, 0.4))
ggsave("figs/pub_new/validation_metrics_CV_means.png", p, width=6, height=4)






# weekly performance ------------------------------------------------------

metric_date_df <- ensCV_df |>
  filter(date >= "2021-05-01") |>
  group_by(date) |>
  summarise(N=n(),
            N_0=sum(licePerFish_rtrt == 0),
            N_True=sum(lice_g05=="TRUE"),
            N_False=sum(lice_g05=="FALSE"),
            mn_lpf=mean(licePerFish_rtrt),
            rmse=rmse_vec(IP_predBlend, truth=licePerFish_rtrt),
            rho=cor(IP_predBlend, licePerFish_rtrt, method="spearman", use="pairwise"),
            ROC_AUC=roc_auc_vec(IP_predBlend, truth=lice_g05, event_level="second")) |>
  ungroup() |>
  mutate(sim="predBlend") |>
  bind_rows(
    ensCV_df |>
      filter(date >= "2021-05-01") |>
      group_by(date) |>
      summarise(N=n(),
                N_0=sum(licePerFish_rtrt == 0),
                N_True=sum(lice_g05=="TRUE"),
                N_False=sum(lice_g05=="FALSE"),
                mn_lpf=mean(licePerFish_rtrt),
                rmse=rmse_vec(IP_predFwk_RMSE, truth=licePerFish_rtrt),
                rho=cor(IP_predFwk_rsq, licePerFish_rtrt, method="spearman", use="pairwise"),
                ROC_AUC=roc_auc_vec(IP_predFwk_ROC, truth=lice_g05, event_level="second")) |>
      ungroup() |>
      mutate(sim="predFcst")
  ) |>
  bind_rows(
    ensCV_df |>
      filter(date >= "2021-05-01") |>
      group_by(date) |>
      summarise(N=n(),
                N_0=sum(licePerFish_rtrt == 0),
                N_True=sum(lice_g05=="TRUE"),
                N_False=sum(lice_g05=="FALSE"),
                mn_lpf=mean(licePerFish_rtrt),
                rmse=rmse_vec(IP_sim_07, truth=licePerFish_rtrt),
                rho=cor(IP_sim_07, licePerFish_rtrt, method="spearman", use="pairwise"),
                ROC_AUC=roc_auc_vec(IP_sim_07, truth=lice_g05, event_level="second")) |>
      ungroup() |>
      mutate(sim="sim_07")
  ) |>
  bind_rows(
    ensCV_df |>
      filter(date >= "2021-05-01") |>
      group_by(date) |>
      summarise(N=n(),
                N_0=sum(licePerFish_rtrt == 0),
                N_True=sum(lice_g05=="TRUE"),
                N_False=sum(lice_g05=="FALSE"),
                mn_lpf=mean(licePerFish_rtrt),
                rmse=rmse_vec(IP_sim_avgAll, truth=licePerFish_rtrt),
                rho=cor(IP_sim_avgAll, licePerFish_rtrt, method="spearman", use="pairwise"),
                ROC_AUC=roc_auc_vec(IP_sim_avgAll, truth=lice_g05, event_level="second")) |>
      ungroup() |>
      mutate(sim="sim_avgAll")
  ) |>
  filter(N >= 10) |>
  pivot_longer(7:9, names_to="metric") |>
  filter(metric != "ROC_AUC" | (N_True > 0 & N_False > 0)) |>
  mutate(metric=factor(metric, levels=c("ROC_AUC", "rho", "rmse"),
                       labels=c("'AUC'['ROC']", "rho", "RMSE"))) |>
  left_join(sim_i) |>
  droplevels()

# Performance summaries: values
all_metrics_df |> 
  filter(sim %in% c("predBlend", "predFcst", "sim_03", "sim_07")) |> 
  arrange(type, metric, value) |>
  print(n=50)


metric_date_df |> 
  ggplot(aes(mn_lpf, value)) +
  geom_point(shape=1) + 
  stat_smooth(method="lm") + 
  facet_grid(metric~sim, scales="free_y") 

metric_date_df |> 
  ggplot(aes(N, value)) +
  geom_point(shape=1) + 
  stat_smooth(method="lm") + 
  facet_grid(metric~sim, scales="free_y") 

metric_date_df |> 
  ggplot(aes(N_True/N, value)) +
  geom_point(shape=1) + 
  stat_smooth(method="lm") + 
  facet_grid(metric~sim, scales="free_y") 

metric_date_df |>
  group_by(metric, sim) |>
  summarise(lpf_r=cor(value, mn_lpf),
            N_r=cor(value, N),
            pG05_r=cor(value, N_True/N)) |>
  pivot_longer(ends_with("_r")) |>
  ggplot(aes(metric, sim, fill=abs(value))) +
  geom_raster() + 
  colorspace::scale_fill_continuous_sequential(palette="heat", limits=c(0,1)) +
  facet_wrap(~name)

# farm performance --------------------------------------------------------

metric_farm_df <- ensCV_df |>
  filter(date >= "2021-05-01") |>
  group_by(sepaSite) |>
  summarise(N=n(),
            N_0=sum(licePerFish_rtrt == 0),
            N_True=sum(lice_g05=="TRUE"),
            N_False=sum(lice_g05=="FALSE"),
            mn_lpf=mean(licePerFish_rtrt),
            rmse=rmse_vec(IP_predBlend, truth=licePerFish_rtrt),
            rho=cor(IP_predBlend, licePerFish_rtrt, method="spearman", use="pairwise"),
            ROC_AUC=roc_auc_vec(IP_predBlend, truth=lice_g05, event_level="second")) |>
  ungroup() |>
  mutate(sim="predBlend") |>
  bind_rows(
    ensCV_df |>
      filter(date > "2021-05-01") |>
      group_by(sepaSite) |>
      summarise(N=n(),
                N_0=sum(licePerFish_rtrt == 0),
                N_True=sum(lice_g05=="TRUE"),
                N_False=sum(lice_g05=="FALSE"),
                mn_lpf=mean(licePerFish_rtrt),
                rmse=rmse_vec(IP_predFf_RMSE, truth=licePerFish_rtrt),
                rho=cor(IP_predFf_rsq, licePerFish_rtrt, method="spearman", use="pairwise"),
                ROC_AUC=roc_auc_vec(IP_predFf_ROC, truth=lice_g05, event_level="second")) |>
      ungroup() |>
      mutate(sim="predFcst")
  ) |>
  bind_rows(
    ensCV_df |>
      filter(date >= "2021-05-01") |>
      group_by(sepaSite) |>
      summarise(N=n(),
                N_0=sum(licePerFish_rtrt == 0),
                N_True=sum(lice_g05=="TRUE"),
                N_False=sum(lice_g05=="FALSE"),
                mn_lpf=mean(licePerFish_rtrt),
                rmse=rmse_vec(IP_sim_07, truth=licePerFish_rtrt),
                rho=cor(IP_sim_07, licePerFish_rtrt, method="spearman", use="pairwise"),
                ROC_AUC=roc_auc_vec(IP_sim_07, truth=lice_g05, event_level="second")) |>
      ungroup() |>
      mutate(sim="sim_07")
  ) |>
  bind_rows(
    ensCV_df |>
      filter(date >= "2021-05-01") |>
      group_by(sepaSite) |>
      summarise(N=n(),
                N_0=sum(licePerFish_rtrt == 0),
                N_True=sum(lice_g05=="TRUE"),
                N_False=sum(lice_g05=="FALSE"),
                mn_lpf=mean(licePerFish_rtrt),
                rmse=rmse_vec(IP_sim_avgAll, truth=licePerFish_rtrt),
                rho=cor(IP_sim_avgAll, licePerFish_rtrt, method="spearman", use="pairwise"),
                ROC_AUC=roc_auc_vec(IP_sim_avgAll, truth=lice_g05, event_level="second")) |>
      ungroup() |>
      mutate(sim="sim_avgAll")
  ) |>
  filter(N >= 10) |>
  pivot_longer(7:9, names_to="metric") |>
  filter(metric != "ROC_AUC" | (N_True > 0 & N_False > 0)) |>
  mutate(metric=factor(metric, levels=c("ROC_AUC", "rho", "rmse"),
                       labels=c("'AUC'['ROC']", "rho", "RMSE"))) |>
  left_join(sim_i) |>
  droplevels()


bind_rows(
  metric_farm_df |> mutate(type="Within farm"),
  metric_date_df |> mutate(type="Within week")
) |>
  group_by(metric, sim, type) |>
  summarise(lpf_r=cor(value, mn_lpf),
            N_r=cor(value, N),
            pG0_r=cor(value, (N-N_0)/N),
            pG05_r=cor(value, N_True/N)) |>
  pivot_longer(ends_with("_r")) |>
  ggplot(aes(metric, sim, fill=value)) +
  geom_tile(colour="grey30") + 
  colorspace::scale_fill_binned_diverging(palette="Blue-Red 3", limits=c(-1,1), 
                                          breaks=seq(-1,1,by=0.25), rev=T) +
  facet_grid(type~name)

p <- bind_rows(
  metric_farm_df |> mutate(type="Within farm"),
  metric_date_df |> mutate(type="Within week")
) |>
  group_by(metric, lab, type) |>
  summarise(lpf_r=cor(value, mn_lpf),
            N_r=cor(value, N),
            pG0_r=cor(value, (N-N_0)/N),
            pG05_r=cor(value, N_True/N)) |>
  pivot_longer(ends_with("_r")) |>
  mutate(lab=if_else(lab=="'3D.7'", "Opt['Param']", lab),
         lab=factor(lab, levels=names(modType3_cols))) |>
  ggplot(aes(name, value, colour=lab, shape=type)) + 
  geom_jitter(size=2, stroke=1, height=0, width=0.1) + 
  scale_shape_manual(values=c(1, 2)) +
  scale_colour_manual(values=modType3_cols, 
                      labels=c(expression(Ens['Fcst']),
                               expression(Ens['Blend']),
                               expression(Ens['Avg']),
                               expression(Opt['Param']))) +
  scale_x_discrete(breaks=c("lpf_r", "N_r", "pG0_r", "pG05_r"),
                   labels=c("Mean lice per fish", "Number of records", 
                            "Proportion of records > 0 lpf", 
                            "Proportion of records > 0.5 lpf") |>
                     str_wrap(width=10)) +
  scale_y_continuous("Pearson's correlation coefficient", limits=c(-1, 1)) +
  facet_grid(.~metric, labeller=label_parsed) +
  theme(panel.grid.major.y=element_line(colour="grey90", linewidth=0.3),
        axis.title.x=element_blank(),
        legend.title=element_blank())
ggsave("figs/pub_new/metric_correlations.png", p, width=12, height=4)


# resampling evaluation ---------------------------------------------------

resampleEval_df <- read_csv("out/ensBlend_resample_performance.csv")


point_df <- resampleEval_df |> 
  filter(nSim != "n20") |>
  # mean rank among weeks or farms
  group_by(metric, type, nSim, modType, name, sample) |> 
  summarise(mnRank=mean(rank), mdRank=median(rank), value=mean(value), skill=mean(skill)) |>
  ungroup() |>
  mutate(modType=factor(modType, 
                        levels=c("ens", "avg", "cand"),
                        labels=c("Ens['Blend']", "Ens['Avg']", "Constituent")),
         type=factor(type, levels=c("date", "site"),
                     labels=paste0("Mean within-", c("week", "farm"))),
         metric=factor(metric, levels=c("RMSE", "r", "ROC_AUC"),
                       labels=c("RMSE", "rho", "'AUC'['ROC']")),
         sample=factor(sample, levels=1:100),
         nSim=factor(nSim, levels=paste0("n", c(3, 5, 10, 15, 20)),
                     labels=paste0("n: ", c(3, 5, 10, 15, 20)))) |>
  arrange(metric, type, nSim, sample, skill, desc(name)) |>
  group_by(metric, type, nSim, sample) |> 
  mutate(bestMod=last(name),
         bestMod=case_when(bestMod=="ens_pred" ~ "Ens['Blend']",
                           bestMod=="avg_pred" ~ "Ens['Avg']",
                           .default="Constituent"),
         bestMod=factor(bestMod, levels=c("Ens['Blend']", "Ens['Avg']", "Constituent"))) |>
  ungroup() |>
  group_by(metric, type, nSim, sample, modType) |>
  mutate(modType2=case_when(modType=="Ens['Blend']" ~ "Ens['Blend']",
                            modType=="Ens['Avg']" ~ "Ens['Avg']",
                            modType=="Constituent" & skill==max(skill) ~ "Opt['Param']",
                            modType=="Constituent" & skill < max(skill) ~ "Other")) |>
  ungroup() |>
  mutate(modType2=factor(modType2, 
                         levels=c("Ens['Blend']", "Ens['Avg']", "Opt['Param']", "Other"),
                         labels=c("Ens['Blend']", "Ens['Avg']", "Opt['Param']", "Other")))

const_z <- point_df |>
  filter(modType=="Constituent") |>
  group_by(nSim, sample, metric, type) |>
  summarise(mn=mean(skill),
            sd=sd(skill)) |>
  ungroup()



p <- resampleEval_df |>
  filter(nSim != "n20") |>
  group_by(metric, type, nSim, modType, name, sample) |> 
  summarise(mnSkill=mean(skill)) |> 
  group_by(metric, type, nSim, sample) |> 
  mutate(rank=min_rank(desc(mnSkill))) |>
  ungroup() |>
  mutate(modType=factor(modType, 
                        levels=c("ens", "avg", "cand"),
                        labels=c("Ens['Blend']", "Ens['Avg']", "Constituent")),
         type=factor(type, levels=c("date", "site"),
                     labels=paste0("Mean within-", c("week", "farm"))),
         metric=factor(metric, levels=c("RMSE", "r", "ROC_AUC"),
                       labels=c("RMSE", "rho", "'AUC'['ROC']")),
         sample=factor(sample, levels=1:100),
         nSim=factor(nSim, levels=paste0("n", c(15, 10, 5, 3)),
                     labels=paste0(c(15, 10, 5, 3)))) |>
  ggplot(aes(rank, nSim, fill=modType)) + 
  ggdist::geom_dots(layout="bar", aes(group=nSim), side="both", slab_colour="grey30", slab_linewidth=0.05) + 
  scale_fill_manual("Model type", values=modType2_cols,
                    labels=lab_expressions) +
  scale_x_continuous("Rank within resample", breaks=seq(1,17,by=2)) +
  scale_y_discrete("Number of constituents") +
  facet_grid(metric~type, labeller=labeller(metric=label_parsed)) +
  theme(panel.grid.major.x=element_line(colour="grey80", linewidth=0.25),
        legend.position="inside",
        legend.position.inside=c(0.875, 0.95),
        legend.background=element_blank(),
        legend.title=element_blank())
ggsave("figs/pub_new/ensBlend_eval_dotbars.png", p, width=6, height=9, dpi=400)




pA <- point_df |>
  filter(type=="Mean within-week") |>
  left_join(const_z, by=join_by(nSim, sample, metric, type)) |>
  mutate(skill_z=(skill - mn)/sd,
         nSim=as.numeric(nSim)) |>
  filter(abs(skill_z) <= 6) |>
  ggplot(aes(skill_z, nSim, fill=modType2, colour=modType2)) +
  geom_vline(xintercept=0, linewidth=0.5, colour="grey85") +
  geom_dots(aes(group=nSim), side="top", layout="hex", alpha=0.9, binwidth = unit(c(0.3, Inf), "mm"),
            overflow="compress", colour=NA, orientation="horizontal", scale=0.8) +
  stat_pointinterval(aes(y=nSim-0.15, group=paste(nSim, modType2)), 
                     .width=c(0.8), linewidth=0.5,
                     position=position_dodge(width=0.25), shape=1, size=0.75) + 
  scale_fill_manual("Model type", values=modType2_cols,
                    labels=lab2_expressions) +
  scale_colour_manual("Model type", values=modType2_cols,
                      labels=lab2_expressions) +
  scale_x_continuous(" ") +
  scale_y_continuous("Number of constituents", 
                     breaks=1:4, labels=c("3", "5", "10", "15")) +
  facet_grid(type~metric, scales="free_x", labeller=labeller(metric=label_parsed)) + 
  theme_bw() + 
  theme(legend.position="inside",
        legend.position.inside=c(0.94, 0.85),
        legend.title=element_blank(),
        legend.key.height=unit(2.5, "mm"),
        legend.background=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.y=element_blank())
pB <- point_df |>
  filter(type=="Mean within-farm") |>
  left_join(const_z, by=join_by(nSim, sample, metric, type)) |>
  mutate(skill_z=(skill - mn)/sd,
         nSim=as.numeric(nSim)) |>
  filter(abs(skill_z) <= 6) |>
  ggplot(aes(skill_z, nSim, fill=modType2, colour=modType2)) +
  geom_vline(xintercept=0, linewidth=0.5, colour="grey85") +
  geom_dots(aes(group=nSim), side="top", layout="hex", alpha=0.9, binwidth = unit(c(0.3, Inf), "mm"),
            overflow="compress", colour=NA, orientation="horizontal", scale=0.8) +
  stat_pointinterval(aes(y=nSim-0.15, group=paste(nSim, modType2)), 
                     .width=c(0.8), linewidth=0.5,
                     position=position_dodge(width=0.25), shape=1, size=0.75) + 
  scale_fill_manual("Model type", values=modType2_cols,
                    labels=lab2_expressions) +
  scale_colour_manual("Model type", values=modType2_cols,
                      labels=lab2_expressions) +
  scale_x_continuous("Z-score of skill within resample") +
  scale_y_continuous("Number of constituents", 
                     breaks=1:4, labels=c("3", "5", "10", "15")) +
  facet_grid(type~metric, scales="free_x", labeller=labeller(metric=label_parsed)) + 
  theme_bw() + 
  theme(legend.position="none",
        legend.background=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.y=element_blank())
p <- cowplot::plot_grid(pA, pB, nrow=2, align="v", axis="lr")
ggsave("figs/pub_new/ensBlend_eval_dotplot2_opt.png", p, width=9, height=6, dpi=400)




p1 <- resampleEval_df |>
  filter(nSim != "n20") |>
  filter(modType %in% c("ens", "cand")) |>
  # mean rank among weeks or farms
  group_by(metric, type, nSim, modType, name, sample) |> 
  summarise(mnSkill=mean(skill)) |> 
  arrange(nSim, type, sample, metric, modType) |>
  group_by(nSim, type, sample, metric) |>
  summarise(nBetterThan=sum(last(mnSkill) > mnSkill, na.rm=T)) |>
  mutate(prEnsBetter=nBetterThan/as.numeric(str_sub(nSim, 2, -1)),
         type=factor(type, levels=c("date", "site"),
                     labels=paste0("Mean within-", c("week", "farm"))),
         metric=factor(metric, levels=c("RMSE", "r", "ROC_AUC"),
                       labels=c("RMSE", "rho", "'AUC'['ROC']")),
         nSim=factor(nSim, levels=paste0("n", c(3,5,10,15)), labels=c(3, 5, 10,15))) |>
  ggplot(aes(nSim, fill=prEnsBetter, group=prEnsBetter)) + 
  geom_hline(yintercept=0.5) +
  geom_bar(position="fill", colour="grey30") +
  scale_fill_gradient2("Ensemble percentile    \nvs. constituents", 
                       midpoint=0.5, limits=c(0, 1), labels=scales::label_percent(suffix="")) +
  scale_y_continuous("Percentage of resamples", 
                     breaks=c(0, 0.5, 1),
                     labels=scales::label_percent()) +
  xlab("Number of constituents per resample") +
  facet_grid(type~metric, labeller=labeller(metric=label_parsed)) +
  # facet_grid(.~metric, labeller=labeller(metric=label_parsed)) +
  ggtitle(expression(Ens['Blend'])) +
  theme(panel.grid.major.y=element_line(colour="grey80"),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())

p2 <- resampleEval_df |>
  filter(nSim != "n20") |>
  filter(modType %in% c("avg", "cand")) |>
  # mean rank among weeks or farms
  group_by(metric, type, nSim, modType, name, sample) |> 
  summarise(mnSkill=mean(skill)) |> 
  group_by(nSim, type, sample, metric) |>
  summarise(nBetterThan=sum(first(mnSkill) > mnSkill, na.rm=T)) |>
  mutate(prEnsBetter=nBetterThan/as.numeric(str_sub(nSim, 2, -1)),
         type=factor(type, levels=c("date", "site"),
                     labels=paste0("Mean within-", c("week", "farm"))),
         metric=factor(metric, levels=c("RMSE", "r", "ROC_AUC"),
                       labels=c("RMSE", "rho", "'AUC'['ROC']")),
         nSim=factor(nSim, levels=paste0("n", c(3,5,10,15)), labels=c(3, 5, 10,15))) |>
  ggplot(aes(nSim, fill=prEnsBetter, group=prEnsBetter)) + 
  geom_hline(yintercept=0.5) +
  geom_bar(position="fill", colour="grey30") +
  scale_fill_gradient2("Ensemble percentile    \nvs. constituents", 
                       midpoint=0.5, limits=c(0, 1), labels=scales::label_percent(suffix="")) +
  scale_y_continuous("Percentage of resamples", 
                     breaks=c(0, 0.5, 1),
                     labels=scales::label_percent()) +
  xlab("Number of constituents per resample") +
  facet_grid(type~metric, labeller=labeller(metric=label_parsed)) +
  # facet_grid(.~metric, labeller=labeller(metric=label_parsed)) +
  ggtitle(expression(Ens['Avg'])) +
  theme(panel.grid.major.y=element_line(colour="grey80"),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())

p <- ggpubr::ggarrange(p1, p2, nrow=1, common.legend=T)
ggsave("figs/pub_new/ens_nBetter_meanSkill_type.png", p, width=9, height=8)


pA <- point_df |>
  filter(type=="Mean within-week") |>
  left_join(const_z, by=join_by(nSim, sample, metric, type)) |>
  mutate(skill_z=(skill - mn)/sd,
         nSim=as.numeric(nSim)) |>
  filter(abs(skill_z) <= 6) |>
  ggplot(aes(skill_z, nSim+0.15, fill=modType, colour=modType)) +
  geom_vline(xintercept=0, linewidth=0.5, colour="grey85") +
  geom_dots(aes(order=modType, group=nSim), side="both", layout="hex", alpha=0.9, binwidth = unit(c(0.3, Inf), "mm"),
            overflow="compress", colour=NA, orientation="horizontal", scale=0.8) +
  stat_pointinterval(aes(y=nSim-0.15, group=paste(nSim, modType)), 
                     .width=c(0.8), linewidth=0.5,
                     position=position_dodge(width=0.25), shape=1, size=0.75) + 
  scale_fill_manual("Model type", values=modType_cols,
                    labels=lab_expressions) +
  scale_colour_manual("Model type", values=modType_cols,
                    labels=lab_expressions) +
  scale_x_continuous(" ") +
  scale_y_continuous("Number of constituents", 
                     breaks=1:4, labels=c("3", "5", "10", "15")) +
  facet_grid(type~metric, scales="free_x", labeller=labeller(metric=label_parsed)) + 
  theme_bw() + 
  theme(legend.position="inside",
        legend.position.inside=c(0.93, 0.83),
        legend.title=element_blank(),
        legend.key.height=unit(2.5, "mm"),
        legend.background=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.y=element_blank())
pB <- point_df |>
  filter(type=="Mean within-farm") |>
  left_join(const_z, by=join_by(nSim, sample, metric, type)) |>
  mutate(skill_z=(skill - mn)/sd,
         nSim=as.numeric(nSim)) |>
  filter(abs(skill_z) <= 6) |>
  ggplot(aes(skill_z, nSim+0.15, fill=modType, colour=modType)) +
  geom_vline(xintercept=0, linewidth=0.5, colour="grey85") +
  geom_dots(aes(order=modType, group=nSim), side="both", layout="hex", alpha=0.9, binwidth = unit(c(0.3, Inf), "mm"),
            overflow="compress", colour=NA, orientation="horizontal", scale=0.8) +
  stat_pointinterval(aes(y=nSim-0.15, group=paste(nSim, modType)), 
                     .width=c(0.8), linewidth=0.5,
                     position=position_dodge(width=0.25), shape=1, size=0.75) + 
  scale_fill_manual("Model type", values=modType_cols,
                    labels=lab_expressions) +
  scale_colour_manual("Model type", values=modType_cols,
                      labels=lab_expressions) +
  scale_x_continuous("Z-score of skill within resample") +
  scale_y_continuous("Number of constituents", 
                     breaks=1:4, labels=c("3", "5", "10", "15")) +
  facet_grid(type~metric, scales="free_x", labeller=labeller(metric=label_parsed)) + 
  theme_bw() + 
  theme(legend.position="none",
        legend.background=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.y=element_blank())
p <- cowplot::plot_grid(pA, pB, nrow=2, align="v", axis="lr")
ggsave("figs/pub_new/ensBlend_eval_dotplot.png", p, width=9, height=6, dpi=400)


pA <- point_df |>
  filter(type=="Mean within-week") |>
  left_join(const_z, by=join_by(nSim, sample, metric, type)) |>
  mutate(skill_z=(skill - mn)/sd,
         nSim=as.numeric(nSim)) |>
  filter(abs(skill_z) <= 6) |>
  ggplot(aes(skill_z, nSim, fill=modType, colour=modType)) +
  geom_vline(xintercept=0, linewidth=0.5, colour="grey85") +
  geom_dots(aes(group=nSim), side="top", layout="hex", alpha=0.9, binwidth = unit(c(0.3, Inf), "mm"),
            overflow="compress", colour=NA, orientation="horizontal", scale=0.8) +
  stat_pointinterval(aes(y=nSim-0.15, group=paste(nSim, modType)), 
                     .width=c(0.8), linewidth=0.5,
                     position=position_dodge(width=0.25), shape=1, size=0.75) + 
  scale_fill_manual("Model type", values=modType_cols,
                    labels=lab_expressions) +
  scale_colour_manual("Model type", values=modType_cols,
                      labels=lab_expressions) +
  scale_x_continuous(" ") +
  scale_y_continuous("Number of constituents", 
                     breaks=1:4, labels=c("3", "5", "10", "15")) +
  facet_grid(type~metric, scales="free_x", labeller=labeller(metric=label_parsed)) + 
  theme_bw() + 
  theme(legend.position="inside",
        legend.position.inside=c(0.93, 0.88),
        legend.title=element_blank(),
        legend.key.height=unit(2.5, "mm"),
        legend.background=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.y=element_blank())
pB <- point_df |>
  filter(type=="Mean within-farm") |>
  left_join(const_z, by=join_by(nSim, sample, metric, type)) |>
  mutate(skill_z=(skill - mn)/sd,
         nSim=as.numeric(nSim)) |>
  filter(abs(skill_z) <= 6) |>
  ggplot(aes(skill_z, nSim, fill=modType, colour=modType)) +
  geom_vline(xintercept=0, linewidth=0.5, colour="grey85") +
  geom_dots(aes(group=nSim), side="top", layout="hex", alpha=0.9, binwidth = unit(c(0.3, Inf), "mm"),
            overflow="compress", colour=NA, orientation="horizontal", scale=0.8) +
  stat_pointinterval(aes(y=nSim-0.15, group=paste(nSim, modType)), 
                     .width=c(0.8), linewidth=0.5,
                     position=position_dodge(width=0.25), shape=1, size=0.75) + 
  scale_fill_manual("Model type", values=modType_cols,
                    labels=lab_expressions) +
  scale_colour_manual("Model type", values=modType_cols,
                      labels=lab_expressions) +
  scale_x_continuous("Z-score of skill within resample") +
  scale_y_continuous("Number of constituents", 
                     breaks=1:4, labels=c("3", "5", "10", "15")) +
  facet_grid(type~metric, scales="free_x", labeller=labeller(metric=label_parsed)) + 
  theme_bw() + 
  theme(legend.position="none",
        legend.background=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.y=element_blank())
p <- cowplot::plot_grid(pA, pB, nrow=2, align="v", axis="lr")
ggsave("figs/pub_new/ensBlend_eval_dotplot2.png", p, width=9, height=6, dpi=400)




pA <- point_df |>
  filter(type=="Mean within-week") |>
  left_join(const_z, by=join_by(nSim, sample, metric, type)) |>
  mutate(skill_z=(skill - mn)/sd,
         nSim=as.numeric(nSim)) |>
  filter(abs(skill_z) <= 6) |>
  ggplot(aes(skill_z, nSim, fill=modType, colour=modType)) +
  geom_vline(xintercept=0, linewidth=0.5, colour="grey85") +
  ggridges::geom_density_ridges(aes(group=paste(nSim, modType)), scale=0.7, 
                                fill=NA, show.legend=F, rel_min_height=0.001) +
  stat_pointinterval(aes(y=nSim-0.15, group=paste(nSim, modType)), 
                     .width=c(0.8), linewidth=0.5,
                     position=position_dodge(width=0.25), shape=1, size=0.75) + 
  scale_fill_manual("Model type", values=modType_cols,
                    labels=lab_expressions) +
  scale_colour_manual("Model type", values=modType_cols,
                      labels=lab_expressions) +
  scale_x_continuous(" ") +
  scale_y_continuous("Number of constituents", 
                     breaks=1:4, labels=c("3", "5", "10", "15")) +
  facet_grid(type~metric, scales="free_x", labeller=labeller(metric=label_parsed)) + 
  theme_bw() + 
  theme(legend.position="inside",
        legend.position.inside=c(0.94, 0.89),
        legend.title=element_blank(),
        legend.text=element_text(size=7),
        legend.key.height=unit(2.5, "mm"),
        legend.key.width=unit(2.5, "mm"),
        legend.background=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.y=element_blank())
pB <- point_df |>
  filter(type=="Mean within-farm") |>
  left_join(const_z, by=join_by(nSim, sample, metric, type)) |>
  mutate(skill_z=(skill - mn)/sd,
         nSim=as.numeric(nSim)) |>
  filter(abs(skill_z) <= 6) |>
  ggplot(aes(skill_z, nSim, fill=modType, colour=modType)) +
  geom_vline(xintercept=0, linewidth=0.5, colour="grey85") +
  ggridges::geom_density_ridges(aes(group=paste(nSim, modType)), scale=0.7, 
                                fill=NA, show.legend=F, rel_min_height=0.001) +
  stat_pointinterval(aes(y=nSim-0.15, group=paste(nSim, modType)), 
                     .width=c(0.8), linewidth=0.5,
                     position=position_dodge(width=0.25), shape=1, size=0.75) + 
  scale_fill_manual("Model type", values=modType_cols,
                    labels=lab_expressions) +
  scale_colour_manual("Model type", values=modType_cols,
                      labels=lab_expressions) +
  scale_x_continuous("Z-score of skill within resample") +
  scale_y_continuous("Number of constituents", 
                     breaks=1:4, labels=c("3", "5", "10", "15")) +
  facet_grid(type~metric, scales="free_x", labeller=labeller(metric=label_parsed)) + 
  theme_bw() + 
  theme(legend.position="none",
        legend.background=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.y=element_blank())
p <- cowplot::plot_grid(pA, pB, nrow=2, align="v", axis="lr")
ggsave("figs/pub_new/ensBlend_eval_ridges.png", p, width=9, height=6, dpi=400)



point_df |>
  group_by(type, nSim, metric, modType, sample) |>
  slice_max(skill) |>
  ungroup() |>
  left_join(const_z, by=join_by(nSim, sample, metric, type)) |>
  mutate(skill_z=(skill - mn)/sd) |>
  group_by(type, nSim, metric, modType) |>
  summarise(mn=mean(skill_z),
            md=median(skill_z),
            sd=sd(skill_z)) |>
  arrange(modType, metric, nSim) |>
  ggplot(aes(nSim, mn, ymin=mn-sd, ymax=mn+sd, colour=modType)) + 
  geom_line(aes(group=modType)) +
  geom_point(position=position_dodge(width=0.2)) + 
  geom_linerange(position=position_dodge(width=0.2)) + 
  scale_colour_manual("Model type", values=modType_cols,
                      labels=lab_expressions) +
  facet_grid(type~metric)






# Blending proportions ------------------------------------------------------

mod <- "n20_D4"
out_ensBlend <- readRDS(glue("out/ensembles/ensBlend_{mod}_FULL_stanfit.rds"))
dat_ensBlend <- readRDS(glue("out/ensembles/ensBlend_{mod}_FULL_standata.rds"))
ensFull_LatLon <- read_csv("out/valid_df_2021-2024.csv") |>
  select(rowNum, date, CV_k, sepaSite, sepaSiteNum, licePerFish_rtrt, starts_with("sim")) |>
  select(-contains("avg")) |>
  mutate(across(starts_with("sim_"), ~.x - mean(.x), .names="c_{.col}")) |>
  left_join(site_i) |>
  select(-sepaSite) |>
  arrange(rowNum)
ensBlend_rec <- make_spline_recipe(ensFull_LatLon, 
                                   as.numeric(str_split_fixed(mod, "D", 2)[,2]), 
                                   sim_i$sim[1:20])
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_bbox <- st_bbox(mesh_fp)
mesh_land <- st_convex_hull(mesh_fp) |>
  st_difference(mesh_fp) |>
  st_crop(site_i |> st_as_sf(coords=c("easting", "northing"), crs=27700) |> st_buffer(10e3))


# blending proportions by site
map_df <- site_i |>
  mutate(sepaSiteNum=row_number()) |>
  bind_cols(ensFull_LatLon |> summarise(across(c(licePerFish_rtrt, contains("sim")), mean))) 
b_p_ls <- rstan::extract(out_ensBlend, "b_p")[[1]]
# b_p_ls <- make_predictions_ensBlend_sLonLat(out_ensBlend, 
#                                             newdata=bake(ensBlend_rec, map_df), 
#                                             iter=3000, mode="b_p") 
site_b_p_post <- map_dfr(1:dim(b_p_ls)[1], 
                      ~as_tibble(b_p_ls[.x,,]) |>
                        set_names(dat_ensBlend$sim_names) |>
                        mutate(rowNum=row_number(),
                               iter=.x)) |> 
  pivot_longer(starts_with("sim"), names_to="sim", values_to="p") |>
  left_join(sim_i) |>
  rename(Simulation=lab)
site_b_p_mns <- site_b_p_post |>
  group_by(rowNum, Simulation) |>
  summarise(p_mn=mean(p)) |>
  ungroup() |>
  rename(sepaSiteNum=rowNum) |>
  inner_join(ensFull_df |> group_by(sepaSite) |> slice_head() |> select(sepaSiteNum, sepaSite)) |>
  inner_join(site_i)
gc()
  
p_a <- site_b_p_post |> 
  ggplot(aes(p, group=rowNum)) + 
  geom_line(alpha=0.2, linewidth=0.25, stat="density", adjust=2) +
  labs(x=expression(paste("Ensemble blending weight (", italic(pi[~~k]), ")")),
       y="log density") +
  scale_y_continuous(transform="log1p") +
  facet_wrap(~Simulation, labeller=label_parsed, scales="free_y", ncol=2, strip.position="right") +
  theme(axis.title=element_text(size=9),
        axis.text=element_text(size=7),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave(glue("figs/pub_new/ensBlend_{mod}_p_sitePosterior.png"), p_a, width=8, height=8)



# blending proportions by parameter: maps
map_df <- expand_grid(easting=seq(min(site_i$easting)-30e3, max(site_i$easting)+30e3, by=4e3),
                      northing=seq(min(site_i$northing)-30e3, max(site_i$northing)+30e3, by=4e3)) |>
  st_as_sf(coords=c("easting", "northing"), crs=27700, remove=F) |>
  st_intersection(st_buffer(mesh_fp, 5e3)) |>
  st_intersection(site_i |> 
                    st_as_sf(coords=c("easting", "northing"), crs=27700) |> 
                    st_buffer(50e3) |> 
                    st_union()) |>
  st_drop_geometry() |>
  mutate(sepaSiteNum=row_number()) |>
  bind_cols(ensFull_LatLon |> summarise(across(c(licePerFish_rtrt, contains("sim")), mean))) 
b_p_ls <- make_predictions_ensBlend_sLonLat(out_ensBlend, 
                                            newdata=bake(ensBlend_rec, map_df), 
                                            iter=1000, mode="b_p") 
b_p_post <- map_dfr(1:dim(b_p_ls)[2], 
                    ~as_tibble(b_p_ls[,.x,]) |>
                      set_names(dat_ensBlend$sim_names) |>
                      mutate(rowNum=row_number(),
                             iter=.x)) |> 
  pivot_longer(starts_with("sim"), names_to="sim", values_to="p") |>
  left_join(sim_i) |>
  rename(Simulation=lab)

sim_p_map_df <- b_p_post |>
  group_by(Simulation, rowNum) |>
  summarise(p_mn=mean(p),
            p_sd=sd(p)) |>
  full_join(map_df |> select(easting, northing) |> mutate(rowNum=row_number()),
            by=join_by(rowNum)) |>
  mutate(log_p=log10(p_mn))

p <- ggplot(sim_p_map_df) + 
  geom_raster(aes(easting, northing, fill=p_mn)) +
  stat_contour(aes(easting, northing, z=p_mn), colour="white", linewidth=0.1) +
  colorspace::scale_fill_continuous_sequential(name="Posterior mean weight (p)",
                                               palette="GnBu",
                                               rev=T,
                                               limits=c(0, 1),
                                               breaks=c(0, 0.5, 1)) +
  geom_sf(data=mesh_land, fill="grey90", colour="grey40", linewidth=0.1) +
  geom_point(data=site_b_p_mns, aes(easting, northing, fill=p_mn),
             shape=21, colour="black", stroke=0.1) +
  scale_x_continuous(limits=range(site_i$easting)*c(0.96, 1), oob=scales::oob_keep) +
  scale_y_continuous(limits=range(site_i$northing), oob=scales::oob_keep) +
  facet_wrap(~Simulation, nrow=4, labeller=label_parsed) +
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        legend.position="bottom",
        legend.title.position="top",
        legend.title=element_text(size=9, hjust=0.5),
        legend.box.margin=margin(0,0,0,0),
        legend.margin=margin(0,0,0,0),
        legend.key.height=unit(0.2, "cm"),
        legend.key.width=unit(0.8, "cm"),
        legend.text=element_text(size=6),
        panel.spacing=unit(0.1, 'cm'))
ggsave(glue("figs/pub_new/ensBlend_{mod}_pSim_map.png"), p, height=12, width=7, dpi=200)


p <- ggplot(sim_p_map_df) + 
  geom_raster(aes(easting, northing, fill=log_p)) +
  stat_contour(aes(easting, northing, z=log_p), colour="white", linewidth=0.1) +
  colorspace::scale_fill_continuous_sequential(name="Posterior mean weight (p)",
                                               palette="GnBu",
                                               rev=T,
                                               limits=c(NA, 0),
                                               breaks=c(-3, -2, -1, 0),
                                               labels=10^(c(-3, -2, -1, 0))) +
  geom_sf(data=mesh_land, fill="grey90", colour="grey40", linewidth=0.1) +
  geom_point(data=site_b_p_mns, aes(easting, northing, fill=log10(p_mn)),
             shape=21, colour="black", stroke=0.1) +
  scale_x_continuous(limits=range(site_i$easting)*c(0.96, 1), oob=scales::oob_keep) +
  scale_y_continuous(limits=range(site_i$northing), oob=scales::oob_keep) +
  facet_wrap(~Simulation, nrow=4, labeller=label_parsed) +
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        legend.position="bottom",
        legend.title.position="top",
        legend.title=element_text(size=9, hjust=0.5),
        legend.box.margin=margin(0,0,0,0),
        legend.margin=margin(0,0,0,0),
        legend.key.height=unit(0.2, "cm"),
        legend.key.width=unit(0.8, "cm"),
        legend.text=element_text(size=6),
        panel.spacing=unit(0.1, 'cm'))
ggsave(glue("figs/pub_new/ensBlend_{mod}_pSim_map-log10.png"), p, height=12, width=7, dpi=200)


selectedSims <-  site_b_p_mns |> select(Simulation, p_mn) |>
  group_by(Simulation) |>
  slice_max(p_mn) |>
  ungroup() |>
  filter(p_mn > 0.1)
  
p <- ggplot(sim_p_map_df |> filter(Simulation %in% selectedSims$Simulation)) + 
  geom_raster(aes(easting, northing, fill=p_mn)) +
  stat_contour(aes(easting, northing, z=p_mn), colour="grey30", alpha=0.5,
               linewidth=0.1, breaks=seq(0, 1, by=0.1)) +
  colorspace::scale_fill_continuous_sequential(name="Posterior mean\nblending weight (p)",
                                               palette="GnBu",
                                               rev=T,
                                               limits=c(0, 1),
                                               breaks=seq(0, 1, by=0.1),
                                               labels=c("0", rep("", 4), "0.5", rep("", 4), "1")) +
  geom_sf(data=mesh_land, fill="grey90", colour="grey40", linewidth=0.1) +
  geom_point(data=site_b_p_mns |> filter(Simulation %in% selectedSims$Simulation), 
             aes(easting, northing, fill=p_mn),
             shape=21, colour="black", stroke=0.1) +
  scale_x_continuous(limits=range(site_i$easting)*c(0.96, 1), oob=scales::oob_keep) +
  scale_y_continuous(limits=range(site_i$northing), oob=scales::oob_keep) +
  facet_wrap(~Simulation, nrow=2, labeller=label_parsed) +
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        legend.position="inside",
        legend.position.inside=c(0.85, 0.25),
        legend.title.position="top",
        legend.title=element_text(size=9, hjust=0),
        legend.box.margin=margin(0,0,0,0),
        legend.margin=margin(0,0,0,0),
        legend.key.height=unit(0.8, "cm"),
        legend.key.width=unit(0.2, "cm"),
        legend.text=element_text(size=6),
        legend.ticks=element_line(colour="grey30", linewidth=0.1),
        panel.spacing=unit(0.1, 'cm'))
ggsave(glue("figs/pub_new/ensBlend_{mod}_pSim_map_selected.png"), p, height=6.5, width=4.5, dpi=200)



sim_key <- read_csv("out/sim_2021-2024/sim_i.csv") |> 
  mutate(sim=paste0("sim_", i)) |>
  select(-outDir) |> 
  inner_join(sim_i) |>
  mutate(across(matches("swim|Thresh"), ~if_else(lab_short=="2D", NA, .x)),
         across(matches("swim"), ~abs(.x)*100),
         fixDepth=as.numeric(fixDepth),
         mortSal_fn=as.numeric(mortSal_fn=="logistic"),
         eggTemp_fn=as.numeric(eggTemp_fn=="logistic")) |>
  rename(Simulation=lab)

param_post_sum <- sim_key |> 
  select(-i, -sim, -lab_short, -Simulation) |> 
  names() |>
  map(~summarise_param_posterior(b_p_post, sim_key, .x)) |>
  reduce(full_join, by=join_by(iter, rowNum))
gc()


# param_post_sum |>
#   pivot_longer(-(1:2), names_to="var", values_to="value") |>
#   ggplot(aes(value)) + 
#   geom_density() + 
#   facet_wrap(~var, scales="free")

param_map_ls <- param_post_sum |>
  pivot_longer(-(1:2), names_to="var", values_to="value") |>
  group_by(rowNum, var) |>
  summarise(post_mn=mean(value),
            post_sd=sd(value)) |>
  ungroup() |>
  inner_join(var_pretty, by=join_by(var)) |>
  arrange(var_order) |>
  full_join(map_df |> select(easting, northing) |> mutate(rowNum=row_number()),
            by=join_by(rowNum)) |>
  group_split(var_order, .keep=TRUE)


param_map_plot_ls <- map(param_map_ls, ~make_param_map_plot(.x, site_i, mesh_land, sim_key))
p <- plot_grid(plotlist=param_map_plot_ls, align="hv", axis="tblr", nrow=2)
ggsave(glue("figs/pub_new/ensBlend_{mod}_pParam_map.png"), p, height=10.2, width=12, dpi=200)




# parameterization performance --------------------------------------------

sim_params <- read_csv("out/sim_2021-2024/sim_i.csv") |> 
  mutate(across(matches("Thresh|swim"), ~if_else(!fixDepth, .x, NA)),
         fixDepth=as.numeric(fixDepth), 
         across(ends_with("fn"), ~as.numeric(if_else(.x=="constant", 0, 1))),
         across(contains("SpeedMean"), ~abs(.x)),
         dSalinity=salinityThreshMax - salinityThreshMin,
         sim=paste0("sim_", i)) |> 
  pivot_longer(-any_of(c("i", "outDir", "sim")), 
               names_to="param_name", values_to="param_val")

mod <- "n20_D4"
out_ensBlend <- readRDS(glue("out/ensembles/ensBlend_{mod}_FULL_stanfit.rds"))
dat_ensBlend <- readRDS(glue("out/ensembles/ensBlend_{mod}_FULL_standata.rds"))

# blending proportions by site
b_p_ls <- rstan::extract(out_ensBlend, "b_p")[[1]]
b_p_post <- map_dfr(1:dim(b_p_ls)[1], 
                    ~as_tibble(b_p_ls[.x,,]) |>
                      set_names(dat_ensBlend$sim_names) |>
                      mutate(rowNum=row_number(),
                             iter=.x)) |> 
  pivot_longer(starts_with("sim"), names_to="sim", values_to="p") |>
  left_join(sim_i) |>
  rename(Simulation=lab)
p_summary <- b_p_post |>
  group_by(sim, Simulation, lab_short) |>
  sevcheck::get_intervals(p) |>
  full_join(sim_params, by="sim")

sim_params |> 
  ggplot(aes(param_val, i, colour=i %in% c("07", "03", "11"))) +
  # geom_hline(yintercept=c("07", "03", "11"), colour="grey") +
  geom_text(aes(label=i)) + 
  facet_wrap(~param_name, scales="free", nrow=3)
sim_params |> 
  group_by(param_name) |>
  mutate(rank=rank(param_val)) |>
  ggplot(aes(rank, param_val, colour=i %in% c("07", "03", "11"))) +
  geom_text(aes(label=i)) + 
  facet_wrap(~param_name, scales="free", nrow=3)

p_summary |> 
  group_by(param_name) |>
  mutate(rank=rank(param_val)) |>
  ggplot(aes(rank, param_val, colour=mn)) +
  geom_text(aes(label=i)) + 
  scale_colour_distiller(palette="Reds", direction=1) +
  facet_wrap(~param_name, scales="free", nrow=3)

p_summary |> 
  group_by(param_name) |>
  mutate(rank=rank(param_val)) |>
  ggplot(aes(mn, rank)) +
  geom_text(aes(label=i)) + 
  facet_wrap(~param_name, scales="free", nrow=3)


param_metrics_df <- inner_join(sim_params, 
           all_metrics_df, 
           by=join_by(sim), relationship="many-to-many") |>
  filter(! metric %in% c("R^2", "r"))

param_metrics_df |>
  ggplot(aes(param_val, value, colour=type)) + 
  geom_smooth(method="loess", span=2, se=F, linewidth=0.5, linetype=3) +
  geom_point() +
  labs(x="Parameter value", y="Cross-validation score") +
  facet_grid(metric~param_name, scales="free")

param_metrics_df |>
  filter(metric=="RMSE") |>
  ggplot(aes(param_val, value, colour=type)) + 
  geom_smooth(method="loess", span=2, se=F, linewidth=0.5, linetype=3) +
  geom_point() +
  facet_wrap(~param_name, scales="free_x")
param_metrics_df |>
  filter(metric=="'AUC'['ROC']") |>
  ggplot(aes(param_val, value, colour=type)) + 
  geom_smooth(method="loess", span=2, se=F, linewidth=0.5, linetype=3) +
  geom_point() +
  facet_wrap(~param_name, scales="free_x")
param_metrics_df |>
  filter(metric=="rho") |>
  ggplot(aes(param_val, value, colour=type)) + 
  geom_smooth(method="loess", span=2, se=F, linewidth=0.5, linetype=3) +
  geom_point() +
  facet_wrap(~param_name, scales="free_x")
param_metrics_df |>
  filter(metric=="'AUC'['PR']") |>
  ggplot(aes(param_val, value, colour=type)) + 
  geom_smooth(method="loess", span=2, se=F, linewidth=0.5, linetype=3) +
  geom_point() +
  facet_wrap(~param_name, scales="free_x")



library(randomForest)
library(vip)
RMSE_df <- param_metrics_df |>
  filter(metric=="RMSE") |>
  select(value, param_name, param_val) |>
  mutate(param_val=if_else(is.na(param_val), 0, param_val)) |>
  pivot_wider(names_from=param_name, values_from=param_val) |>
  rename(RMSE=value) 
rf_RMSE <- randomForest(RMSE ~ ., data=RMSE_df, importance=T)
varImpPlot(rf_RMSE)

rho_df <- param_metrics_df |>
  filter(metric=="rho") |>
  select(value, param_name, param_val) |>
  mutate(param_val=if_else(is.na(param_val), 0, param_val)) |>
  pivot_wider(names_from=param_name, values_from=param_val) |>
  rename(rho=value)
rf_rho <- randomForest(rho ~ ., data=rho_df, importance=T)
varImpPlot(rf_rho)

ROC_df <- param_metrics_df |>
  filter(metric=="'AUC'['ROC']") |>
  select(sim, value, param_name, param_val) |>
  mutate(param_val=if_else(is.na(param_val), 0, param_val)) |>
  pivot_wider(names_from=param_name, values_from=param_val) |>
  rename(ROC=value) |>
  select(-sim)
rf_ROC <- randomForest(ROC ~ ., data=ROC_df, importance=T)
varImpPlot(rf_ROC)

PR_df <- param_metrics_df |>
  filter(metric=="'AUC'['PR']") |>
  select(sim, value, param_name, param_val) |>
  mutate(param_val=if_else(is.na(param_val), 0, param_val)) |>
  pivot_wider(names_from=param_name, values_from=param_val) |>
  rename(PR=value) |>
  select(-sim)
rf_PR <- randomForest(PR ~ ., data=PR_df, importance=T)
varImpPlot(rf_PR)





# ensFcst explainers ------------------------------------------------------

library(DALEX); library(DALEXtra); library(tidymodels)

wf_fitted <- readRDS("out/ensembles/licePerFish_best_fitted_1wk_rmse.rds")
wf_explain <- explain_tidymodels(wf_fitted, 
                                 data=ensFull_df |> left_join(site_i), 
                                 y=ensFull_df$licePerFish_rtrt)
plot(model_parts(wf_explain))

mods <- c(paste0("sim_", str_pad(1:20, 2, "left", "0")))
acdep_df <- model_profile(wf_explain, variables=mods, type="accumulated", center=F)

acdep_df$agr_profiles |>
  left_join(ensFull_df |>
              summarise(across(starts_with("sim"), mean)) |>
              pivot_longer(everything(), names_to="_vname_", values_to="_mean_")) |>
  mutate(`_x_`=`_x_` + `_mean_`) |>
  ggplot(aes(`_x_`, `_yhat_`^4, colour=`_vname_`)) + 
  geom_line() + 
  # ylim(0, NA) +
  scale_colour_viridis_d(option="turbo", begin=0.1) +
  facet_wrap(~`_vname_`) +
  labs(x=expression("AEIP"), y="Conditional mean predicted lice per fish")
  



# scatterplots ------------------------------------------------------------

ensCV_df <- read_csv("out/ensemble_CV.csv")

preds_df <- ensCV_df |>
  select(rowNum, licePerFish_rtrt, lice_g05, IP_sim_avg2D, IP_sim_avg3D,
         IP_predBlend, IP_predFwk_RMSE) |>
  pivot_longer(starts_with("IP")) |>
  mutate(name=factor(name, 
                     levels=paste0("IP_", c("null", "sim_avg2D", "sim_avg3D", 
                                            "predFwk_RMSE", "predBlend")),
                     labels=c("Null", "Mean['2D']", "Mean['3D']",
                              "Ens['Fcst']",  "Ens['Blend']"))) 
p <- preds_df |>
  filter(licePerFish_rtrt > 0) |>
  ggplot(aes(value, licePerFish_rtrt)) + 
  geom_abline(colour="grey") + 
  geom_hline(colour="grey", linetype=3, yintercept=0.5^0.25) +
  geom_vline(colour="grey", linetype=3, xintercept=0.5^0.25) +
  stat_smooth(data=preds_df, linetype=2, linewidth=0.5, alpha=0.5,
              colour="cadetblue", fill="cadetblue") +
  geom_point(size=0.75, shape=1, alpha=0.1) +
  geom_jitter(data=preds_df |> filter(licePerFish_rtrt==0), 
              size=0.75, shape=1, alpha=0.1, width=0, height=0.025) + 
  coord_equal() +
  facet_wrap(~name, labeller=label_parsed, ncol=2) +
  labs(x=expression(paste("Predicted (Mean"~~italic("L. salmonis")~~"per fish)"^0.25)),
       y=expression(paste("(Mean"~~italic("L. salmonis")~~"per fish)"^0.25)))
ggsave("figs/pub_new/predictions_CV_scatterplot.png", p, width=5, height=7)


walk(unique(ensCV_df$sepaSite), 
     ~(ensCV_df |>
         filter(sepaSite==.x) |>
         pivot_longer(starts_with("IP_")) |> 
         filter(name %in% c("IP_sim_07", "IP_predBlend", "IP_sim_avgAll",
                            "IP_predFf_RMSE", "IP_predFf_rsq")) |>
         ggplot(aes(value, licePerFish_rtrt)) + 
         geom_abline() + 
         geom_hline(yintercept=0.5^0.25, linetype=2) + 
         geom_vline(xintercept=0.5^0.25, linetype=2) + 
         geom_point(alpha=0.5, shape=1) + 
         geom_line(stat="smooth", method="lm", formula=y~x, se=F, colour="dodgerblue") +
         facet_wrap(~name, nrow=2) + 
         xlim(0, 2) + ylim(0, 2.25) + coord_equal()) |>
       ggsave(glue("figs/siteScatter/{.x}.png"), plot=_, width=8, height=6))

walk(unique(ensCV_df$date), 
     ~(ensCV_df |>
         filter(date==.x) |>
         pivot_longer(starts_with("IP_")) |> 
         filter(name %in% c("IP_sim_07", "IP_predBlend", "IP_sim_avgAll",
                            "IP_predFwk_RMSE", "IP_predFwk_rsq")) |>
         ggplot(aes(value, licePerFish_rtrt)) + 
         geom_abline() + 
         geom_hline(yintercept=0.5^0.25, linetype=2) + 
         geom_vline(xintercept=0.5^0.25, linetype=2) + 
         geom_point(alpha=0.5, shape=1) + 
         geom_line(stat="smooth", method="lm", formula=y~x, se=F, colour="dodgerblue") +
         facet_wrap(~name, nrow=2) + 
         xlim(0, 2) + ylim(0, 2.25) + coord_equal()) |>
       ggsave(glue("figs/dateScatter/{.x}.png"), plot=_, width=5, height=7))


walk(unique(ensCV_df$sepaSite), 
     ~(ensCV_df |>
         filter(sepaSite==.x) |>
         pivot_longer(starts_with("IP_")) |> 
         filter(name %in% c("IP_sim_07", "IP_predBlend", "IP_sim_avgAll",
                            "IP_predFwk_ROC", "IP_predFwk_PR")) |>
         ggplot(aes(value, as.numeric(lice_g05))) +
         stat_smooth(method="glm", method.args=list(family="binomial"), 
                   formula=y~x, se=T, colour="dodgerblue", fullrange=T) +
         geom_dots(aes(side=lice_g05), scale=0.4) +
         scale_side_mirrored(guide="none") +
         coord_cartesian(ylim = c(0, 1)) +
         facet_wrap(~name, nrow=2)) |>
       ggsave(glue("figs/sitePr/{.x}.png"), plot=_, width=7, height=5))

walk(unique(ensCV_df$date), 
     ~(ensCV_df |>
         filter(date==.x) |>
         pivot_longer(starts_with("IP_")) |> 
         filter(name %in% c("IP_sim_07", "IP_predBlend", "IP_sim_avgAll",
                            "IP_predFwk_ROC", "IP_predFwk_PR")) |>
         ggplot(aes(value, as.numeric(lice_g05))) +
         stat_smooth(method="glm", method.args=list(family="binomial"), 
                     formula=y~x, se=T, colour="dodgerblue", fullrange=T) +
         geom_dots(aes(side=lice_g05), scale=0.4) +
         scale_side_mirrored(guide="none") +
         coord_cartesian(ylim = c(0, 1)) +
         facet_wrap(~name, nrow=2)) |>
       ggsave(glue("figs/datePr/{.x}.png"), plot=_, width=7, height=5))



# IP sLL ------------------------------------------------------------------

set.seed(1003)
mod <- "n20_D4"
out_ensBlend <- readRDS(glue("out/ensembles/ensBlend_{mod}_FULL_stanfit.rds"))
dat_ensBlend <- readRDS(glue("out/ensembles/ensBlend_{mod}_FULL_standata.rds"))
IP_LatLon <- read_csv("out/valid_df_2021-2024.csv") |>
  select(rowNum, date, CV_k, sepaSite, sepaSiteNum, licePerFish_rtrt, starts_with("sim")) |>
  select(-contains("avg")) |>
  mutate(across(starts_with("sim_"), ~.x - mean(.x), .names="c_{.col}")) |>
  group_by(sepaSite, sepaSiteNum) |>
  summarise(rowNum=first(rowNum),
            date=first(date),
            CV_k=first(CV_k),
            licePerFish_rtrt=mean(licePerFish_rtrt),
            across(contains("sim"), ~mean(.x, na.rm=T))) |>
  left_join(site_i) |>
  arrange(sepaSiteNum)

site_LU <- IP_LatLon |> group_by(sepaSite) |> slice_head(n=1) |> ungroup() |> select(contains("sepa"))
b_p_ls <- rstan::extract(out_ensBlend, "b_p")[[1]]
site_p_post <- map_dfr(1:dim(b_p_ls)[1], 
                       ~as_tibble(b_p_ls[.x,,]) |>
                         set_names(dat_ensBlend$sim_names) |>
                         mutate(sepaSiteNum=row_number(),
                                iter=.x)) |> 
  pivot_longer(starts_with("sim"), names_to="sim", values_to="p") |>
  left_join(site_LU, by=join_by(sepaSiteNum)) |>
  nest(p=c(iter, p))
gc()
rm(b_p_ls); rm(out_ensBlend); gc()

influx_df <- readRDS("out/sim_2021-2024/processed/connectivity_day.rds") |>
  select(sepaSite, date, sim, influx_m3) |>
  mutate(sim=paste0("sim_", sim),
         influx_m3_4rt=replace_na(influx_m3, 0)^0.25)

date_seq <- sort(unique(influx_df$date))
ens_ls <- vector("list", length(date_seq))
for(i in seq_along(date_seq)) {
  ens_ls[[i]] <- influx_df |>
    filter(date==date_seq[i]) |>
    inner_join(site_p_post, by=join_by("sim", "sepaSite")) |>
    unnest(p) |>
    mutate(wtIP=influx_m3_4rt * p) |>
    group_by(sepaSite, date, iter) |>
    summarise(ens_IP=sum(wtIP), .groups="keep") |>
    group_by(sepaSite, date) |>
    summarise(lice_mn=mean(ens_IP),
              lice_q005=quantile(ens_IP, probs=0.005),
              lice_q025=quantile(ens_IP, probs=0.025),
              lice_q975=quantile(ens_IP, probs=0.975),
              lice_q995=quantile(ens_IP, probs=0.995),
              .groups="keep") |>
    ungroup() |>
    mutate(across(starts_with("lice_"), ~.x^4))
  if(i %% 14 == 0) {
    cat("Finished", as.character(date_seq[i]), "\n")
  }
  gc()
}
ens_ls |>
  reduce(bind_rows) |>
  saveRDS("out/sim_2021-2024/processed/influx_ens.rds")


influx_ens <- readRDS("out/sim_2021-2024/processed/influx_ens.rds")

thresholds <- c(0, 1e-4, 1e-3, 1e-2, 1e-1, 1)

fig_influx <- influx_ens |>
  filter(date >= "2021-05-01") |>
  group_by(date) |>
  rename(lice=lice_mn) |>
  summarise(lt_t1=mean(lice==0),
            lt_t2=mean(lice > thresholds[1] & lice < thresholds[2]),
            lt_t3=mean(between(lice, thresholds[2], thresholds[3])),
            lt_t4=mean(between(lice, thresholds[3], thresholds[4])),
            lt_t5=mean(between(lice, thresholds[4], thresholds[5])),
            lt_t6=mean(between(lice, thresholds[5], thresholds[6])),
            lt_t7=mean(lice > thresholds[6])) |>
  ungroup() |>
  pivot_longer(starts_with("lt_"), names_to="threshold", values_to="propSites") |>
  mutate(threshold=factor(threshold, 
                          labels=c("0", 
                                   paste(thresholds[1:5], "-", thresholds[2:6]),
                                   paste(">", thresholds[6]))),
         threshold_num=as.numeric(threshold)) |>
  filter(threshold_num != 1) |>
  ggplot(aes(date, propSites, fill=threshold_num, group=threshold_num)) +
  geom_hline(yintercept=c(0, 1), colour="grey", linewidth=0.2) +
  geom_area(colour="grey30", linewidth=0.05, outline.type="both") +
  scale_x_date(date_breaks="1 year", 
               date_labels="%Y", expand=expansion(mult=c(0.05, 0.05))) +
  scale_y_continuous("Proportion of active farms", limits=c(0,1)) +
  scale_fill_viridis_b(expression(paste("Ensemble mean daily copepodids" %.% "m"^"-3" %.% "h"^"-1")),
                       option="turbo", begin=0.05,
                       breaks=c(2.5, 3.5, 4.5, 5.5, 6.5),
                       labels=c("0.0001", "0.001", "0.01", "0.1", "1")) +
  theme(panel.grid.major.x=element_line(colour="grey", linewidth=0.6),
        panel.grid.minor.x=element_line(colour="grey", linewidth=0.2),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=11),
        legend.position="bottom", 
        legend.key.width=unit(1.5, "cm"), 
        legend.key.height=unit(0.2, "cm"),
        strip.text=element_text(size=11))
ggsave("figs/pub_new/ensBlend_influx_daily.png", fig_influx, width=7, height=4, dpi=400)



influx_ens |>
  filter(date >= "2021-05-01") |>
  group_by(date) |>
  rename(lice=lice_mn) |>
  summarise(lt_t1=mean(lice==0),
            lt_t2=mean(lice > thresholds[1] & lice < thresholds[2]),
            lt_t3=mean(between(lice, thresholds[2], thresholds[3])),
            lt_t4=mean(between(lice, thresholds[3], thresholds[4])),
            lt_t5=mean(between(lice, thresholds[4], thresholds[5])),
            lt_t6=mean(between(lice, thresholds[5], thresholds[6])),
            lt_t7=mean(lice > thresholds[6])) |>
  ungroup() |>
  pivot_longer(starts_with("lt_"), names_to="threshold", values_to="propSites") |>
  mutate(threshold=factor(threshold, 
                          labels=c("0", 
                                   paste(thresholds[1:5], "-", thresholds[2:6]),
                                   paste(">", thresholds[6]))),
         threshold_num=as.numeric(threshold)) |>
  filter(threshold_num != 1) |>
  ggplot(aes(date, propSites, fill=threshold_num, group=threshold_num)) +
  geom_hline(yintercept=0, colour="grey", linewidth=0.2) +
  geom_area(colour="grey30", linewidth=0.05, outline.type="both") +
  scale_x_date(date_breaks="1 year", 
               date_labels="%Y", expand=expansion(mult=c(0.05, 0.05))) +
  scale_y_continuous("Proportion of active farms") +
  scale_fill_viridis_b(expression(paste("Ensemble mean daily copepodids" %.% "m"^"-3" %.% "h"^"-1")),
                       option="turbo", begin=0.05,
                       breaks=c(2.5, 3.5, 4.5, 5.5, 6.5),
                       labels=c("0.0001", "0.001", "0.01", "0.1", "1")) +
  facet_grid(threshold~.) +
  theme(panel.grid.major.x=element_line(colour="grey", linewidth=0.6),
        panel.grid.minor.x=element_line(colour="grey", linewidth=0.2),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=11),
        legend.position="bottom", 
        legend.key.width=unit(1.5, "cm"), 
        legend.key.height=unit(0.2, "cm"),
        strip.text=element_text(size=11))


thresholds <- c(0, 1e-3, 1e-2, 1e-1, 1)

# fig_influx <- influx_ens |>
#   group_by(date) |>
#   rename(lice=lice_mn) |>
#   summarise(lt_t1=mean(lice==0),
#             lt_t2=mean(lice > thresholds[1] & lice < thresholds[2]),
#             lt_t3=mean(between(lice, thresholds[2], thresholds[3])),
#             lt_t4=mean(between(lice, thresholds[3], thresholds[4])),
#             lt_t5=mean(between(lice, thresholds[4], thresholds[5])),
#             lt_t6=mean(lice > thresholds[5])) |>
#   ungroup() |>
#   pivot_longer(starts_with("lt_"), names_to="threshold", values_to="propSites") |>
#   mutate(threshold=factor(threshold, 
#                           labels=c("0", 
#                                    paste(thresholds[1:4], "-", thresholds[2:5]),
#                                    paste(">", thresholds[5]))),
#          threshold_num=as.numeric(threshold)) |>
#   filter(threshold_num != 1) |>
#   ggplot(aes(date, propSites, fill=threshold_num, group=threshold_num)) +
#   geom_area(colour="grey30", linewidth=0.05, outline.type="both") +
#   scale_x_date(date_breaks="1 year", date_minor_breaks="3 months", date_labels="%Y") +
#   scale_y_continuous("Proportion of active farms", limits=c(0,1)) +
#   scale_fill_viridis_b(expression(paste("Copepodids" %.% "m"^"-2" %.% "h"^"-1")),
#                        option="turbo", begin=0.05,
#                        breaks=c(2.5, 3.5, 4.5, 5.5),
#                        labels=c("0.001", "0.01", "0.1", "1")) +
#   theme(panel.grid.major.x=element_line(colour="grey", linewidth=0.6),
#         panel.grid.minor.x=element_line(colour="grey", linewidth=0.2),
#         axis.title.x=element_blank(),
#         axis.title.y=element_text(size=14),
#         legend.position="bottom", 
#         legend.title=element_text(size=12),
#         legend.key.width=unit(1.5, "cm"), 
#         legend.key.height=unit(0.2, "cm"),
#         strip.text=element_text(size=14))
# ggsave("figs/talk/ens_influx_daily_sLonLatD4.png", fig_influx, width=6, height=4, dpi=400)




# Tidal cycles?
fig_influx <- influx_ens |>
  filter(date >= "2021-05-01") |>
  group_by(date) |>
  rename(lice=lice_mn) |>
  summarise(lt_t1=mean(lice==0),
            lt_t2=mean(lice > thresholds[1] & lice < thresholds[2]),
            lt_t3=mean(between(lice, thresholds[2], thresholds[3])),
            lt_t4=mean(between(lice, thresholds[3], thresholds[4])),
            lt_t5=mean(between(lice, thresholds[4], thresholds[5])),
            lt_t6=mean(between(lice, thresholds[5], thresholds[6])),
            lt_t7=mean(lice > thresholds[6])) |>
  ungroup() |>
  pivot_longer(starts_with("lt_"), names_to="threshold", values_to="propSites") |>
  mutate(threshold=factor(threshold, 
                          labels=c("0", 
                                   paste(thresholds[1:5], "-", thresholds[2:6]),
                                   paste(">", thresholds[6]))),
         threshold_num=as.numeric(threshold)) |>
  filter(threshold_num != 1) |>
  ggplot(aes(date, propSites, fill=threshold_num, group=threshold_num)) +
  geom_area(colour="grey30", linewidth=0.05, outline.type="both", alpha=0.5) +
  scale_x_date(date_breaks="14 days", 
               date_labels="%j", expand=expansion(mult=c(0.05, 0.05))) +
  scale_y_continuous("Proportion of active farms", limits=c(0,1)) +
  scale_fill_viridis_b(expression(paste("Ensemble mean daily copepodids" %.% "m"^"-3" %.% "h"^"-1")),
                       option="turbo", begin=0.05,
                       breaks=c(2.5, 3.5, 4.5, 5.5, 6.5),
                       labels=c("0.0001", "0.001", "0.01", "0.1", "1")) +
  theme(panel.grid.major.x=element_line(colour="grey", linewidth=0.6),
        panel.grid.minor.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=11),
        legend.position="bottom", 
        legend.key.width=unit(1.5, "cm"), 
        legend.key.height=unit(0.2, "cm"),
        strip.text=element_text(size=11)) 
ggsave("figs/pub_new/ens_influx_SpringNeap_14d.png", fig_influx, width=10, height=4, dpi=400)






influx_ens |>
  group_by(date) |>
  summarise(lice_CI95width=median(lice_q975-lice_q025)) |>
  ggplot(aes(date, lice_CI95width)) + geom_point(alpha=0.25)

influx_ens |>
  # filter(sepaSite=="AAC3") |>
  mutate(week=round_date(date, "month")) |>
  group_by(sepaSite, week) |>
  summarise(across(where(is.numeric), mean)) |>
  ggplot(aes(week, lice_mn^0.25, ymin=lice_q025^0.25, ymax=lice_q975^0.25)) + 
  geom_ribbon(alpha=0.5, colour=NA) + 
  geom_line() + 
  facet_wrap(~sepaSite)

influx_ens |>
  # filter(sepaSite=="AAC3") |>
  mutate(week=round_date(date, "month")) |>
  group_by(sepaSite, week) |>
  summarise(across(where(is.numeric), mean)) |>
  mutate(year=year(week),
         week_std=ymd(paste(2023, month(week), day(week), sep="-"))) |>
  ggplot(aes(week, lice_mn^0.25, group=week)) + 
  stat_slabinterval(normalize="groups") +
  theme_bw() 




# IP by simulation --------------------------------------------------------

thresholds <- c(0, 1e-4, 1e-3, 1e-2, 1e-1, 1)
influx_ens <- readRDS("out/sim_2021-2024/processed/influx_ens.rds")
c_daily <- readRDS(glue("out/sim_2021-2024/processed/connectivity_day.rds")) |>
  select(sepaSite, sim, date, influx_m3) |> 
  mutate(sim=paste0("sim_", sim)) |>
  filter(sepaSite %in% unique(influx_ens$sepaSite)) |>
  bind_rows(influx_ens |> select(sepaSite, date, lice_mn) |>
              rename(influx_m3=lice_mn) |>
              mutate(sim="predBlend")) |>
  left_join(sim_i |> mutate(lab=fct_relabel(lab, ~gsub("'Blend'", "Blend", .x))), by="sim") |>
  mutate(lice=influx_m3,
         tl=case_when(lice==0 ~ 1,
                      lice > thresholds[1] & lice < thresholds[2] ~ 2,
                      between(lice, thresholds[2], thresholds[3]) ~ 3,
                      between(lice, thresholds[3], thresholds[4]) ~ 4,
                      between(lice, thresholds[4], thresholds[5]) ~ 5,
                      between(lice, thresholds[5], thresholds[6]) ~ 6,
                      lice > thresholds[6] ~ 7)) |>
  filter(tl > 1) |>
  droplevels()
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")

sites <- sort(unique(c_daily$sepaSite))
site_groups <- split(sites, ceiling(seq_along(sites)/(3*5)))
mod_labs <- c(expression('2D.4'), expression('2D.3'), expression('2D.2'), expression('2D.1'),
              expression('3D.16'), expression('3D.15'), expression('3D.14'), expression('3D.13'),
              expression('3D.12'), expression('3D.11'), expression('3D.10'), expression('3D.9'),
              expression('3D.8'), expression('3D.7'), expression('3D.6'), expression('3D.5'),
              expression('3D.4'), expression('3D.3'), expression('3D.2'), expression('3D.1'),
              expression(Ens[Blend]))

for(i in seq_along(site_groups)) {
  p <- c_daily |>
    filter(sepaSite %in% site_groups[[i]]) |>
    ggplot(aes(date, lab, fill=tl)) + 
    geom_raster() + 
    scale_fill_viridis_b(expression(paste("Daily copepodids" %.% "m"^"-3" %.% "h"^"-1")),
                         option="turbo", begin=0.05,
                         breaks=c(2.5, 3.5, 4.5, 5.5, 6.5),
                         labels=c("0.0001", "0.001", "0.01", "0.1", "1")) +
    scale_y_discrete(limits=rev(levels(c_daily$lab)),
                     labels=mod_labs) +
    facet_wrap(~sepaSite, ncol=3, strip.position="right", axes="all_x", axis.labels="margins") +
    theme_classic() +
    theme(panel.grid.major.x=element_line(colour="grey90"),
          legend.position="bottom",
          legend.key.width=unit(1.5, "cm"), 
          legend.key.height=unit(0.2, "cm"),
          axis.title=element_blank(),
          axis.text.y=element_text(size=7)) 
  ggsave(glue("figs/pub_new/IP_by_site_{i}.png"), p, 
         height=270*ceiling(length(site_groups[[i]])/3)/5, width=190, units="mm")
}






# density sLL -------------------------------------------------------------

set.seed(1003)
mod <- "n20_D4"
out_ensBlend <- readRDS(glue("out/ensembles/ensBlend_{mod}_FULL_stanfit.rds"))
dat_ensBlend <- readRDS(glue("out/ensembles/ensBlend_{mod}_FULL_standata.rds"))

# from 03a*.R -- (weekly copepodid IP)^0.25 in each grid cell
f <- dirf("out/sim_2021-2024/processed/weekly", "Mature")
ens_ls <- vector("list", length(f))
ps_lims <- tibble(ens_mn=c(0,0),
                  ens_CL005=c(0,0),
                  ens_CL025=c(0,0),
                  ens_CL975=c(0,0),
                  ens_CL995=c(0,0),
                  ens_CI95width=c(0,0),
                  ens_CI99width=c(0,0),
                  sim_sd=c(0,0),
                  sim_CV=c(0,0),
                  ens_mn_orig=c(0,0),
                  ens_CL005_orig=c(0,0),
                  ens_CL025_orig=c(0,0),
                  ens_CL975_orig=c(0,0),
                  ens_CL995_orig=c(0,0),
                  ens_CI95width_orig=c(0,0),
                  ens_CI99width_orig=c(0,0))
p_dir <- "out/ensembles/p_meshCentroids/"

library(doFuture)
plan(multicore, workers=15)
for(i in 1:length(f)) {
  timestep <- ymd("2021-01-01") + dhours(as.numeric(str_sub(str_split_fixed(f[i], "_t_", 2)[,2], 1, -5)))
  ps_i <- readRDS(f[i]) |>
    select(i, all_of(dat_ensBlend$sim_names))
  ps_mx <- as.matrix(ps_i |> select(-i))
  
  # Calculate ensIP in parallel -- all on 4th rt scale
  ensIP <- foreach(j=1:nrow(ps_i), .combine=rbind, .inorder=TRUE, 
                   .options.future=list(globals=structure(TRUE, add=c("ps_mx", "p_dir", "ps_i")))) %dofuture% {
    ens_j <- ps_mx[j,,drop=F] %*% readRDS(glue("{p_dir}/i_{as.integer(ps_i$i[j])}.rds"))
    c(mean(ens_j), quantile(ens_j, probs=c(0.005, 0.025, 0.975, 0.995)))
                   }
  sim_mn <- rowMeans(ps_mx)
  sim_sd <- apply(ps_mx, 1, sd)
  ensIP <- cbind(ensIP, 
                 sim_sd,
                 sim_sd/sim_mn)
  gc()
  # Ensemble values on 4th rt scale, then back-transformed to original scale
  ens_ls[[i]] <- tibble(i=ps_i$i,
                        ens_mn=ensIP[,1],
                        ens_CL005=ensIP[,2],
                        ens_CL025=ensIP[,3],
                        ens_CL975=ensIP[,4],
                        ens_CL995=ensIP[,5],
                        ens_CI95width=ens_CL975-ens_CL025,
                        ens_CI99width=ens_CL995-ens_CL005,
                        sim_sd=ensIP[,6],
                        sim_CV=ensIP[,7]) |>
    mutate(ens_mn_orig=ens_mn^4,
           ens_CL005_orig=ens_CL005^4,
           ens_CL025_orig=ens_CL025^4,
           ens_CL975_orig=ens_CL975^4,
           ens_CL995_orig=ens_CL995^4,
           ens_CI95width_orig=ens_CL975_orig - ens_CL025_orig,
           ens_CI99width_orig=ens_CL995_orig - ens_CL005_orig)
  ps_lims$ens_mn <- range(c(ps_lims$ens_mn, range(ens_ls[[i]]$ens_mn)))
  ps_lims$ens_CL005 <- range(c(ps_lims$ens_CL005, range(ens_ls[[i]]$ens_CL005)))
  ps_lims$ens_CL025 <- range(c(ps_lims$ens_CL025, range(ens_ls[[i]]$ens_CL025)))
  ps_lims$ens_CL975 <- range(c(ps_lims$ens_CL975, range(ens_ls[[i]]$ens_CL975)))
  ps_lims$ens_CL995 <- range(c(ps_lims$ens_CL995, range(ens_ls[[i]]$ens_CL995)))
  ps_lims$ens_CI95width <- range(c(ps_lims$ens_CI95width, range(ens_ls[[i]]$ens_CI95width))) 
  ps_lims$ens_CI99width <- range(c(ps_lims$ens_CI99width, range(ens_ls[[i]]$ens_CI99width))) 
  ps_lims$sim_sd <- range(c(ps_lims$sim_sd, range(ens_ls[[i]]$sim_sd))) 
  ps_lims$sim_CV <- range(c(ps_lims$sim_CV, range(ens_ls[[i]]$sim_CV)))
  ps_lims$ens_mn_orig <- range(c(ps_lims$ens_mn_orig, range(ens_ls[[i]]$ens_mn_orig)))
  ps_lims$ens_CL005_orig <- range(c(ps_lims$ens_CL005_orig, range(ens_ls[[i]]$ens_CL005_orig)))
  ps_lims$ens_CL025_orig <- range(c(ps_lims$ens_CL025_orig, range(ens_ls[[i]]$ens_CL025_orig)))
  ps_lims$ens_CL975_orig <- range(c(ps_lims$ens_CL975_orig, range(ens_ls[[i]]$ens_CL975_orig)))
  ps_lims$ens_CL995_orig <- range(c(ps_lims$ens_CL995_orig, range(ens_ls[[i]]$ens_CL995_orig)))
  ps_lims$ens_CI95width_orig <- range(c(ps_lims$ens_CI95width_orig, range(ens_ls[[i]]$ens_CI95width_orig))) 
  ps_lims$ens_CI99width_orig <- range(c(ps_lims$ens_CI99width_orig, range(ens_ls[[i]]$ens_CI99width_orig))) 
  cat("Finished", as.character(timestep), "\n")
  gc()
}
saveRDS(ps_lims, "out/sim_2021-2024/processed/ps_lims.rds")

timesteps <- ymd("2019-04-01") + dhours(as.numeric(str_sub(str_split_fixed(f, "_t_", 2)[,2], 1, -5)))
ens_df <- map2_dfr(ens_ls, timesteps, ~.x |> mutate(date=.y))
saveRDS(ens_df, "out/sim_2021-2024/processed/ens_weekly.rds")

ens_avg <- ens_df |>
  filter(date >= "2021-05-01") |>
  group_by(i) |>
  summarise(across(where(is.numeric), .fn=list(mn=mean, md=median))) |>
  ungroup()
saveRDS(ens_avg, "out/sim_2021-2024/processed/ens_avg_n20_D4.rds")




# maps --------------------------------------------------------------------

# Left side: Ensemble mean(copepodid density)
# Right side: Ensemble mean(weekly CI width)
# ens_df <- readRDS("out/sim_2021-2024/processed/ens_weekly.rds")
ens_avg <- readRDS("out/sim_2021-2024/processed/ens_avg_n20_D4.rds")

# WeStCOMS mesh
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_sf <- st_read("data/WeStCOMS2_mesh.gpkg") |> select(i, geom)
linnhe_mesh <- mesh_sf |> 
  # st_crop(c(xmin=150000, xmax=220000, ymin=710000, ymax=785000))
  st_crop(c(xmin=130000, xmax=220000, ymin=710000, ymax=785000))
skye_mesh <- mesh_sf |> 
  st_crop(c(xmin=100000, xmax=198000, ymin=780000, ymax=920000))

westcoms_panel <- ggplot() +
  geom_sf(data=mesh_fp, fill="grey", colour="grey30", linewidth=0.2) +
  guides(fill=guide_colourbar(title.position="top", direction="horizontal")) +
  scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), ".0\u00B0W")) +
  scale_y_continuous(breaks=c(54, 56, 58), labels=paste0(c(54, 56, 58), ".0\u00B0N")) +
  theme(legend.position=c(0.285, 0.1),
        legend.background=element_blank(),
        legend.key.height=unit(0.1, "cm"),
        legend.key.width=unit(0.43, "cm"),
        legend.title=element_blank(),
        legend.text=element_text(size=7),
        axis.title=element_blank())
linnhe_panel <- ggplot() +
  geom_sf(data=mesh_fp, fill="grey", colour="grey30", linewidth=0.2) +
  guides(fill=guide_colourbar(title.position="top", direction="horizontal")) +
  scale_x_continuous(limits=c(143000, 216000), breaks=c(-5.8, -5.4, -5)) +
  scale_y_continuous(limits=c(720000, 778000), breaks=c(56.4, 56.7)) +
  theme(legend.position="none",
        axis.title=element_blank())
skye_panel <- ggplot() +
  geom_sf(data=mesh_fp, fill="grey", colour="grey30", linewidth=0.2) +
  guides(fill=guide_colourbar(title.position="top", direction="horizontal")) +
  scale_x_continuous(limits=c(110000, 194000), breaks=c(-6.5, -6, -5.5)) +
  scale_y_continuous(limits=c(786000, 899000), breaks=c(57, 57.5)) +
  theme(legend.position="none",
        axis.title=element_blank())

mn_lims <- c(0, 0.002^0.25)
mn_breaks <- c(0, 0.0001, 0.001)
mn_labs <- c("0", "1e-4", "1e-3")

ci_lims <- c(0, 0.15)
ci_breaks <- c(0, 1e-6, 1e-4)
ci_labs <- c("0", "1e-6", "1e-4")


# mn_lims <- c(0, 0.1)
# mn_breaks <- c(0, 0.01, 0.05, 0.1)

ens_avg <- ens_avg |>
  mutate(ens_mn_mn=pmin(ens_mn_mn/5, mn_lims[2]),
         ens_CI95width_mn=pmin(ens_CI95width_mn, ci_lims[2]))

ens_map <- vector("list", 6)
# Mean densities
ens_map[[1]] <- westcoms_panel + 
  geom_sf(data=ens_avg |> right_join(mesh_sf, y=_),
          aes(fill=ens_mn_mn), colour=NA) + 
  scale_fill_viridis_c("",
                       option="turbo", limits=mn_lims,
                       breaks=mn_breaks^0.25, labels=mn_labs) +
  annotate("text", x=79000, y=547000, label="Ensemble mean", size=3) +
  annotate("text", x=79000, y=525000, 
           label=expression("cop." %.% "m"^"-3" %.% "h"^"-1"), parse=T, size=3)
ens_map[[2]] <- linnhe_panel + 
  geom_sf(data=ens_avg |> right_join(linnhe_mesh, y=_),
          aes(fill=ens_mn_mn), colour=NA) + 
  scale_fill_viridis_c(option="turbo", limits=mn_lims,
                       breaks=mn_breaks^0.25, labels=mn_labs) 
ens_map[[3]] <- skye_panel + 
  geom_sf(data=ens_avg |> right_join(skye_mesh, y=_),
          aes(fill=ens_mn_mn), colour=NA) + 
  scale_fill_viridis_c(option="turbo", limits=mn_lims,
                       breaks=mn_breaks^0.25, labels=mn_labs) 
ens_map[[4]] <- westcoms_panel + 
  geom_sf(data=ens_avg |> right_join(mesh_sf, y=_),
          aes(fill=ens_CI95width_mn), colour=NA) + 
  scale_fill_scico(palette="imola", limits=ci_lims,
                   breaks=ci_breaks^0.25, labels=ci_labs) +
  annotate("text", x=79000, y=547000, label="95% CI width", size=3) +
  annotate("text", x=79000, y=525000, 
           label=expression("cop." %.% "m"^"-3" %.% "h"^"-1"), parse=T, size=3)
ens_map[[5]] <- linnhe_panel + 
  geom_sf(data=ens_avg |> right_join(linnhe_mesh, y=_),
          aes(fill=ens_CI95width_mn), colour=NA) + 
  scale_fill_scico(palette="imola", limits=ci_lims,
                   breaks=ci_breaks^0.25, labels=ci_labs)
ens_map[[6]] <- skye_panel + 
  geom_sf(data=ens_avg |> right_join(skye_mesh, y=_),
          aes(fill=ens_CI95width_mn), colour=NA) + 
  scale_fill_scico(palette="imola", limits=ci_lims,
                   breaks=ci_breaks^0.25, labels=ci_labs)

plot_grid(plotlist=ens_map, ncol=2, nrow=3, labels="auto", byrow=FALSE,
          rel_heights=c(2.1, 0.8, 1.23), rel_widths=c(1, 1)) |>
  ggsave("figs/pub_new/ens_map_n20_D4.png", plot=_, width=4.75, height=9.1, dpi=600)

ggsave("figs/talk/ens_map_WeStCOMS.png", ens_map[[1]], width=3.25, height=7, dpi=300)





# fig overview inset ------------------------------------------------------

ens_avg <- readRDS("out/sim_2019-2023/processed/ens_avg_n20_D4_all.rds")

# WeStCOMS mesh
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
examp_mesh <- st_read("data/WeStCOMS2_mesh.gpkg") |> select(i, geom) |> 
  st_crop(c(xmin=90000, xmax=200000, ymin=820000, ymax=920000))
ggplot() + 
  geom_sf(data=mesh_fp, fill="grey", colour="grey30", linewidth=0.2) + 
  scale_x_continuous(limits=c(100000, 194000)) + 
  scale_y_continuous(limits=c(850000, 899000)) + 
  theme(axis.title=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        legend.position="none") + 
  geom_sf(data=ens_avg |> right_join(examp_mesh, y=_), aes(fill=ens_mn), colour=NA) + 
  scale_fill_viridis_c(option="turbo")
ggsave("figs/pub_new/fig_overview_example_map.png", width=5.5, height=3)

mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
examp_mesh <- st_read("data/WeStCOMS2_mesh.gpkg") |> select(i, geom) |> 
  st_crop(c(xmin=12000, xmax=247000, ymin=470000, ymax=1000000))
ggplot() + 
  geom_sf(data=mesh_fp, fill="grey", colour="grey30", linewidth=0.2) + 
  scale_x_continuous(limits=c(50000, 250000)) + 
  scale_y_continuous(limits=c(650000, 900000)) + 
  theme(axis.title=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        legend.position="none",
        panel.border=element_blank()) + 
  geom_sf(data=ens_avg |> right_join(examp_mesh, y=_), aes(fill=ens_mn), colour=NA) + 
  scale_fill_viridis_c(option="turbo")
ggsave("figs/pub_new/fig_overview_example_map2.png", width=5.5, height=7)




# rmse stormcloud ---------------------------------------------------------

# ensCV_df <- read_csv("out/ensemble_CV.csv")

# sim_07 would be selected as 'optimal'
farm_r.df <- metrics_by_farm |>
  filter(N >= 10) |>
  filter(!grepl("null", sim)) |>
  # filter(sim %in% c("predFcst", "predBlend", "sim_07", "sim_avgAll")) |>
  left_join(sim_i) |>
  droplevels() |>
  select(sepaSite, sim, N, prop_g05, rmse, rho, ROC_AUC, lab, lab_short) |>
  pivot_longer(any_of(c("rmse", "rho", "ROC_AUC")), names_to="metric") |>
  mutate(type="By farm") |>
  mutate(metric=factor(metric, levels=c("ROC_AUC", "rho", "rmse"),
                     labels=c("'AUC'['ROC']", "rho", "RMSE"))) |>
  arrange(lab) |>
  mutate(type=factor(type, 
                     levels=c("global", "By farm", "By week"),
                     labels=c("Global", "'Within farm'", "'Within week'")),
         lab=if_else(sim %in% sim_i$sim[c(1:6,8:20)], "Other", lab),
         lab=if_else(sim=="sim_07", "Opt['Param']", lab),
         lab_short=if_else(sim %in% sim_i$sim[c(1:6,8:20)], "Other", lab_short),
         lab_short=if_else(sim=="sim_07", "Opt['Param']", lab_short)) |>
  filter(!is.na(value)) |>
  mutate(lab=factor(lab, levels=names(modType3_cols)),
         lab_short=factor(lab_short, levels=names(modType3_cols)))
  # mutate(lab=lvls_revalue(lab, c("Ens['Fcst']", "Ens['Blend']", "Ens['Avg']", "Opt['Param']")),
  #        lab_short=lvls_revalue(lab_short, c("Ens['Fcst']", "Ens['Blend']", "Ens['Avg']", "Opt['Param']")))
week_r.df <- metrics_by_week |>
  filter(N >= 10) |>
  filter(!grepl("null", sim)) |>
  # filter(sim %in% c("predFcst", "predBlend", "sim_07", "sim_avgAll")) |>
  left_join(sim_i) |>
  droplevels() |>
  select(date, sim, N, prop_g05, rmse, rho, ROC_AUC, lab, lab_short) |>
  pivot_longer(any_of(c("rmse", "rho", "ROC_AUC")), names_to="metric") |>
  mutate(type="By week") |>
  mutate(metric=factor(metric, levels=c("ROC_AUC", "rho", "rmse"),
                       labels=c("'AUC'['ROC']", "rho", "RMSE"))) |>
  arrange(lab) |>
  mutate(type=factor(type, 
                     levels=c("global", "By farm", "By week"),
                     labels=c("Global", "'Within farm'", "'Within week'")),
         lab=if_else(sim %in% sim_i$sim[c(1:6,8:20)], "Other", lab),
         lab=if_else(sim=="sim_07", "Opt['Param']", lab),
         lab_short=if_else(sim %in% sim_i$sim[c(1:6,8:20)], "Other", lab_short),
         lab_short=if_else(sim=="sim_07", "Opt['Param']", lab_short)) |>
  filter(!is.na(value)) |>
  mutate(lab=factor(lab, levels=names(modType3_cols)),
         lab_short=factor(lab_short, levels=names(modType3_cols)))
  # mutate(lab=lvls_revalue(lab, c("Ens['Fcst']", "Ens['Blend']", "Ens['Avg']", "Opt['Param']")),
  #        lab_short=lvls_revalue(lab_short, c("Ens['Fcst']", "Ens['Blend']", "Ens['Avg']", "Opt['Param']")))
mn_ci <- bind_rows(farm_r.df, week_r.df) |>
  droplevels() |>
  group_by(sim, lab, lab_short, metric, type) |>
  summarise(mn=mean(value, na.rm=T),
            md=median(value, na.rm=T),
            N=sum(!is.na(value)),
            se=sd(value, na.rm=T)/sqrt(N),
            ci_lo=mn - qt(0.975, N-1)*se,
            ci_hi=mn + qt(0.975, N-1)*se)

all_metrics_mns <- all_metrics_df |>
  filter(grepl("(Ens|3D|2D)", lab)) |>
  droplevels() |>
  bind_rows(all_metrics_df |>
              filter(sim=="sim_07") |>
              mutate(lab_short="Opt['Param']",
                     lab="Opt['Param']")) |>
  mutate(type=paste0("'", type, "'")) |>
  mutate(lab_short=factor(lab_short, levels=levels(farm_r.df$lab_short))) |>
  mutate(type=factor(type, 
                     levels=c("global", "'By farm'", "'By week'"),
                     labels=c("Global", "'Within farm'", "'Within week'")))

# Maps among weeks are much more stable
# More variability among farms in predicting time series
pA_df <- bind_rows(farm_r.df, week_r.df) |>
  filter(metric=="RMSE") |>
  droplevels()
pA <- pA_df |>
  ggplot(aes(value, lab, fill=lab_short, colour=lab_short)) + 
  geom_dots(data=pA_df |> filter(lab != "Other"), side="bottom", scale=0.5) + 
  stat_slab(normalize="xy", scale=0.5, colour=NA, fill_type="gradient",
            aes(slab_alpha=after_stat(-pmax(abs(1-2*cdf), 0.25)))) +
  stat_pointinterval(.width=c(0.5, 0.8), colour="black", fatten_point=1.25, point_interval="mean_qi") +
  geom_rug(data=all_metrics_mns |> filter(metric=="RMSE") |> filter(!lab_short %in% c("2D", "3D")),
           aes(x=value), sides="b", length=unit(0.05, "npc"), linewidth=0.5, alpha=0.75) +
  geom_rug(data=all_metrics_mns |> filter(metric=="RMSE") |> filter(lab_short %in% c("2D", "3D")),
           aes(x=value), sides="b", length=unit(0.035, "npc"), linewidth=0.2) +
  scale_slab_alpha_continuous(range=c(0.01, 0.75), guide="none") +
  scale_fill_manual(values=modType3_cols) +
  scale_colour_manual(values=modType3_cols) +
  scale_y_discrete(breaks=levels(mn_ci$lab), labels=parse(text=levels(mn_ci$lab)), limits=levels(mn_ci$lab)) +
  labs(x="Cross validation score") +
  facet_grid(metric~type, labeller="label_parsed") +
  theme(panel.grid.major=element_line(colour="grey90", linewidth=0.2),
        strip.text=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=9),
        legend.position="none")
pB_df <- bind_rows(farm_r.df, week_r.df) |>
  filter(metric=="rho") |>
  droplevels()
pB <- pB_df |>
  ggplot(aes(value, lab, fill=lab_short, colour=lab_short)) + 
  geom_dots(data=pB_df |> filter(lab != "Other"), side="bottom", scale=0.5) + 
  stat_slab(normalize="xy", scale=0.5, colour=NA, fill_type="gradient",
            aes(slab_alpha=after_stat(-pmax(abs(1-2*cdf), 0.25)))) +
  stat_pointinterval(.width=c(0.5, 0.8), colour="black", fatten_point=1.25, point_interval="mean_qi") +
  geom_rug(data=all_metrics_mns |> filter(metric=="rho") |> filter(!lab_short %in% c("2D", "3D")),
           aes(x=value), sides="b", length=unit(0.05, "npc"), linewidth=0.5, alpha=0.75) +
  geom_rug(data=all_metrics_mns |> filter(metric=="rho") |> filter(lab_short %in% c("2D", "3D")),
           aes(x=value), sides="b", length=unit(0.035, "npc"), linewidth=0.2) +
  scale_slab_alpha_continuous(range=c(0.01, 0.75), guide="none") +
  
  scale_fill_manual(values=modType3_cols) +
  scale_colour_manual(values=modType3_cols) +
  scale_y_discrete(breaks=levels(mn_ci$lab), labels=parse(text=levels(mn_ci$lab)), limits=levels(mn_ci$lab)) +
  labs(x="Cross validation score") +
  facet_grid(metric~type, labeller="label_parsed") +
  theme(panel.grid.major=element_line(colour="grey90", linewidth=0.2),
        strip.text=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=9),
        legend.position="none")
pC_df <- bind_rows(farm_r.df, week_r.df) |>
  filter(metric=="'AUC'['ROC']") |>
  filter(prop_g05 > 0 & prop_g05 < 1) |>
  droplevels()
pC <- pC_df |>
  ggplot(aes(value, lab, fill=lab_short, colour=lab_short)) + 
  geom_dots(data=pC_df |> filter(lab != "Other"), side="bottom", scale=0.5) + 
  stat_slab(normalize="xy", scale=0.5, colour=NA, fill_type="gradient",
            aes(slab_alpha=after_stat(-pmax(abs(1-2*cdf), 0.25)))) +
  stat_pointinterval(.width=c(0.5, 0.8), colour="black", fatten_point=1.25, point_interval="mean_qi") +
  geom_rug(data=all_metrics_mns |> filter(metric=="'AUC'['ROC']") |> filter(!lab_short %in% c("2D", "3D")),
           aes(x=value), sides="b", length=unit(0.05, "npc"), linewidth=0.5, alpha=0.75) +
  geom_rug(data=all_metrics_mns |> filter(metric=="'AUC'['ROC']") |> filter(lab_short %in% c("2D", "3D")),
           aes(x=value), sides="b", length=unit(0.035, "npc"), linewidth=0.2) +
  scale_slab_alpha_continuous(range=c(0.01, 0.75), guide="none") +
  
  scale_fill_manual(values=modType3_cols) +
  scale_colour_manual(values=modType3_cols) +
  scale_y_discrete(breaks=levels(mn_ci$lab), labels=parse(text=levels(mn_ci$lab)), limits=levels(mn_ci$lab)) +
  labs(x="Cross validation score") +
  facet_grid(metric~type, labeller="label_parsed") +
  theme(panel.grid.major=element_line(colour="grey90", linewidth=0.2),
        strip.text=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=9),
        legend.position="none")
p <- cowplot::plot_grid(pA, pB, pC, ncol=1, align="hv", axis="tblr", labels="auto")
ggsave("figs/pub_new/metric_stormclouds.png", p, width=9, height=9, dpi=300)





# rank stormcloud ---------------------------------------------------------

# ensCV_df <- read_csv("out/ensemble_CV.csv")

p <- metric_ranks |>
  filter(N >= 10) |>
  mutate(lab=factor(lab, 
                    levels=c("Null['farm']", "Null['time']", 
                             paste0("'2D.", c(1, 2, 3, 4), "'"),
                             paste0("'3D.", c(13, 4, 2, 6, 11, 8, 14, 5,
                                              15, 10, 16, 12, 1, 9, 3, 7), "'"),
                             paste0("Ens['", c("Avg", "Blend", "Fcst"), "']")) |>
                      rev())) |>
  arrange(lab) |>
  mutate(type=factor(type, 
                     levels=c("global", "byFarm", "byWeek"),
                     labels=c("Global", "Mean within farm", "Mean within week"))) |>
  ggplot(aes(rank-0.5, lab, fill=lab_short, colour=lab_short)) + 
  stat_histinterval(normalize="xy", scale=0.5, alpha=0.5, 
                    outline_bars=T, slab_colour="black", slab_linewidth=0.2,
                    breaks=breaks_fixed(width=1)) +
  stat_pointinterval(.width=c(0.5, 0.8), colour="black", fatten_point=1.2) +
  # scale_fill_scico_d(palette="glasgow", guide="none", end=0.8) +
  # scale_colour_scico_d(palette="glasgow", guide="none", end=0.8) +
  scale_fill_manual(values=modType3_cols) +
  scale_colour_manual(values=modType3_cols) +
  labs(x="Rank") +
  scale_y_discrete(labels=label_parsed) +
  facet_grid(.~type) +
  theme(legend.position="none",
        panel.grid.major=element_line(colour="grey90", linewidth=0.2),
        strip.text=element_text(size=11),
        axis.title.x=element_text(size=9),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=9))
ggsave("figs/pub_new/RMSE_stormclouds_ranks_all.png", p, width=10, height=10, dpi=300)

metric_ranks |>
  filter(N >= 10) |>
  arrange(lab) |>
  mutate(type=factor(type, 
                     levels=c("global", "byFarm", "byWeek"),
                     labels=c("Global", "'By farm'", "'By week'"))) |>
  group_by(lab, type) |>
  summarise(mnRank=mean(rank),
            sdRank=sd(rank),
            q25=quantile(rank, probs=0.25),
            q75=quantile(rank, probs=0.75),
            mdRank=median(rank)) |>
  arrange(type, mnRank) |>
  print(n=28)

metric_ranks |>
  filter(N >= 10) |>
  arrange(lab) |>
  mutate(type=factor(type, 
                     levels=c("global", "byFarm", "byWeek"),
                     labels=c("Global", "'By farm'", "'By week'"))) |>
  group_by(lab) |>
  summarise(mnRank=mean(rank),
            sdRank=sd(rank),
            q25=quantile(rank, probs=0.25),
            q75=quantile(rank, probs=0.75),
            mdRank=median(rank)) |>
  arrange(mnRank) |>
  print(n=28)
  


# maps of rmse + IDW ---------------------------------------------------------

library(terra)
#ensCV_df <- read_csv("out/ensemble_CV.csv")
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_rast <- st_read("data/WeStCOMS2_meshFootprint.gpkg") |>
  rast(resolution=500)
locs <- st_read("figs/place_names.gpkg")

rmse_info <- tibble(breaks=seq(0.1, 0.75, by=0.05)) |>
  mutate(break_labs=as.character(round(breaks, 1)),
         break_labs=if_else(row_number() %% 2 == 0, "", break_labs),
         letter=letters[row_number()],
         mdpt=(breaks + (lead(breaks)-breaks)/2))
farm_rmse.df <- ensCV_df |>
  group_by(sepaSite) |>
  summarise(across(starts_with("IP"),
                   ~yardstick::rmse_vec(.x, truth=licePerFish_rtrt))) |>
  inner_join(site_i)

p_ls <- vector("list", 4)
mods <- c("IP_predFf_RMSE", "IP_predBlend", "IP_sim_avgAll", "IP_sim_07")
col_labs <- c(expression(Ens['Fcst']), expression(Ens['Blend']), 
              expression(Ens['Avg']), expression(Opt['Param']))
for(i in seq_along(mods)) {
  col_lab <- col_labs[i]
  map_interp <- interpIDW(mesh_rast,
                          farm_rmse.df |>
                            rename_with(~"predColumn", .cols=matches(mods[i])) |>
                            select(easting, northing, predColumn) |>
                            drop_na() |>
                            as.matrix(),
                          radius=1000e3) |>
    mask(mesh_fp)
  farm_rmse.df_i <- farm_rmse.df |>
    rename_with(~"predColumn", .cols=matches(mods[i])) |>
    filter(!is.na(predColumn)) |>
    select(sepaSite, predColumn, easting, northing) |>
    mutate(letter=cut(predColumn,
                      breaks=rmse_info$breaks,
                      labels=letters[1:(length(rmse_info$breaks)-1)]))
  farm_rmse_count <- farm_rmse.df_i |>
    count(letter) |>
    mutate(scaled=n/max(n)) |>
    full_join(rmse_info |> select(letter, mdpt) |> drop_na()) |>
    mutate(n=replace_na(n, 0),
           scaled=replace_na(scaled, 0)) |>
    arrange(letter)
  low_polygon <- tibble(x=c(81000, 96000, 96000, 81000)+4000,
                        y=rep(c(0, 54800/nrow(farm_rmse_count)), each=2) + 652000)
  x_rng <- diff(range(low_polygon$x))
  y_rng <- diff(range(low_polygon$y))
  farm_rmse_bar.df <- map_dfr(1:nrow(farm_rmse_count),
                           ~low_polygon |>
                             mutate(mdpt=farm_rmse_count$mdpt[.x],
                                    x=if_else(x==max(x),
                                              x,
                                              max(x)-x_rng*farm_rmse_count$scaled[.x]),
                                    y=y + (y_rng*(.x-1)))
  )
  farm_rmse_count_labs <- farm_rmse_bar.df |>
    group_by(mdpt) |>
    summarise(x=min(x), y=mean(y)) |>
    ungroup() |>
    left_join(farm_rmse_count) |>
    mutate(prop=paste0(round(n/sum(n)*100), "%"))

  p_ls[[i]] <- as_tibble(map_interp) |>
    bind_cols(crds(map_interp)) |>
    ggplot() +
    geom_sf(data=mesh_fp, fill="grey90", colour="grey", size=0.1) +
    geom_raster(aes(x, y, fill=lyr.1)) +
    geom_point(data=farm_rmse.df_i, aes(easting, northing, fill=predColumn),
               shape=21, size=1, stroke=0.25, colour="grey10") +
    geom_sf_text(data=locs, aes(label=Label), colour="magenta", fontface="bold", 
                 size=3.5, 
                 nudge_x=c(0, 5, -2, -8, 7, -6, 5, 2) * 1e3, 
                 nudge_y=c(10, -1, 8, 0, 3, 0, 3, -4) * 1e3) +
    geom_polygon(data=farm_rmse_bar.df, aes(x, y, fill=mdpt, group=mdpt),
                 colour="grey10", linewidth=0.15) +
    geom_text(data=farm_rmse_count_labs, aes(x, y, label=prop),
              size=2, hjust=1, vjust=0.5, nudge_x=-1000) +
    colorspace::scale_fill_binned_diverging(
      name=col_lab, palette="Blue-Red 3", rev=F,
      limits=c(0.1, 0.75), mid=0.435, breaks=rmse_info$breaks, labels=rmse_info$break_labs) +
    scale_y_continuous(limits=c(630000, 955000), oob=scales::oob_keep,
                       breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
    scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
                       limits=c(75000, 235000), oob=scales::oob_keep) +
    theme(legend.position="inside",
          legend.position.inside=c(c(0.194,0.183,0.196,0.18)[i], 0.204),
          legend.background=element_blank(),
          legend.key.height=unit(0.43, "cm"),
          legend.key.width=unit(0.0, "cm"),
          legend.text=element_text(size=7),
          legend.title=element_text(size=10, vjust=1, hjust=1),
          legend.ticks=element_line(colour="grey10", linewidth=0.25),
          legend.ticks.length=unit(0.04, "cm"),
          axis.title=element_blank())
  if(i > 1) p_ls[[i]] <- p_ls[[i]] + theme(axis.text.y=element_blank())
}

p_farms <- plot_grid(plotlist=p_ls, labels="auto", align="hv", axis="tblr", nrow=1)
# p_farms |>
#   ggsave("figs/pub_new/ens_farm-rmse+IDW_map.png", plot=_ , width=13, height=6, dpi=300)

p_weeklyA <- metric_date_df |>
  filter(metric=="RMSE") |>
  filter(sim %in% c("predFcst", "predBlend", "sim_avgAll", "sim_07")) |>
  mutate(lab=if_else(sim=="sim_07", "Opt['Param']", lab),
         lab=factor(lab, levels=names(modType3_cols)[1:4])) |>
  ggplot(aes(date, value, colour=lab)) + 
  geom_point(shape=1) + 
  geom_line(stat="smooth", method="gam", formula=y~s(x, k=15), se=F) +
  scale_x_date(date_breaks="1 year", 
               date_labels="%Y", expand=expansion(mult=c(0.05, 0.1))) +
  ylab("Weekly RMSE") +
  scale_colour_manual(values=modType3_cols[1:4],
                      labels=c(bquote(Ens['Fcst']), 
                               bquote(Ens['Blend']),
                               bquote(Ens['Avg']),
                               bquote(Opt['Param']))) +
  theme_bw() + 
  guides(colour=guide_legend(override.aes=list(size=1), title=NULL),
         linetype=guide_legend(title=NULL),
         shape=guide_legend(title=NULL)) +
  theme(legend.position="inside",
        legend.position.inside=c(0.92, 0.85),
        legend.background=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=9),
        panel.grid.major.y=element_line(colour="grey85", linewidth=0.4),
        panel.grid.minor.y=element_line(colour="grey90", linewidth=0.2),
        axis.text=element_text(size=8))


p_weeklyB <- metric_date_df |>
  filter(metric=="RMSE") |>
  filter(sim %in% c("predFcst", "predBlend", "sim_avgAll", "sim_07")) |>
  mutate(lab=if_else(sim=="sim_07", "Opt['Param']", lab),
         lab=factor(lab, levels=names(modType3_cols)[1:4])) |>
  mutate(date_std=ymd("2020-12-31")+yday(date),
         year=year(date)) |>
  ggplot(aes(date_std, value, colour=lab)) + 
  geom_point(alpha=0.8, size=0.9, shape=1) + 
  geom_line(stat="smooth", method="gam", formula=y~s(x, k=10, bs="cc"), se=F) +
  scale_x_date(date_breaks="1 month", date_labels="%b") +
  ylab("Weekly RMSE") + 
  scale_colour_manual(values=modType3_cols[1:4],
                      labels=c(bquote(Ens['Fcst']), 
                               bquote(Ens['Blend']),
                               bquote(Ens['Avg']),
                               bquote(Opt['Param']))) +
  theme_bw() + 
  guides(colour=guide_legend(override.aes=list(size=1), title=NULL),
         linetype=guide_legend(title=NULL),
         shape=guide_legend(title=NULL)) +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=9),
        panel.grid.major.y=element_line(colour="grey85", linewidth=0.2),
        panel.grid.minor.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        axis.text=element_text(size=8))

# p_weekly |>
#   ggsave("figs/pub_new/validation_metrics_byWeek.png", plot=_, width=6.5, height=4)

plot_grid(p_farms, 
          plot_grid(p_weeklyA, p_weeklyB, nrow=1, rel_widths=c(1, 0.5), align="h", axis="tb", labels=c("e", "f")), 
          nrow=2, rel_heights=c(1, 0.68)) |>
  ggsave("figs/pub_new/ens_RMSE_farm_week.png", plot=_, width=13, height=10)

# ggsave("figs/talk/rmse+IDW_ensFcst.png", p_ls[[1]], width=3, height=5.5)
# ggsave("figs/talk/rmse+IDW_ensBlend.png", p_ls[[2]], width=3, height=5.5)
# ggsave("figs/talk/rmse+IDW_ens3D.png", p_ls[[3]], width=3, height=5.5)
# ggsave("figs/talk/rmse+IDW_ens2D.png", p_ls[[4]], width=3, height=5.5)



# maps of rho + IDW ---------------------------------------------------------

library(terra)
#ensCV_df <- read_csv("out/ensemble_CV.csv")
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_rast <- st_read("data/WeStCOMS2_meshFootprint.gpkg") |>
  rast(resolution=500)

rho_info <- tibble(breaks=seq(-1, 1, by=0.2)) |>
  mutate(break_labs=as.character(round(breaks, 1)),
         break_labs=if_else(row_number() %% 2 == 0, "", break_labs),
         letter=letters[row_number()],
         mdpt=(breaks + (lead(breaks)-breaks)/2))
farm_rho.df <- ensCV_df |>
  group_by(sepaSite) |>
  summarise(across(starts_with("IP"),
                   ~cor(.x, licePerFish_rtrt, method="spearman", use="pairwise"))) |>
  inner_join(site_i)

p_ls <- vector("list", 4)
mods <- c("IP_predFf_RMSE", "IP_predBlend", "IP_sim_avgAll", "IP_sim_07")
col_labs <- c(expression(Ens['Fcst']), expression(Ens['Blend']), 
              expression(Ens['Avg']), expression(Opt['Param']))
for(i in seq_along(mods)) {
  col_lab <- col_labs[i]
  map_interp <- interpIDW(mesh_rast,
                          farm_rho.df |>
                            rename_with(~"predColumn", .cols=matches(mods[i])) |>
                            select(easting, northing, predColumn) |>
                            drop_na() |>
                            as.matrix(),
                          radius=1000e3) |>
    mask(mesh_fp)
  farm_rho.df_i <- farm_rho.df |>
    rename_with(~"predColumn", .cols=matches(mods[i])) |>
    filter(!is.na(predColumn)) |>
    select(sepaSite, predColumn, easting, northing) |>
    mutate(letter=cut(predColumn,
                      breaks=rho_info$breaks,
                      labels=letters[1:(length(rho_info$breaks)-1)]))
  farm_rho_count <- farm_rho.df_i |>
    count(letter) |>
    mutate(scaled=n/max(n)) |>
    full_join(rho_info |> select(letter, mdpt) |> drop_na()) |>
    mutate(n=replace_na(n, 0),
           scaled=replace_na(scaled, 0)) |>
    arrange(letter)
  low_polygon <- tibble(x=c(81000, 96000, 96000, 81000)+4000,
                        y=rep(c(0, 54800/nrow(farm_rho_count)), each=2) + 652000)
  x_rng <- diff(range(low_polygon$x))
  y_rng <- diff(range(low_polygon$y))
  farm_rho_bar.df <- map_dfr(1:nrow(farm_rho_count),
                              ~low_polygon |>
                                mutate(mdpt=farm_rho_count$mdpt[.x],
                                       x=if_else(x==max(x),
                                                 x,
                                                 max(x)-x_rng*farm_rho_count$scaled[.x]),
                                       y=y + (y_rng*(.x-1)))
  )
  farm_rho_count_labs <- farm_rho_bar.df |>
    group_by(mdpt) |>
    summarise(x=min(x), y=mean(y)) |>
    ungroup() |>
    left_join(farm_rho_count) |>
    mutate(prop=paste0(round(n/sum(n)*100), "%"))
  
  p_ls[[i]] <- as_tibble(map_interp) |>
    bind_cols(crds(map_interp)) |>
    ggplot() +
    geom_sf(data=mesh_fp, fill="grey90", colour="grey", size=0.1) +
    geom_raster(aes(x, y, fill=lyr.1)) +
    geom_point(data=farm_rho.df_i, aes(easting, northing, fill=predColumn),
               shape=21, size=1, stroke=0.25, colour="grey10") +
    geom_sf_text(data=locs, aes(label=Label), colour="magenta", fontface="bold", 
                 size=3.5, 
                 nudge_x=c(0, 5, -2, -8, 7, -6, 5, 2) * 1e3, 
                 nudge_y=c(10, -1, 8, 0, 3, 0, 3, -4) * 1e3) +
    geom_polygon(data=farm_rho_bar.df, aes(x, y, fill=mdpt, group=mdpt),
                 colour="grey10", linewidth=0.15) +
    geom_text(data=farm_rho_count_labs, aes(x, y, label=prop),
              size=2, hjust=1, vjust=0.5, nudge_x=-1000) +
    colorspace::scale_fill_binned_diverging(
      name=col_lab, palette="Blue-Red 3", rev=T,
      limits=c(-1, 1), mid=0, breaks=rho_info$breaks, labels=rho_info$break_labs) +
    scale_y_continuous(limits=c(630000, 955000), oob=scales::oob_keep,
                       breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
    scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
                       limits=c(75000, 235000), oob=scales::oob_keep) +
    theme(legend.position="inside",
          legend.position.inside=c(c(0.212,0.202,0.218,0.2)[i], 0.204),
          legend.background=element_blank(),
          legend.key.height=unit(0.43, "cm"),
          legend.key.width=unit(0.0, "cm"),
          legend.text=element_text(size=7),
          legend.title=element_text(size=10, vjust=1, hjust=1),
          legend.ticks=element_line(colour="grey10", linewidth=0.25),
          legend.ticks.length=unit(0.04, "cm"),
          axis.title=element_blank())
  if(i > 1) p_ls[[i]] <- p_ls[[i]] + theme(axis.text.y=element_blank())
}

p_farms <- plot_grid(plotlist=p_ls, labels="auto", align="hv", axis="tblr", nrow=1)
# p_farms |>
#   ggsave("figs/pub_new/ens_farm-rho+IDW_map.png", plot=_ , width=13, height=6, dpi=300)

p_weeklyA <- metric_date_df |>
  filter(metric=="rho") |>
  filter(sim %in% c("predFcst", "predBlend", "sim_avgAll", "sim_07")) |>
  mutate(lab=if_else(sim=="sim_07", "Opt['Param']", lab),
         lab=factor(lab, levels=names(modType3_cols)[1:4])) |>
  ggplot(aes(date, value, colour=lab)) + 
  geom_point(shape=1) + 
  geom_line(stat="smooth", method="gam", formula=y~s(x, k=15), se=F) +
  scale_x_date(date_breaks="1 year", #date_minor_breaks="3 months", 
               date_labels="%Y", expand=expansion(mult=c(0.05, 0.1))) +
  ylab("Weekly rho") +
  scale_colour_manual(values=modType3_cols[1:4],
                      labels=c(bquote(Ens['Fcst']), 
                               bquote(Ens['Blend']),
                               bquote(Ens['Avg']),
                               bquote(Opt['Param']))) +
  theme_bw() + 
  guides(colour=guide_legend(override.aes=list(size=1), title=NULL),
         linetype=guide_legend(title=NULL),
         shape=guide_legend(title=NULL)) +
  theme(legend.position="inside",
        legend.position.inside=c(0.92, 0.24),
        legend.background=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=9),
        panel.grid.major.y=element_line(colour="grey85", linewidth=0.4),
        panel.grid.minor.y=element_line(colour="grey90", linewidth=0.2),
        axis.text=element_text(size=8))

p_weeklyB <- metric_date_df |>
  filter(metric=="rho") |>
  filter(sim %in% c("predFcst", "predBlend", "sim_avgAll", "sim_07")) |>
  mutate(lab=if_else(sim=="sim_07", "Opt['Param']", lab),
         lab=factor(lab, levels=names(modType3_cols)[1:4])) |>
  mutate(date_std=ymd("2020-12-31")+yday(date)) |>
  ggplot(aes(date_std, value, colour=lab)) + 
  geom_point(alpha=0.8, size=0.9, shape=1) + 
  geom_line(stat="smooth", method="gam", formula=y~s(x, k=10, bs="cc"), se=F) +
  scale_x_date(date_breaks="1 month", date_labels="%b") +
  ylab("Weekly rho") + 
  scale_colour_manual(values=modType3_cols[1:4],
                      labels=c(bquote(Ens['Fcst']), 
                               bquote(Ens['Blend']),
                               bquote(Ens['Avg']),
                               bquote(Opt['Param']))) +
  theme_bw() + 
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=9),
        panel.grid.major.y=element_line(colour="grey85", linewidth=0.2),
        panel.grid.minor.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        axis.text=element_text(size=8))
# p_weekly |>
#   ggsave("figs/pub_new/validation_metrics_byWeek.png", plot=_, width=6.5, height=4)

plot_grid(p_farms, 
          plot_grid(p_weeklyA, p_weeklyB, nrow=1, rel_widths=c(1, 0.5), align="h", axis="tb", labels=c("e", "f")), 
          nrow=2, rel_heights=c(1, 0.68)) |>
  ggsave("figs/pub_new/ens_rho_farm_week.png", plot=_, width=13, height=10)



# ggsave("figs/talk/rho+IDW_ensFcst.png", p_ls[[1]], width=3, height=5.5)
# ggsave("figs/talk/rho+IDW_ensBlend.png", p_ls[[2]], width=3, height=5.5)
# ggsave("figs/talk/rho+IDW_ens3D.png", p_ls[[3]], width=3, height=5.5)
# ggsave("figs/talk/rho+IDW_ens2D.png", p_ls[[4]], width=3, height=5.5)




# maps of ROC + IDW ---------------------------------------------------------

library(terra)
#ensCV_df <- read_csv("out/ensemble_CV.csv")
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_rast <- st_read("data/WeStCOMS2_meshFootprint.gpkg") |>
  rast(resolution=500)

ROC_info <- tibble(breaks=seq(0, 1, by=0.1)) |>
  mutate(break_labs=as.character(round(breaks, 1)),
         break_labs=if_else(row_number() %% 2 == 0, "", break_labs),
         letter=letters[row_number()],
         mdpt=(breaks + (lead(breaks)-breaks)/2))
farm_ROC.df <- ensCV_df |>
  group_by(sepaSite) |>
  summarise(across(starts_with("IP"),
                   ~yardstick::roc_auc_vec(.x, truth=lice_g05, event_level="second"))) |>
  inner_join(site_i)

p_ls <- vector("list", 4)
mods <- c("IP_predFf_RMSE", "IP_predBlend", "IP_sim_avgAll", "IP_sim_07")
col_labs <- c(expression(Ens['Fcst']), expression(Ens['Blend']), 
              expression(Ens['Avg']), expression(Opt['Param']))
for(i in seq_along(mods)) {
  col_lab <- col_labs[i]
  map_interp <- interpIDW(mesh_rast,
                          farm_ROC.df |>
                            rename_with(~"predColumn", .cols=matches(mods[i])) |>
                            select(easting, northing, predColumn) |>
                            drop_na() |>
                            as.matrix(),
                          radius=1000e3) |>
    mask(mesh_fp)
  farm_ROC.df_i <- farm_ROC.df |>
    rename_with(~"predColumn", .cols=matches(mods[i])) |>
    filter(!is.na(predColumn)) |>
    mutate(predColumn=pmin(pmax(predColumn, 1e-3), 1-1e-3)) |>
    select(sepaSite, predColumn, easting, northing) |>
    mutate(letter=cut(predColumn,
                      breaks=ROC_info$breaks,
                      labels=letters[1:(length(ROC_info$breaks)-1)]))
  farm_ROC_count <- farm_ROC.df_i |>
    count(letter) |>
    mutate(scaled=n/max(n)) |>
    full_join(ROC_info |> select(letter, mdpt) |> drop_na()) |>
    mutate(n=replace_na(n, 0),
           scaled=replace_na(scaled, 0)) |>
    arrange(letter)
  low_polygon <- tibble(x=c(81000, 96000, 96000, 81000)+4000,
                        y=rep(c(0, 54800/nrow(farm_ROC_count)), each=2) + 652000)
  x_rng <- diff(range(low_polygon$x))
  y_rng <- diff(range(low_polygon$y))
  farm_ROC_bar.df <- map_dfr(1:nrow(farm_ROC_count),
                             ~low_polygon |>
                               mutate(mdpt=farm_ROC_count$mdpt[.x],
                                      x=if_else(x==max(x),
                                                x,
                                                max(x)-x_rng*farm_ROC_count$scaled[.x]),
                                      y=y + (y_rng*(.x-1)))
  )
  farm_ROC_count_labs <- farm_ROC_bar.df |>
    group_by(mdpt) |>
    summarise(x=min(x), y=mean(y)) |>
    ungroup() |>
    left_join(farm_ROC_count) |>
    mutate(prop=paste0(round(n/sum(n)*100), "%"))
  
  p_ls[[i]] <- as_tibble(map_interp) |>
    bind_cols(crds(map_interp)) |>
    ggplot() +
    geom_sf(data=mesh_fp, fill="grey90", colour="grey", size=0.1) +
    geom_raster(aes(x, y, fill=lyr.1)) +
    geom_point(data=farm_ROC.df_i, aes(easting, northing, fill=predColumn),
               shape=21, size=1, stroke=0.25, colour="grey10") +
    geom_sf_text(data=locs, aes(label=Label), colour="magenta", fontface="bold", 
                 size=3.5, 
                 nudge_x=c(0, 5, -2, -8, 7, -6, 5, 2) * 1e3, 
                 nudge_y=c(10, -1, 8, 0, 3, 0, 3, -4) * 1e3) +
    geom_polygon(data=farm_ROC_bar.df, aes(x, y, fill=mdpt, group=mdpt),
                 colour="grey10", linewidth=0.15) +
    geom_text(data=farm_ROC_count_labs, aes(x, y, label=prop),
              size=2, hjust=1, vjust=0.5, nudge_x=-1000) +
    colorspace::scale_fill_binned_diverging(
      name=col_lab, palette="Blue-Red 3", rev=T,
      limits=c(0, 1), mid=0.5, breaks=ROC_info$breaks, labels=ROC_info$break_labs) +
    scale_y_continuous(limits=c(630000, 955000), oob=scales::oob_keep,
                       breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
    scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
                       limits=c(75000, 235000), oob=scales::oob_keep) +
    theme(legend.position="inside",
          legend.position.inside=c(c(0.195,0.182,0.199,0.18)[i], 0.204),
          legend.background=element_blank(),
          legend.key.height=unit(0.43, "cm"),
          legend.key.width=unit(0.0, "cm"),
          legend.text=element_text(size=7),
          legend.title=element_text(size=10, vjust=1, hjust=1),
          legend.ticks=element_line(colour="grey10", linewidth=0.25),
          legend.ticks.length=unit(0.04, "cm"),
          axis.title=element_blank())
  if(i > 1) p_ls[[i]] <- p_ls[[i]] + theme(axis.text.y=element_blank())
}

p_farms <- plot_grid(plotlist=p_ls, labels="auto", align="hv", axis="tblr", nrow=1)
# p_farms |>
#   ggsave("figs/pub_new/ens_farm-ROC+IDW_map.png", plot=_ , width=13, height=6, dpi=300)

p_weeklyA <- metric_date_df |>
  filter(metric=="'AUC'['ROC']") |>
  filter(sim %in% c("predFcst", "predBlend", "sim_avgAll", "sim_07")) |>
  mutate(lab=if_else(sim=="sim_07", "Opt['Param']", lab),
         lab=factor(lab, levels=names(modType3_cols)[1:4])) |>
  ggplot(aes(date, value, colour=lab)) + 
  geom_point(shape=1) + 
  geom_line(stat="smooth", method="gam", formula=y~s(x, k=15), se=F) +
  scale_x_date(date_breaks="1 year", #date_minor_breaks="3 months", 
               date_labels="%Y", expand=expansion(mult=c(0.05, 0.1))) +
  ylab("Weekly AUC") +
  scale_colour_manual(values=modType3_cols[1:4],
                      labels=c(bquote(Ens['Fcst']), 
                               bquote(Ens['Blend']),
                               bquote(Ens['Avg']),
                               bquote(Opt['Param']))) +
  theme_bw() + 
  guides(colour=guide_legend(override.aes=list(size=1), title=NULL),
         linetype=guide_legend(title=NULL),
         shape=guide_legend(title=NULL)) +
  theme(legend.position="inside",
        legend.position.inside=c(0.93, 0.15),
        legend.background=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=9),
        panel.grid.major.y=element_line(colour="grey85", linewidth=0.4),
        panel.grid.minor.y=element_line(colour="grey90", linewidth=0.2),
        axis.text=element_text(size=8))

p_weeklyB <- metric_date_df |>
  filter(metric=="'AUC'['ROC']") |>
  filter(sim %in% c("predFcst", "predBlend", "sim_avgAll", "sim_07")) |>
  mutate(lab=if_else(sim=="sim_07", "Opt['Param']", lab),
         lab=factor(lab, levels=names(modType3_cols)[1:4])) |>
  mutate(date_std=ymd("2020-12-31")+yday(date)) |>
  ggplot(aes(date_std, value, colour=lab)) + 
  geom_point(alpha=0.8, size=0.9, shape=1) + 
  geom_line(stat="smooth", method="gam", formula=y~s(x, k=10, bs="cc"), se=F) +
  scale_x_date(date_breaks="1 month", date_labels="%b") +
  ylab("Weekly AUC") + 
  scale_colour_manual(values=modType3_cols[1:4],
                      labels=c(bquote(Ens['Fcst']), 
                               bquote(Ens['Blend']),
                               bquote(Ens['Avg']),
                               bquote(Opt['Param']))) +
  theme_bw() + 
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=9),
        panel.grid.major.y=element_line(colour="grey85", linewidth=0.2),
        panel.grid.minor.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        axis.text=element_text(size=8))

# p_weekly |>
#   ggsave("figs/pub_new/validation_metrics_byWeek.png", plot=_, width=6.5, height=4)

plot_grid(p_farms, 
          plot_grid(p_weeklyA, p_weeklyB, nrow=1, rel_widths=c(1, 0.5), align="h", axis="tb", labels=c("e", "f")), 
          nrow=2, rel_heights=c(1, 0.68)) |>
  ggsave("figs/pub_new/ens_ROC_farm_week.png", plot=_, width=13, height=10)









# bad performers ----------------------------------------------------------

metrics_by_farm |> 
  filter(N >= 10) |>
  filter(grepl("pred", sim)) |>
  arrange(desc(rmse)) |>
  print(n=20)
metrics_by_farm |>
  ggplot(aes(rmse, sepaSite)) +
  geom_boxplot()
worst_farms <- c("AMM1", "FFMC41", "FFMC32", "SAR1", "INV1", "RTRI1", 
                 "OLD1", "FFMC27", "TMR1", "SHUI1", "FFMC59", "HEL1")
ensCV_df |> 
  filter(sepaSite %in% worst_farms) |>
  ggplot(aes(IP_predBlend, licePerFish_rtrt)) + 
  geom_point(alpha=0.5, shape=1) + 
  facet_wrap(~sepaSite, nrow=3, scales="free")

ensCV_df |> 
  filter(sepaSite %in% worst_farms) |>
  ggplot(aes(IP_predBlend^4, licePerFish_rtrt^4)) + 
  geom_point(alpha=0.5, shape=1) + 
  facet_wrap(~sepaSite, nrow=3)






# vertical distributions --------------------------------------------------


z_df <- dir("D:/sealice_ensembling/out/sim_2024-MarMay/processed/", full.names=T) |>
  map_dfr(readRDS) |>
  filter(z < 30)

sites <- sort(unique(z_df$sepaSite))
site_groups <- split(sites, ceiling(seq_along(sites)/8))

rng2wk <- ymd_h(c("2024-05-01 0", "2024-05-14 23"))
z2wk_df <- z_df |>
  filter(between(time, rng2wk[1], rng2wk[2]))

rngMay <- ymd_h(c("2024-05-01 0", "2024-05-31 23"))
zMay_df <- z_df |>
  filter(between(time, rngMay[1], rngMay[2]))

lims <- range(c(z_df$sim_avg, z_df$sim_ens_mn)^0.25)
breaks <- c(0, 0.1, 0.5, 1, 2.5, 5)
lims2wk <- range(c(z2wk_df$sim_avg, z2wk_df$sim_ens_mn)^0.25)
breaks2wk <- c(0, 0.1, 0.5, 1, 2.5, 5)
limsMay <- range(c(zMay_df$sim_avg, zMay_df$sim_ens_mn)^0.25)
breaksMay <- c(0, 0.1, 0.5, 1, 2.5, 5)
for(i in seq_along(site_groups)) {
  p <- z_df |> 
    filter(sepaSite %in% site_groups[[i]]) |>
    filter(sim_ens_mn > 0) |>
    ggplot(aes(time, z, fill=sim_ens_mn^0.25)) +
    geom_raster() +
    scale_y_reverse("Depth (m)", limits=c(30, 0)) +
    scale_x_datetime("2024", date_breaks="7 days", date_minor_breaks="1 day", 
                     date_labels="%d-%b", limits=ymd_h(c("2024-03-31 23", "2024-06-01 1"))) +
    scale_fill_viridis_c(expression(paste("Ensemble copepodids" %.% "m"^"-3" %.% "h"^"-1  ")),
                         option="turbo", end=0.95, labels=breaks, breaks=breaks^0.25,
                         limits=lims) +
    facet_wrap(~sepaSite, ncol=1, strip.position="right",
               axes="all_x", axis.labels="margins") +
    theme_classic() +
    theme(panel.grid.major.x=element_line(colour="grey50"),
          panel.grid.minor.x=element_line(colour="grey90"),
          axis.title.x=element_blank(),
          legend.position="bottom",
          legend.key.width=unit(1.5, "cm"), 
          legend.key.height=unit(0.2, "cm")) 
  ggsave(glue("figs/pub_new/vertDist_EnsBlend_by_site_{i}.png"), p, 
         height=270*ceiling(length(site_groups[[i]])/8), width=190, units="mm")
  p <- z_df |> 
    filter(sepaSite %in% site_groups[[i]]) |>
    filter(sim_avg > 0) |>
    ggplot(aes(time, z, fill=sim_avg^0.25)) +
    geom_raster() +
    scale_y_reverse("Depth (m)", limits=c(30, 0)) +
    scale_x_datetime("2024", date_breaks="7 days", date_minor_breaks="1 day", 
                     date_labels="%d-%b", limits=ymd_h(c("2024-03-31 23", "2024-06-01 1"))) +
    scale_fill_viridis_c(expression(paste("Ensemble copepodids" %.% "m"^"-3" %.% "h"^"-1  ")),
                         option="turbo", end=0.95, labels=breaks, breaks=breaks^0.25,
                         limits=lims) +
    facet_wrap(~sepaSite, ncol=1, strip.position="right",
               axes="all_x", axis.labels="margins") +
    theme_classic() +
    theme(panel.grid.major.x=element_line(colour="grey50"),
          panel.grid.minor.x=element_line(colour="grey90"),
          axis.title.x=element_blank(),
          legend.position="bottom",
          legend.key.width=unit(1.5, "cm"), 
          legend.key.height=unit(0.2, "cm")) 
  ggsave(glue("figs/pub_new/vertDist_EnsAvg_by_site_{i}.png"), p, 
         height=270*ceiling(length(site_groups[[i]])/8), width=190, units="mm")
  
  p <- z2wk_df |> 
    filter(sepaSite %in% site_groups[[i]]) |>
    filter(sim_ens_mn > 0) |>
    ggplot(aes(time, z, fill=sim_ens_mn^0.25)) +
    geom_raster() +
    scale_y_reverse("Depth (m)", limits=c(30, 0)) +
    scale_x_datetime("2024", date_breaks="2 days", date_minor_breaks="1 day", 
                     date_labels="%d-%b", limits=rng2wk) +
    scale_fill_viridis_c(expression(paste("Ensemble copepodids" %.% "m"^"-3" %.% "h"^"-1  ")),
                         option="turbo", end=0.95, breaks=breaks2wk^0.25, labels=breaks2wk,
                         limits=lims2wk) +
    facet_wrap(~sepaSite, ncol=1, strip.position="right",
               axes="all_x", axis.labels="margins") +
    theme_classic() +
    theme(panel.grid.major.x=element_line(colour="grey50"),
          panel.grid.minor.x=element_line(colour="grey90"),
          axis.title.x=element_blank(),
          legend.position="bottom",
          legend.key.width=unit(1.5, "cm"), 
          legend.key.height=unit(0.2, "cm")) 
  ggsave(glue("figs/pub_new/vertDist2wk_EnsBlend_by_site_{i}.png"), p, 
         height=270*ceiling(length(site_groups[[i]])/8), width=190, units="mm")
  p <- z2wk_df |> 
    filter(sepaSite %in% site_groups[[i]]) |>
    filter(sim_avg > 0) |>
    ggplot(aes(time, z, fill=sim_avg^0.25)) +
    geom_raster() +
    scale_y_reverse("Depth (m)", limits=c(30, 0)) +
    scale_x_datetime("2024", date_breaks="2 days", date_minor_breaks="1 day", 
                     date_labels="%d-%b", limits=rng2wk) +
    scale_fill_viridis_c(expression(paste("Ensemble copepodids" %.% "m"^"-3" %.% "h"^"-1  ")),
                         option="turbo", end=0.95, breaks=breaks2wk^0.25, labels=breaks2wk,
                         limits=lims2wk) +
    facet_wrap(~sepaSite, ncol=1, strip.position="right",
               axes="all_x", axis.labels="margins") +
    theme_classic() +
    theme(panel.grid.major.x=element_line(colour="grey50"),
          panel.grid.minor.x=element_line(colour="grey90"),
          axis.title.x=element_blank(),
          legend.position="bottom",
          legend.key.width=unit(1.5, "cm"), 
          legend.key.height=unit(0.2, "cm")) 
  ggsave(glue("figs/pub_new/vertDist2wk_EnsAvg_by_site_{i}.png"), p, 
         height=270*ceiling(length(site_groups[[i]])/8), width=190, units="mm")
  
}

illustrative_sites <- c("DHR1", "MCLN1", 
                        "GRE1",
                        "TAN2", "VUM1")
p <- z_df |> 
  filter(sepaSite %in% illustrative_sites) |>
  filter(sim_ens_mn > 0) |>
  ggplot(aes(time, z, fill=sim_ens_mn^0.25)) +
  geom_raster() +
  scale_y_reverse("Depth (m)", limits=c(30, 0)) +
  scale_x_datetime("2024", date_breaks="7 days", date_minor_breaks="1 day", 
                   date_labels="%d-%b", limits=ymd_h(c("2024-03-31 23", "2024-06-01 1"))) +
  scale_fill_viridis_c(expression(paste("Ensemble copepodids" %.% "m"^"-3" %.% "h"^"-1  ")),
                       option="turbo", end=0.95, breaks=breaks2wk^0.25, labels=breaks2wk,
                       limits=lims2wk) +
  facet_wrap(~sepaSite, ncol=1, strip.position="right",
             axes="all_x", axis.labels="margins") +
  theme_classic() +
  theme(panel.grid.major.x=element_line(colour="grey80"),
        panel.grid.minor.x=element_line(colour="grey90"),
        axis.title.x=element_blank(),
        legend.position="bottom",
        legend.key.width=unit(1.5, "cm"), 
        legend.key.height=unit(0.2, "cm")) 
ggsave(glue("figs/pub_new/vertDist_illustrative_sites.png"), p, 
       height=270, width=190, units="mm")

p <- zMay_df |> 
  filter(sepaSite %in% illustrative_sites) |>
  filter(sim_ens_mn > 0) |>
  ggplot(aes(time, z, fill=sim_ens_mn^0.25)) +
  geom_raster() +
  scale_y_reverse("Depth (m)", limits=c(30, 0)) +
  scale_x_datetime("2024", date_breaks="4 days", date_minor_breaks="1 day", 
                   date_labels="%d-%b", limits=rngMay) +
  scale_fill_viridis_c(expression(paste("Ensemble copepodids" %.% "m"^"-3" %.% "h"^"-1  ")),
                       option="turbo", end=0.95, breaks=breaksMay^0.25, labels=breaksMay,
                       limits=limsMay) +
  facet_wrap(~sepaSite, ncol=1, strip.position="right",
             axes="all_x", axis.labels="margins") +
  theme_classic() +
  theme(panel.grid.major.x=element_line(colour="grey80"),
        panel.grid.minor.x=element_line(colour="grey90"),
        axis.title.x=element_blank(),
        legend.position="bottom",
        legend.key.width=unit(1.5, "cm"), 
        legend.key.height=unit(0.2, "cm")) 
ggsave(glue("figs/pub_new/vertDist_May_illustrative_sites.png"), p, 
       height=270, width=190, units="mm", dpi=500)






# hourly animation --------------------------------------------------------

# WeStCOMS mesh
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_sf <- st_read("data/WeStCOMS2_mesh.gpkg") |> select(i, geom)
linnhe_mesh <- mesh_sf |> 
  st_crop(c(xmin=150000, xmax=220000, ymin=710000, ymax=785000))
skye_mesh <- mesh_sf |> 
  st_crop(c(xmin=100000, xmax=198000, ymin=780000, ymax=920000))
site_i <- read_csv("data/farm_sites_2023.csv") |> 
  st_as_sf(coords=c("easting", "northing"), crs=27700)

westcoms_panel <- ggplot() +
  geom_sf(data=mesh_fp, fill="grey", colour="grey30", linewidth=0.2) +
  guides(fill=guide_colourbar(title.position="top", direction="horizontal")) +
  scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), ".0\u00B0W")) +
  scale_y_continuous(breaks=c(54, 56, 58), labels=paste0(c(54, 56, 58), ".0\u00B0N")) + 
  theme_classic() +
  theme(legend.position="bottom",
        legend.key.height=unit(0.2, "cm"),
        legend.key.width=unit(1, "cm"))
linnhe_panel <- ggplot() +
  geom_sf(data=mesh_fp, fill="grey", colour="grey30", linewidth=0.2) +
  guides(fill=guide_colourbar(title.position="top", direction="horizontal")) +
  scale_x_continuous(limits=c(160000, 216000), breaks=c(-5.8, -5.4, -5)) +
  scale_y_continuous(limits=c(720000, 778000), breaks=c(56.4, 56.7)) + 
  theme_classic() +
  theme(legend.position="bottom",
        legend.key.height=unit(0.2, "cm"),
        legend.key.width=unit(1, "cm"))
skye_panel <- ggplot() +
  geom_sf(data=mesh_fp, fill="grey", colour="grey30", linewidth=0.2) +
  guides(fill=guide_colourbar(title.position="top", direction="horizontal")) +
  scale_x_continuous(limits=c(110000, 194000), breaks=c(-6.5, -6, -5.5)) +
  scale_y_continuous(limits=c(786000, 899000), breaks=c(57, 57.5)) + 
  theme_classic() +
  theme(legend.position="bottom",
        legend.key.height=unit(0.2, "cm"),
        legend.key.width=unit(1, "cm"))

# fig_temp_dir <- "~/OffAqua/sealice_ensembling/figs/temp/"
fig_temp_dir <- "D:/sealice_ensembling/figs/temp/"
f <- dirf("D:/sealice_ensembling/out/sim_2023-MarMay/processed/hourly", "Mature")[1:12]
# lims <- readRDS("out/sim_2023-MarMay/processed/hourly_Mature_pslims.rds")
# lims <- tibble(ens_mn=c(0, 1),
#                ens_CI99width=c(0, 0.5))
# lims_N <- tibble(ens_mn=c(0, 1),
#                  ens_CI99width=c(0, 0.02))
lims <- tibble(mn=c(0, 1),
               mn_N=c(0, 1),
               CI99width=c(0, 0.5^0.25),
               CI99width_N=c(0, 0.5))
breaks <- list(mn=c(0, 0.01, 0.1, 0.25, 0.5, 1),
               mn_N=seq(0, 1, by=0.25),
               CI99width=c(0, 0.01, 0.1, 0.25, 0.5),
               CI99width_N=seq(0, 0.5, by=0.1))

library(doFuture)

foreach(i=seq_along(f), 
        .options.future=list(globals=structure(TRUE, add=c("fig_temp_dir")))) %dofuture% {

  timestep <- ymd_hms("2023-03-17 00:00:00") + 
    dhours(as.numeric(str_sub(str_split_fixed(f[i], "_t_", 2)[,2], 1, -5))-1)
  
  if(file.exists(glue("{fig_temp_dir}westcoms_{format(timestep, '%F_%H')}.png"))) {
    next
  }
  
  ps_i <- readRDS(f[i]) |> filter(ens_mn > 0) |>
    mutate(ens_mn=pmin(ens_mn, lims$mn[2]))
  if(nrow(ps_i)==0) { 
    next
  }
  
  # WeStCOMS
  fig_a <- westcoms_panel + 
    geom_sf(data=mesh_sf |> inner_join(ps_i), aes(fill=ens_mn), colour=NA) + 
    geom_sf(data=site_i, colour="violet", shape=1, size=0.5) + 
    scale_fill_viridis_c("Ensemble mean cop./m2/h", option="turbo", limits=lims$mn,
                         breaks=breaks$mn^0.25, labels=breaks$mn) +
    ggtitle(format(timestep, "%b-%d %H:%M")) 
  fig_b <- westcoms_panel + 
    geom_sf(data=mesh_sf |> inner_join(ps_i), aes(fill=ens_CI99width), colour=NA) + 
    geom_sf(data=site_i, colour="violet", shape=1, size=0.5) +
    scale_fill_viridis_c("Ensemble 99% CI width", option="turbo", limits=lims$CI99width,
                         breaks=breaks$CI99width^0.25, labels=breaks$CI99width) +
    ggtitle(format(timestep, "%b-%d %H:%M"))
  ggpubr::ggarrange(fig_a, fig_b, nrow=1, common.legend=FALSE) |>
    ggsave(filename=glue("{fig_temp_dir}westcoms_{format(timestep, '%F_%H')}.png"), 
           plot=_, width=6.5, height=8)
  gc()
  fig_a <- westcoms_panel + 
    geom_sf(data=mesh_sf |> inner_join(ps_i), aes(fill=ens_mn^4), colour=NA) + 
    scale_fill_viridis_c("Ensemble mean cop./m2/h", option="turbo", limits=lims$mn_N,
                         breaks=breaks$mn_N, labels=breaks$mn_N) +
    ggtitle(format(timestep, "%b-%d %H:%M")) 
  fig_b <- westcoms_panel + 
    geom_sf(data=mesh_sf |> inner_join(ps_i), aes(fill=ens_CI99width^4), colour=NA) + 
    scale_fill_viridis_c("Ensemble 99% CI width", option="turbo", limits=lims$CI99width_N,
                         breaks=breaks$CI99width_N, labels=breaks$CI99width_N) +
    ggtitle(format(timestep, "%b-%d %H:%M"))
  ggpubr::ggarrange(fig_a, fig_b, nrow=1) |>
    ggsave(filename=glue("{fig_temp_dir}westcoms-N_{format(timestep, '%F_%H')}.png"), 
           plot=_, width=6.5, height=8)
  gc()
  
  # Linnhe
  ps_linnhe <- ps_i |> filter(i %in% linnhe_mesh$i)
  if(nrow(ps_linnhe) > 0) {
    fig_a <- linnhe_panel + 
      geom_sf(data=mesh_sf |> inner_join(ps_linnhe), aes(fill=ens_mn), colour=NA) + 
      geom_sf(data=site_i, colour="violet", shape=1, size=0.5) +
      scale_fill_viridis_c("Ensemble mean cop./m2/h", option="turbo", limits=lims$mn,
                           breaks=breaks$mn^0.25, labels=breaks$mn) +
      ggtitle(format(timestep, "%b-%d %H:%M")) 
    fig_b <- linnhe_panel + 
      geom_sf(data=mesh_sf |> inner_join(ps_linnhe), aes(fill=ens_CI99width), colour=NA) + 
      geom_sf(data=site_i, colour="violet", shape=1, size=0.5) +
      scale_fill_viridis_c("Ensemble 99% CI width", option="turbo", limits=lims$CI99width,
                           breaks=breaks$CI99width^0.25, labels=breaks$CI99width) +
      ggtitle(format(timestep, "%b-%d %H:%M"))
    ggpubr::ggarrange(fig_a, fig_b, nrow=1) |>
      ggsave(filename=glue("{fig_temp_dir}linnhe_{format(timestep, '%F_%H')}.png"), 
             plot=_, width=7, height=4.5)
    gc()
    fig_a <- linnhe_panel + 
      geom_sf(data=mesh_sf |> inner_join(ps_linnhe), aes(fill=ens_mn^4), colour=NA) + 
      geom_sf(data=site_i, colour="violet", shape=1, size=0.5) +
      scale_fill_viridis_c("Ensemble mean cop./m2/h", option="turbo", limits=lims$mn_N,
                           breaks=breaks$mn_N, labels=breaks$mn_N) +
      ggtitle(format(timestep, "%b-%d %H:%M")) 
    fig_b <- linnhe_panel + 
      geom_sf(data=mesh_sf |> inner_join(ps_linnhe), aes(fill=ens_CI99width^4), colour=NA) + 
      geom_sf(data=site_i, colour="violet", shape=1, size=0.5) +
      scale_fill_viridis_c("Ensemble 99% CI width", option="turbo", limits=lims$CI99width_N,
                           breaks=breaks$CI99width_N, labels=breaks$CI99width_N) +
      ggtitle(format(timestep, "%b-%d %H:%M"))
    ggpubr::ggarrange(fig_a, fig_b, nrow=1) |>
      ggsave(filename=glue("{fig_temp_dir}linnhe-N_{format(timestep, '%F_%H')}.png"), 
             plot=_, width=7, height=4.5)
    gc()
  }
  
  # Skye
  ps_skye <- ps_i |> filter(i %in% skye_mesh$i)
  if(nrow(ps_skye) > 0) {
    fig_a <- skye_panel + 
      geom_sf(data=mesh_sf |> inner_join(ps_skye), aes(fill=ens_mn), colour=NA) + 
      geom_sf(data=site_i, colour="violet", shape=1, size=0.5) +
      scale_fill_viridis_c("Ensemble mean cop./m2/h", option="turbo", limits=lims$mn,
                           breaks=breaks$mn^0.25, labels=breaks$mn) +
      ggtitle(format(timestep, "%b-%d %H:%M")) 
    fig_b <- skye_panel + 
      geom_sf(data=mesh_sf |> inner_join(ps_skye), aes(fill=ens_CI99width), colour=NA) + 
      geom_sf(data=site_i, colour="violet", shape=1, size=0.5) +
      scale_fill_viridis_c("Ensemble 99% CI width", option="turbo", limits=lims$CI99width,
                           breaks=breaks$CI99width^0.25, labels=breaks$CI99width) +
      ggtitle(format(timestep, "%b-%d %H:%M"))
    ggpubr::ggarrange(fig_a, fig_b, nrow=1) |>
      ggsave(filename=glue("{fig_temp_dir}skye_{format(timestep, '%F_%H')}.png"), 
             plot=_, width=8, height=6)
    gc()
    fig_a <- skye_panel + 
      geom_sf(data=mesh_sf |> inner_join(ps_skye), aes(fill=ens_mn^4), colour=NA) + 
      geom_sf(data=site_i, colour="violet", shape=1, size=0.5) +
      scale_fill_viridis_c("Ensemble mean cop./m2/h", option="turbo", limits=lims$mn_N,
                           breaks=breaks$mn_N, labels=breaks$mn_N) +
      ggtitle(format(timestep, "%b-%d %H:%M")) 
    fig_b <- skye_panel + 
      geom_sf(data=mesh_sf |> inner_join(ps_skye), aes(fill=ens_CI99width^4), colour=NA) + 
      geom_sf(data=site_i, colour="violet", shape=1, size=0.5) +
      scale_fill_viridis_c("Ensemble 99% CI width", option="turbo", limits=lims$CI99width_N,
                           breaks=breaks$CI99width_N, labels=breaks$CI99width_N) +
      ggtitle(format(timestep, "%b-%d %H:%M"))
    ggpubr::ggarrange(fig_a, fig_b, nrow=1) |>
      ggsave(filename=glue("{fig_temp_dir}/skye-N_{format(timestep, '%F_%H')}.png"), 
             plot=_, width=8, height=6)
    gc()
  }
  gc() 
}

library(av)
sets <- c("westcoms_", "linnhe_", "skye_", 
          "westcoms-N_", "linnhe-N_", "skye-N_")
for(i in sets) {
  dirf(fig_temp_dir, glue("{i}.*png")) |>
    av_encode_video(glue("figs/hourly_anim_{i}2023-MarMay.mp4"),
                      framerate=12)   
}






# site conditions ---------------------------------------------------------

f <- dir("out/siteEnv_2019-2023/sim_01", "siteConditions")

siteEnv_df <- map_dfr(f, 
                      ~read_csv(glue("out/siteEnv_2019-2023/sim_01/{.x}"), 
                                show_col_types=F, col_select=c(1,16:21)) |>
                        mutate(date=ymd(str_split_fixed(.x, "_", 3)[,2])))

siteEnv_df |> 
  pivot_longer(2:7) |> 
  ggplot(aes(date, value, group=site)) + 
  geom_line(alpha=0.1) + 
  facet_wrap(~name, scales="free_y") +
  scale_x_date(date_breaks="1 year", date_minor_breaks="3 months", date_labels="%Y") +
  scale_y_continuous(breaks=0) +
  theme(panel.grid.major=element_line(colour="grey90", linewidth=0.5),
        panel.grid.minor.x=element_line(colour="grey90", linewidth=0.2))
  
siteEnv_df |>
  group_by(date) |>
  mutate(across(2:7, ~c(scale(.x)))) |>
  pivot_longer(2:7) |> 
  ggplot(aes(date, value, group=site)) + 
  geom_line(alpha=0.1) + 
  facet_wrap(~name, scales="free_y") +
  scale_x_date(date_breaks="1 year", date_minor_breaks="3 months", date_labels="%Y") +
  scale_y_continuous(breaks=0) +
  theme(panel.grid.major=element_line(colour="grey90", linewidth=0.5),
        panel.grid.minor.x=element_line(colour="grey90", linewidth=0.2))

# There are not really 'hot' and 'cold' sites consistently, but rather sites are
# relatively hot/cold to other sites with different seasonality
siteEnv_df |>
  select(site, date, temperature) |>
  group_by(date) |>
  mutate(temperature=c(scale(temperature))) |>
  ggplot(aes(date, temperature)) + 
  geom_line() + 
  facet_wrap(~site) +
  scale_x_date(date_breaks="1 year", date_minor_breaks="3 months", date_labels="%Y") +
  scale_y_continuous(breaks=0) +
  theme(panel.grid.major=element_line(colour="grey90", linewidth=0.5),
        panel.grid.minor.x=element_line(colour="grey90", linewidth=0.2))
  
# Salinity is more consistent, with high and low salinity sites
# Mean salinity is a reasonable way to characterize sites
siteEnv_df |>
  select(site, date, salinity) |>
  group_by(date) |>
  mutate(salinity=c(scale(salinity))) |>
  ggplot(aes(date, salinity)) + 
  geom_line() + 
  facet_wrap(~site) +
  scale_x_date(date_breaks="1 year", date_minor_breaks="3 months", date_labels="%Y") +
  scale_y_continuous(breaks=0) +
  theme(panel.grid.major=element_line(colour="grey90", linewidth=0.5),
        panel.grid.minor.x=element_line(colour="grey90", linewidth=0.2))

# UV is also more consistent, with fast and slow sites
# Mean current speed is a reasonable way to characterize sites
siteEnv_df |>
  select(site, date, uv) |>
  group_by(date) |>
  mutate(uv=c(scale(uv))) |>
  ggplot(aes(date, uv)) + 
  geom_line() + 
  facet_wrap(~site) +
  scale_x_date(date_breaks="1 year", date_minor_breaks="3 months", date_labels="%Y") +
  scale_y_continuous(breaks=0) +
  theme(panel.grid.major=element_line(colour="grey90", linewidth=0.5),
        panel.grid.minor.x=element_line(colour="grey90", linewidth=0.2))

siteEnv_df |>
  select(site, date, u) |>
  group_by(date) |>
  mutate(u=c(scale(u))) |>
  ggplot(aes(date, u)) + 
  geom_line() + 
  facet_wrap(~site) +
  scale_x_date(date_breaks="1 year", date_minor_breaks="3 months", date_labels="%Y") +
  scale_y_continuous(breaks=0) +
  theme(panel.grid.major=element_line(colour="grey90", linewidth=0.5),
        panel.grid.minor.x=element_line(colour="grey90", linewidth=0.2))

siteEnv_df |>
  select(site, date, v) |>
  group_by(date) |>
  mutate(v=c(scale(v))) |>
  ggplot(aes(date, v)) + 
  geom_line() + 
  facet_wrap(~site) +
  scale_x_date(date_breaks="1 year", date_minor_breaks="3 months", date_labels="%Y") +
  scale_y_continuous(breaks=0) +
  theme(panel.grid.major=element_line(colour="grey90", linewidth=0.5),
        panel.grid.minor.x=element_line(colour="grey90", linewidth=0.2))
