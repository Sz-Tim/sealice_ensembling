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

cmr <- readRDS("../../00_misc/cmr_cmaps.RDS")

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
    tibble(sim=c("predF1", "predBlend", 
                 "sim_avg3D", "sim_avg2D", "nullTime", "nullFarm"),
           lab_short=c("Ens['Fcst']", "Ens['Blend']", 
                       "Mean3D", "Mean2D", "Null['time']", "Null['farm']"),
           lab=c("Ens['Fcst']", "Ens['Blend']", 
                 "Mean3D", "Mean2D", "Null['time']", "Null['farm']"))
  ) |>
  mutate(lab=factor(lab, 
                    levels=c("Ens['Fcst']", "Ens['Blend']", 
                             "Mean3D", "Mean2D", 
                             paste0("'3D.", 1:20, "'"), paste0("'2D.", 1:20, "'"),
                             "Null['time']", "Null['farm']")),
         lab_short=factor(lab_short, 
                          levels=c("Ens['Fcst']", "Ens['Blend']", 
                                   "Mean3D", "Mean2D", 
                                   "3D", "2D", 
                                   "Null['time']", "Null['farm']")))



# simulation settings -----------------------------------------------------

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
  mutate(across(where(is.numeric), ~signif(.x, 3))) |> 
  select(Simulation, D_h, D_hVert, eggTemp_fn, mortSal_fn, viableDegreeDays,
         swimUpSpeedMean, swimDownSpeedMean,
         salinityThresh, lightThreshNauplius, lightThreshCopepodid) |>
  write_csv("ms/table_1.csv")
  
  



# ensemble results --------------------------------------------------------

ensCV_df <- ensFull_df |> 
  select(rowNum, sepaSite, year, date, licePerFish_rtrt, lice_g05) |>
  left_join(read_csv("out/candidates/CV_candidate_predictions.csv")) |>
  left_join(read_csv("out/ensembles/CV_avg_predictions.csv")) |>
  left_join(read_csv("out/ensembles/CV_ensMix_predictions.csv") |>
              select(rowNum, IP_sLonLatD4_all) |> rename(IP_predBlend=IP_sLonLatD4_all)) |>
  left_join(read_csv("out/ensembles/CV_ensFc-1_rmse.csv") |>
              select(rowNum, .pred) |> rename(IP_predF1=.pred)) |>
  left_join(read_csv("out/ensembles/CV_ensFc-1_roc_auc.csv") |>
              select(rowNum, .pred_TRUE) |> rename(pr_predF1=.pred_TRUE))

years <- unique(ensFull_df$year)
ensNull_time <- ensNull_farm <- vector("list", length(years))
for(yr in seq_along(years)) {
  ensNull_time[[yr]] <- ensFull_df |>
    filter(year != years[yr]) |>
    mutate(week=floor(week(date)/2)) |>
    group_by(week) |>
    summarise(IP_nullTime=mean(licePerFish_rtrt)) |>
    ungroup() |>
    mutate(year=years[yr])
  ensNull_farm[[yr]] <- ensFull_df |>
    filter(year != years[yr]) |>
    group_by(sepaSite) |>
    summarise(IP_nullFarm=mean(licePerFish_rtrt)) |>
    ungroup() |>
    mutate(year=years[yr])
}
ensCV_df <- ensCV_df |>
  mutate(week=floor(week(date)/2)) |>
  left_join(reduce(ensNull_time, bind_rows)) |>
  select(-week) |>
  left_join(reduce(ensNull_farm, bind_rows))


write_csv(ensCV_df, "out/ensemble_CV.csv")



# performance plot --------------------------------------------------------

ensCV_df <- read_csv("out/ensemble_CV.csv") |>
  filter(date >= "2021-05-01") |>
  mutate(lice_g05=factor(lice_g05))

# Mean within site
metrics_by_farm <- full_join(
  ensCV_df |>
    pivot_longer(starts_with("IP_"), names_to="sim") |>
    mutate(sim=str_remove(sim, "IP_")) |>
    group_by(sepaSite, sim) |>
    summarise(rho=cor(value, licePerFish_rtrt, use="pairwise", method="spearman"),
              r=cor(value, licePerFish_rtrt, use="pairwise", method="pearson"),
              rmse=rmse_vec(value, truth=licePerFish_rtrt),
              N=n(),
              prop_g05=mean(lice_g05=="TRUE"),
              prop_0=mean(licePerFish_rtrt==0)) |>
    ungroup(),
  ensCV_df |>
    select(-starts_with("IP_predF")) |> rename_with(~str_replace(.x, "pr_pred", "IP_pred")) |>
    pivot_longer(starts_with("IP_"), names_to="sim") |>
    filter(!is.na(value)) |>
    mutate(sim=str_remove(sim, "IP_")) |>
    group_by(sepaSite, sim) |>
    summarise(ROC_AUC=roc_auc_vec(value, truth=lice_g05, event_level="second"),
              PR_AUC=average_precision_vec(value, truth=lice_g05, event_level="second")) |>
    ungroup()
)

# Mean among site
metrics_by_week <- full_join(
  ensCV_df |>
    pivot_longer(starts_with("IP_"), names_to="sim") |>
    mutate(sim=str_remove(sim, "IP_")) |>
    group_by(date, sim) |>
    summarise(rho=cor(value, licePerFish_rtrt, use="pairwise", method="spearman"),
              r=cor(value, licePerFish_rtrt, use="pairwise", method="pearson"),
              rmse=rmse_vec(value, truth=licePerFish_rtrt),
              N=n(),
              prop_g05=mean(lice_g05=="TRUE"),
              prop_0=mean(licePerFish_rtrt==0)) |>
    ungroup(),
  ensCV_df |>
    select(-starts_with("IP_predF")) |> rename_with(~str_replace(.x, "pr_pred", "IP_pred")) |>
    pivot_longer(starts_with("IP_"), names_to="sim") |>
    filter(!is.na(value)) |>
    mutate(sim=str_remove(sim, "IP_")) |>
    group_by(date, sim) |>
    summarise(ROC_AUC=roc_auc_vec(value, truth=lice_g05, event_level="second"),
              PR_AUC=average_precision_vec(value, truth=lice_g05, event_level="second")) |>
    ungroup()
)

# Medians
metrics_by_farm_mn <- metrics_by_farm |>
  filter(N >= 30) |>
  group_by(sim) |>
  summarise(rho=median(rho, na.rm=T),
            r=median(r, na.rm=T),
            rmse=median(rmse, na.rm=T),
            ROC_AUC=median(ROC_AUC, na.rm=T),
            PR_AUC=median(PR_AUC, na.rm=T),
            N=mean(N, na.rm=T),
            prop_g05=mean(prop_g05),
            prop_0=mean(prop_0, na.rm=T)) |>
  ungroup()
metrics_by_week_mn <- metrics_by_week |>
  filter(N >= 30) |>
  group_by(sim) |>
  summarise(rho=median(rho, na.rm=T),
            r=median(r, na.rm=T),
            rmse=median(rmse, na.rm=T),
            ROC_AUC=median(ROC_AUC, na.rm=T),
            PR_AUC=median(PR_AUC, na.rm=T),
            N=mean(N, na.rm=T),
            prop_g05=mean(prop_g05),
            prop_0=mean(prop_0, na.rm=T)) |>
  ungroup()

metric_ranks <- bind_rows(
  metrics_by_farm |>
    select(sepaSite, sim, N, r, rho, rmse, ROC_AUC, PR_AUC) |>
    pivot_longer(4:8, names_to="metric", values_to="value") |>
    mutate(value_lowGood=if_else(metric=="rmse", value, -value),
           type="byFarm") |>
    drop_na() |>
    group_by(sepaSite, metric) |>
    mutate(rank=min_rank(value_lowGood)) |>
    ungroup(),
  metrics_by_week |>
    select(date, sim, N, r, rho, rmse, ROC_AUC, PR_AUC) |>
    pivot_longer(4:8, names_to="metric", values_to="value") |>
    mutate(value_lowGood=if_else(metric=="rmse", value, -value),
           type="byWeek") |>
    drop_na() |>
    group_by(date, metric) |>
    mutate(rank=min_rank(value_lowGood)) |>
    ungroup()
  ) |>
  left_join(sim_i)

# metric_ranks |> 
#   # filter(N > 10) |>
#   group_by(lab, lab_short, metric, type) |> 
#   summarise(mn=median(rank, na.rm=T)) |> 
#   ungroup() |>
#   mutate(metric=factor(metric, levels=c("ROC_AUC", "PR_AUC", "rsq", "rho", "r", "rmse"),
#                        labels=c("'AUC'['ROC']", "'AUC'['PR']", "R^2", "rho", "r", "RMSE"))) |>
#   mutate(type=factor(type, 
#                      levels=c("global", "byFarm", "byWeek"),
#                      labels=c("Global", "By farm (median)", "By week (median)"))) |>
#   ggplot(aes(mn, lab)) + 
#   geom_point(aes(colour=lab_short, shape=lab_short), size=3) + 
#   geom_rug(aes(colour=lab_short), sides="b") +
#   scale_colour_manual(values=c("black", "red",
#                                scico(2, begin=0.2, end=0.7, palette="broc", direction=1),
#                                scico(2, begin=0.2, end=0.7, palette="broc", direction=1),
#                                "grey50", "grey50")) +
#   scale_shape_manual(values=c(19, 1, 5, 5, 1, 1, 3, 4)) +
#   scale_size_manual(values=c(rep(1.5, 4), rep(0.5, 2), 1, 1)) +
#   scale_alpha_manual(values=c(1, 1, 1, 1, 0.5, 0.5, 1, 1)) +
#   scale_x_continuous("Median rank", breaks=seq(1, 27, by=3),
#                      minor_breaks=seq(1, 27, by=1), limits=c(1, 27), oob=scales::oob_keep) +
#   scale_y_discrete(labels=label_parsed) +
#   facet_grid(type~metric) + 
#   theme_bw() +
#   theme(panel.grid.major=element_line(linewidth=0.5, colour="grey80"),
#         panel.grid.minor=element_line(linewidth=0.1, colour="grey95"))

all_metrics_df <- bind_rows(
  metrics_by_farm_mn |> mutate(type="byFarm"),
  metrics_by_week_mn |> mutate(type="byWeek")
) |>
  pivot_longer(any_of(c("rmse", "rsq", "rho", "r", "ROC_AUC", "PR_AUC")), names_to="metric") |>
  filter(metric %in% c("rmse", "r", "rho", "ROC_AUC")) |>
  mutate(metric=factor(metric, levels=c("ROC_AUC", "PR_AUC", "rsq", "r", "rho", "rmse"),
                       labels=c("'AUC'['ROC']", "'AUC'['PR']", "R^2", "r", "rho", "RMSE"))) |>
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
         grepl("(Null|Mean|Ens)", lab_short)) |>
  arrange(lab) |>
  mutate(label=c("Ens['Fcst']", "Ens['Blend']", "'3D'", "'2D'", "Null['time']", "Null['farm']"),
         value=seq(0.975, 0.775, length.out=6)) %>%
  bind_rows(., 
            . |> 
              filter(grepl("Mean", lab_short)) |>
              mutate(lab=c("3D.1", "2D.1"),
                     lab_short=c("3D", "2D"),
                     label=c(NA, NA)))



talk_rmse <- all_metrics_df |> filter(metric=="RMSE") |>
  metric_plot_base(theme="talk") + 
  scale_y_continuous("Median cross-validation score", limits=c(0.21, 0.4), oob=scales::oob_keep, 
                     breaks=seq(0, 1, by=0.05), minor_breaks=seq(0, 1, by=0.025))
talk_rho <- all_metrics_df |> filter(metric=="rho") |>
  metric_plot_base(theme="talk") + 
  scale_y_continuous("Median cross-validation score)", limits=c(0, 1), oob=scales::oob_keep,
                     breaks=seq(0, 1, by=0.5), minor_breaks=seq(0, 1, by=0.1)) +
  theme(axis.title.y=element_blank())
talk_auc <-  all_metrics_df |> filter(metric=="'AUC'['ROC']") |>
  metric_plot_base(theme="talk") + 
  scale_y_continuous("Median cross-validation score", limits=c(0.5, 1), oob=scales::oob_keep,
                     breaks=seq(0, 1, by=0.1), minor_breaks=seq(0, 1, by=0.05)) +
  theme(axis.title.y=element_blank())
talk_legend <- all_metrics_labs |>
  mutate(label=factor(label, levels=unique(label)),
         lab_short=factor(lab_short, levels=unique(lab_short))) |>
  ggplot() +
  geom_text(aes(type, value, label=label, colour=lab_short),
            hjust=0, nudge_x=-0.15, vjust=0.5, size=2.5, parse=T) +
  geom_point(position=position_nudge(x=-0.35),
             aes(type, value, colour=lab_short, shape=lab_short, size=lab_short), alpha=1) +
  scale_colour_manual(values=c("black", "red",
                               scico(2, begin=0.2, end=0.7, palette="broc", direction=1),
                               "grey50", "grey50",
                               scico(2, begin=0.2, end=0.7, palette="broc", direction=1))) +
  scale_shape_manual(values=c(1, 1, 5, 5, 3, 4, 1, 1)) +
  scale_size_manual(values=c(rep(2.5, 4), 1.5, 1.5, rep(1, 2))) +
  scale_alpha_manual(values=c(1, 1, 1, 1, 1, 1, 0.5, 0.5)) +
  ylim(0.575, 1.175) +
  theme(legend.position="none",
        plot.margin=margin(t=0, b=0, l=0, r=0),
        panel.border=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())
p <- plot_grid(talk_rmse, talk_auc,
          align="h", axis="tb", nrow=1, rel_widths=c(1.12,1))
ggsave("figs/talk/validation_metrics_CV_medians_sLonLatD4_all_n30.png", p, width=3.25, height=3.75)



ms_rmse <- all_metrics_df |> filter(metric=="RMSE") |>
  metric_plot_base(theme="ms") + 
  scale_y_continuous("Median cross-validation score", limits=c(0.22, 0.4), oob=scales::oob_keep, 
                     breaks=seq(0, 1, by=0.05), minor_breaks=seq(0, 1, by=0.01))
ms_rho <- all_metrics_df |> filter(metric=="rho") |>
  metric_plot_base(theme="ms") + 
  scale_y_continuous(limits=c(0, 1), oob=scales::oob_keep,
                     breaks=seq(0, 1, by=0.5), minor_breaks=seq(0, 1, by=0.1)) +
  theme(axis.title.y=element_blank())
ms_auc <-  all_metrics_df |> filter(metric=="'AUC'['ROC']") |>
  metric_plot_base(theme="ms") + 
  scale_y_continuous(limits=c(0.5, 1), oob=scales::oob_keep,
                     breaks=seq(0, 1, by=0.1), minor_breaks=seq(0, 1, by=0.05)) +
  theme(axis.title.y=element_blank())
ms_legend <- all_metrics_labs |>
  mutate(label=factor(label, levels=unique(label)),
         lab_short=factor(lab_short, levels=unique(lab_short))) |>
  ggplot() +
  geom_text(aes(type, value, label=label, colour=lab_short),
            hjust=0, nudge_x=-0.15, vjust=0.5, size=2.5, parse=T) +
  geom_point(position=position_nudge(x=-0.35),
             aes(type, value, colour=lab_short, shape=lab_short, size=lab_short), alpha=1) +
  scale_colour_manual(values=c("black", "red",
                               scico(2, begin=0.2, end=0.7, palette="broc", direction=1),
                               "grey50", "grey50",
                               scico(2, begin=0.2, end=0.7, palette="broc", direction=1))) +
  scale_shape_manual(values=c(1, 1, 5, 5, 3, 4, 1, 1)) +
  # scale_shape_manual(values=c(19, 19, 5, 5, 3, 4, 1, 1)) +
  scale_size_manual(values=c(rep(2.5, 4), 1.5, 1.5, rep(1, 2))) +
  scale_alpha_manual(values=c(1, 1, 1, 1, 1, 1, 0.5, 0.5)) +
  ylim(0.575, 1.175) +
  theme(legend.position="none",
        plot.margin=margin(t=0, b=0, l=0, r=0),
        panel.border=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())
p <- plot_grid(ms_rmse, ms_auc, ms_legend, 
               align="h", axis="tb", nrow=1, rel_widths=c(1.12,1,0.4))
ggsave("figs/pub/validation_metrics_CV_medians_sLonLatD4_all_n30.png", p, width=4.25, height=4)






# weekly performance ------------------------------------------------------

metric_date_df <- ensCV_df |>
  filter(date >= "2021-05-01") |>
  group_by(date) |>
  summarise(r=cor(IP_predBlend, licePerFish_rtrt, use="pairwise", method="spearman"),
            rmse=rmse_vec(IP_predBlend, truth=licePerFish_rtrt),
            ROC_AUC=roc_auc_vec(IP_predBlend, truth=lice_g05, event_level="second")) |>
  ungroup() |>
  mutate(sim="predBlend") |>
  bind_rows(
    ensCV_df |>
      filter(date >= "2021-05-01") |>
      group_by(date) |>
      summarise(r=cor(IP_predF1, licePerFish_rtrt, use="pairwise", method="spearman"),
                rmse=rmse_vec(IP_predF1, truth=licePerFish_rtrt),
                ROC_AUC=roc_auc_vec(pr_predF1, truth=lice_g05, event_level="second")) |>
      ungroup() |>
      mutate(sim="predF1")
  ) |>
  bind_rows(
    ensCV_df |>
      filter(date >= "2021-05-01") |>
      group_by(date) |>
      summarise(r=cor(IP_sim_avg2D, licePerFish_rtrt, use="pairwise", method="spearman"),
                rmse=rmse_vec(IP_sim_avg2D, truth=licePerFish_rtrt),
                ROC_AUC=roc_auc_vec(IP_sim_avg2D, truth=lice_g05, event_level="second")) |>
      ungroup() |>
      mutate(sim="sim_avg2D")
  ) |>
  bind_rows(
    ensCV_df |>
      filter(date >= "2021-05-01") |>
      group_by(date) |>
      summarise(r=cor(IP_sim_avg3D, licePerFish_rtrt, use="pairwise", method="spearman"),
                rmse=rmse_vec(IP_sim_avg3D, truth=licePerFish_rtrt),
                ROC_AUC=roc_auc_vec(IP_sim_avg3D, truth=lice_g05, event_level="second")) |>
      ungroup() |>
      mutate(sim="sim_avg3D")
  ) |>
  pivot_longer(2:4, names_to="metric") |>
  mutate(metric=factor(metric, levels=c("ROC_AUC", "PR_AUC", "rsq", "r", "rmse"),
                       labels=c("'AUC'['ROC']", "'AUC'['PR']", "R^2", "rho", "RMSE"))) |>
  left_join(sim_i) |>
  droplevels()

metric_date_labs <- expand_grid(date=max(metric_date_df$date),
                                metric=levels(metric_date_df$metric),
                                lab=levels(metric_date_df$lab)) |>
  mutate(value=c(0.88, 0.82, 0.745, 0.71,
                 0.85, 0.73, 0.685, 0.611,
                 0.21, 0.32, 0.305, 0.33)) |>
  filter(metric %in% c("RMSE", "'AUC'['ROC']"))

p <- metric_date_df |>
  filter(metric %in% c("RMSE", "'AUC'['ROC']")) |>
  ggplot(aes(date, value, colour=lab)) + 
  geom_point(size=0.5, shape=1) + 
  geom_line(stat="smooth", method="gam", formula=y~s(x), se=F) +
  scale_x_date(date_breaks="1 year", #date_minor_breaks="3 months", 
               date_labels="%Y", expand=expansion(mult=c(0.05, 0.1))) +
  ylab("Cross validation metric by week") +
  scale_colour_manual("Model", 
                      values=c("black", "red",
                               scico(2, begin=0.2, end=0.7, palette="broc", direction=1)), 
                      labels=c(bquote(Ens['Fcst']), 
                               bquote(Ens['Blend']), 
                               "Mean: 3D",
                               "Mean: 2D")) +
  facet_grid(metric~., scales="free_y", labeller=label_parsed) +
  theme_bw() + 
  guides(colour=guide_legend(override.aes=list(size=1))) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size=9),
        panel.grid.major.y=element_line(colour="grey85", linewidth=0.4),
        panel.grid.minor.y=element_line(colour="grey90", linewidth=0.2),
        axis.text=element_text(size=8))
ggsave("figs/pub/validation_metrics_byWeek.png", p, width=6.5, height=4)


p <- metric_date_df |>
  filter(metric %in% c("RMSE", "'AUC'['ROC']")) |>
  ggplot(aes(date, value, colour=lab, shape=lab)) + 
  geom_point(size=0.75) + 
  geom_line(stat="smooth", method="gam", formula=y~s(x), se=F) +
  # geom_line(stat="smooth", method="loess", formula=y~x, se=F, span=0.25) +
  scale_x_date(date_breaks="1 year", #date_minor_breaks="3 months", 
               date_labels="%Y") +
  ylab("Cross validation metric by week") +
  scale_colour_manual(values=c("black", "red",
                               scico(2, begin=0.2, end=0.7, palette="broc", direction=1))) +
  scale_shape_manual(values=c(19, 19, 5, 5)) +
  facet_grid(metric~., scales="free_y", labeller=label_parsed) +
  theme_bw() + 
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        panel.grid.major.y=element_line(colour="grey85", linewidth=0.4),
        panel.grid.minor.y=element_line(colour="grey90", linewidth=0.2),
        strip.text=element_text(size=12),
        axis.text=element_text(size=10))
ggsave("figs/talk/validation_metrics_byWeek.png", p, width=5.5, height=4)

metric_date_df |>
  filter(metric %in% c("RMSE", "'AUC'['ROC']")) |>
  mutate(lab=forcats::lvls_reorder(lab, c(3, 1, 2, 4))) |>
  arrange(date, metric, lab) |>
  group_by(date, metric) |>
  mutate(diff_v_3D=value - first(value)) |>
  ungroup() |>
  mutate(lab=forcats::lvls_reorder(lab, c(2, 3, 1, 4))) |>
  ggplot(aes(date, diff_v_3D, colour=lab, shape=lab)) + 
  geom_point(size=0.75) + 
  geom_line(stat="smooth", method="loess", formula=y~x, se=F, span=0.3) +
  scale_x_date(date_breaks="1 year", date_minor_breaks="3 months", 
               date_labels="%Y") +
  ylab("Cross validation metric by week") +
  scale_colour_manual(values=c("black", "red",
                               scico(2, begin=0.2, end=0.7, palette="broc", direction=1))) +
  scale_shape_manual(values=c(19, 19, 5, 5)) +
  facet_grid(metric~., scales="free_y", labeller=label_parsed) +
  theme_bw() + 
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        panel.grid.major.y=element_line(colour="grey85", linewidth=0.4),
        panel.grid.minor.y=element_line(colour="grey90", linewidth=0.2),
        strip.text=element_text(size=12),
        axis.text=element_text(size=10))



# 
# p_b <- metric_date_df |>
#   ggplot(aes(date, value, colour=metric)) +
#   geom_point(size=0.5, shape=1) + 
#   geom_line(stat="smooth", method="loess", formula=y~x, se=F, span=0.4) +
#   geom_text(data=metric_date_df |> slice_max(date), 
#             aes(date, value, label=metric),
#             hjust=0, nudge_x=20, vjust=0.5, size=2.5, parse=T) +
#   scale_x_date(date_breaks="1 year", date_minor_breaks="3 months", 
#                date_labels="%Y", expand=expansion(mult=c(0.05, 0.1))) +
#   scale_y_continuous(expression("Ens"["Blend"]~"metrics by week"), 
#                      limits=c(0, 1), breaks=c(0, 0.5, 1)) +
#   scale_colour_manual(values=colorspace::divergingx_hcl("Earth", n=6)[c(1,2,5,6)]) +
#   theme_bw() + 
#   theme(legend.position="none",
#         axis.title.x=element_blank(),
#         axis.title.y=element_text(size=9),
#         panel.grid.major.y=element_line(colour="grey85", linewidth=0.4),
#         panel.grid.minor.y=element_line(colour="grey90", linewidth=0.2),
#         axis.text=element_text(size=6))
# 
# p <- cowplot::plot_grid(p_a, p_b, nrow=2, align="hv", axis="lr", 
#                    rel_heights=c((1+sqrt(5))/2, 1), labels="auto")
# ggsave("figs/pub/validation_metrics_CV_medians_sLonLatD4_all.png", p, width=6, height=6)
# 
# 
# 
# 
# dummy_lims <- expand_grid(lab="Ens['Fcst']", 
#                           date=ymd("2022-01-01"),
#                           metric=unique(all_metrics_df$metric),
#                           lim=c("min", "max")) |>
#   mutate(value=c(0, 0.4, 0, 1, 0.5, 1, 0, 1))
# 
# p <- metrics_by_week |>
#   filter(grepl("pred", sim)) |>
#   filter(date > "2021-06-01") |>
#   select(date, sim, r, rmse, ROC_AUC, PR_AUC) |>
#   pivot_longer(3:6, names_to="metric") |>
#   mutate(metric=factor(metric, levels=c("ROC_AUC", "PR_AUC", "rsq", "r", "rmse"),
#                        labels=c("'AUC'['ROC']", "'AUC'['PR']", "R^2", "rho", "RMSE"))) |>
#   left_join(sim_i) |>
#   ggplot(aes(date, value, group=lab, colour=lab)) +
#   geom_point(data=dummy_lims, colour="white", alpha=0) +
#   geom_line(data=metrics_by_week |> filter(date > "2021-06-01" & sim=="nullFarm"),
#             aes(y=1-prop_0), group=NA, colour="grey", linewidth=0.2) +
#   geom_point(size=0.5, shape=1, alpha=0.5) +
#   geom_line(stat="smooth", method="loess", formula=y~x, se=F, span=0.4, alpha=0.5) +
#   scale_colour_manual(values=scico::scico(6, palette="glasgow", end=0.8)[1:3]) +
#   scale_x_date(date_breaks="1 year", date_minor_breaks="3 months", 
#                                               date_labels="%Y", expand=expansion(mult=c(0.05, 0.1))) +
#   scale_y_continuous("Metric value by week") +
#   facet_grid(metric~., scales="free_y", labeller=label_parsed) +
#   theme_bw() + 
#   theme(legend.position="none",
#         axis.title.x=element_blank(),
#         axis.title.y=element_text(size=9),
#         panel.grid.major.y=element_line(colour="grey85", linewidth=0.4),
#         panel.grid.minor.y=element_line(colour="grey90", linewidth=0.2),
#         axis.text=element_text(size=6))
# ggsave("figs/pub/validation_metrics_CV_dates.png", p, width=6, height=6)



# Performance summaries: values
all_metrics_df |> filter(sim=="predF1")








# farm performance --------------------------------------------------------

metric_farm_df <- ensCV_df |>
  filter(date >= "2021-05-01") |>
  group_by(sepaSite) |>
  summarise(r=cor(IP_predBlend, licePerFish_rtrt, use="pairwise", method="spearman"),
            rmse=rmse_vec(IP_predBlend, truth=licePerFish_rtrt),
            ROC_AUC=roc_auc_vec(IP_predBlend, truth=lice_g05, event_level="second")) |>
  ungroup() |>
  mutate(sim="predBlend") |>
  bind_rows(
    ensCV_df |>
      filter(date > "2021-05-01") |>
      group_by(sepaSite) |>
      summarise(r=cor(IP_predF1, licePerFish_rtrt, use="pairwise", method="spearman"),
                rmse=rmse_vec(IP_predF1, truth=licePerFish_rtrt),
                ROC_AUC=roc_auc_vec(pr_predF1, truth=lice_g05, event_level="second")) |>
      ungroup() |>
      mutate(sim="predF1")
  ) |>
  bind_rows(
    ensCV_df |>
      filter(date >= "2021-05-01") |>
      group_by(sepaSite) |>
      summarise(r=cor(IP_sim_avg2D, licePerFish_rtrt, use="pairwise", method="spearman"),
                rmse=rmse_vec(IP_sim_avg2D, truth=licePerFish_rtrt),
                ROC_AUC=roc_auc_vec(IP_sim_avg2D, truth=lice_g05, event_level="second")) |>
      ungroup() |>
      mutate(sim="sim_avg2D")
  ) |>
  bind_rows(
    ensCV_df |>
      filter(date >= "2021-05-01") |>
      group_by(sepaSite) |>
      summarise(r=cor(IP_sim_avg3D, licePerFish_rtrt, use="pairwise", method="spearman"),
                rmse=rmse_vec(IP_sim_avg3D, truth=licePerFish_rtrt),
                ROC_AUC=roc_auc_vec(IP_sim_avg3D, truth=lice_g05, event_level="second")) |>
      ungroup() |>
      mutate(sim="sim_avg3D")
  ) |>
  pivot_longer(2:4, names_to="metric") |>
  mutate(metric=factor(metric, levels=c("ROC_AUC", "PR_AUC", "rsq", "r", "rmse"),
                       labels=c("'AUC'['ROC']", "'AUC'['PR']", "R^2", "rho", "RMSE"))) |>
  left_join(sim_i) |>
  droplevels()




# Blending proportions ------------------------------------------------------

mod <- "all_sLonLatD4"
# Blending ensemble
out_ensBlend <- readRDS(glue("out/ensembles/ensMix_{mod}_FULL_stanfit.rds"))
dat_ensBlend <- readRDS(glue("out/ensembles/ensMix_{mod}_FULL_standata.rds"))

b_p_post <- rstan::extract(out_ensBlend, pars="b_p")[[1]] |>
  as_tibble(.name_repair="minimal") |>
  set_names(sim_i$lab[match(dat_ensBlend$sim_names, sim_i$sim)]) |>
  mutate(iter=row_number()) |>
  pivot_longer(-iter, names_to="Simulation", values_to="p") |>
  mutate(lab_short=str_sub(Simulation, 1, 2),
         Simulation=factor(Simulation, levels=levels(sim_i$lab)))
p_a <- ggplot(b_p_post, aes(p, Simulation, colour=lab_short, fill=lab_short)) + 
  stat_slab(normalize="xy", scale=0.7, colour=NA, fill="dodgerblue4",
            aes(slab_alpha=after_stat(-pmax(abs(1-2*cdf), 0.5))), fill_type="gradient") +
  stat_pointinterval(.width=c(0.5, 0.8, 0.95), colour="dodgerblue4") +
  scale_y_discrete(limits=rev(sort(unique(b_p_post$Simulation)))) + 
  scale_slab_alpha_continuous(range=c(0.01, 0.8)) +
  xlim(0, 1) +
  labs(x=expression(paste("Ensemble blending weight (", italic(pi[~~k]), ")")),
       y="Simulation variant") +
  theme(legend.position="none",
        axis.title=element_text(size=9),
        axis.text=element_text(size=7))


sim_key <- read_csv("out/sim_2019-2023/sim_i.csv") |> 
  mutate(sim=paste0("sim_", i)) |>
  select(sim, salinityMort, eggTemp, swimSpeed, salinityThresh, fixDepth) |> 
  inner_join(sim_i) |>
  mutate(swimSpeed=if_else(lab_short=="2D", NA, swimSpeed),
         salinityThresh=if_else(lab_short=="2D", NA, salinityThresh)) |>
  rename(Simulation=lab)
swim_post <- b_p_post |>
  inner_join(sim_key |> select(Simulation, swimSpeed)) |>
  group_by(swimSpeed, iter) |>
  summarise(p=sum(p)) |>
  ungroup() |>
  arrange(swimSpeed) |>
  mutate(swimSpeed=factor(swimSpeed, levels=c(0.001, 0.0001), 
                          labels=paste(c("0.10", "0.01"), "cm/s")))
eggTemp_post <- b_p_post |>
  inner_join(sim_key |> select(Simulation, eggTemp)) |>
  group_by(eggTemp, iter) |>
  summarise(p=sum(p)) |>
  ungroup() |>
  arrange(eggTemp) |>
  mutate(eggTemp=factor(eggTemp, levels=c(TRUE, FALSE), 
                        labels=c("f(Temp.)", "Constant")))
salMort_post <- b_p_post |>
  inner_join(sim_key |> select(Simulation, salinityMort)) |>
  group_by(salinityMort, iter) |>
  summarise(p=sum(p)) |>
  ungroup() |>
  arrange(salinityMort) |>
  mutate(salinityMort=factor(salinityMort, levels=c(TRUE, FALSE), 
                             labels=c("f(Sal.)", "Constant")))
salThresh_post <- b_p_post |>
  inner_join(sim_key |> select(Simulation, salinityThresh)) |>
  group_by(salinityThresh, iter) |>
  summarise(p=sum(p)) |>
  ungroup() |>
  arrange(salinityThresh) |>
  mutate(salinityThresh=factor(salinityThresh, levels=c("high", "low"), 
                               labels=c("23-31 psu", "20-23 psu")))
fixDepth_post <- b_p_post |>
  inner_join(sim_key |> select(Simulation, fixDepth)) |>
  group_by(fixDepth, iter) |>
  summarise(p=sum(p)) |>
  ungroup() |>
  arrange(fixDepth) |>
  mutate(fixDepth=factor(fixDepth, levels=c(TRUE, FALSE), 
                         labels=c("2D", "3D")))

p_b <- swim_post |> mutate(var="Vertical swim speed") |> rename(name=swimSpeed) |>
  bind_rows(eggTemp_post |> mutate(var="Egg production") |> rename(name=eggTemp)) |>
  bind_rows(salMort_post |> mutate(var="Larval mortality") |> rename(name=salinityMort)) |>
  bind_rows(salThresh_post |> mutate(var="Sinking threshold") |> rename(name=salinityThresh)) |>
  bind_rows(fixDepth_post |> mutate(var="Dimensionality") |> rename(name=fixDepth)) |>
  drop_na() |>
  mutate(var=factor(var, levels=c("Dimensionality", "Egg production", "Larval mortality",
                                  "Vertical swim speed", "Sinking threshold")),
         name=lvls_reorder(name, c(1, 2, 3, 5, 6, 7, 4, 9, 8))) |>
  ggplot(aes(p, name)) + 
  stat_slab(normalize="xy", scale=0.7, colour=NA, fill="dodgerblue4",
            aes(slab_alpha=after_stat(-pmax(abs(1-2*cdf), 0.5))), fill_type="gradient") +
  stat_pointinterval(.width=c(0.5, 0.8, 0.95), shape=1, colour="dodgerblue4") +
  scale_slab_alpha_continuous(range=c(0.01, 0.8)) +
  xlim(0, 1) +
  labs(x=expression(paste("Total blending weight (", italic(Sigma~~pi[~~k]), ")"))) +
  facet_grid(var~., scales="free_y", switch="y") +
  scale_y_discrete(position="right") +
  theme(legend.position="none",
        axis.title=element_text(size=9),
        axis.text=element_text(size=7),
        axis.title.y=element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(size=7))

p_full <- cowplot::plot_grid(p_a, p_b, ncol=2, labels=c("a", "b"))
ggsave(glue("figs/pub/ens_Blend_p_{mod}.png"), p_full, width=7, height=6)

# Posterior summaries
bind_rows(
  eggTemp_post |> mutate(var="eggTemp") |> rename(name=eggTemp) |> 
    drop_na() |>
    group_by(var, name) |> sevcheck::get_intervals(p),
  salMort_post |> mutate(var="salinityMort") |> rename(name=salinityMort) |> 
    drop_na() |>
    group_by(var, name) |> sevcheck::get_intervals(p),
  swim_post |> mutate(var="swimSpeed") |> rename(name=swimSpeed) |> 
    drop_na() |>
    group_by(var, name) |> sevcheck::get_intervals(p),
  salThresh_post |> mutate(var="salinityThresh") |> rename(name=salinityThresh) |> 
    drop_na() |>
    group_by(var, name) |> sevcheck::get_intervals(p),
  fixDepth_post |> mutate(var="fixDepth") |> rename(name=fixDepth) |> 
    drop_na() |>
    group_by(var, name) |> sevcheck::get_intervals(p)
)


# Comparisons
bind_rows(
  eggTemp_post |> mutate(var="eggTemp") |> 
    drop_na() |>
    group_by(var, iter) |> summarise(d=first(p)-last(p)) |>
    group_by(var) |> sevcheck::get_intervals(d),
  salMort_post |> mutate(var="salinityMort") |>
    drop_na() |>
    group_by(var, iter) |> summarise(d=first(p)-last(p)) |>
    group_by(var) |> sevcheck::get_intervals(d),
  swim_post |> mutate(var="swimSpeed") |>
    drop_na() |>
    group_by(var, iter) |> summarise(d=first(p)-last(p)) |>
    group_by(var) |> sevcheck::get_intervals(d),
  salThresh_post |> mutate(var="salinityThresh") |> 
    drop_na() |>
    group_by(var, iter) |> summarise(d=first(p)-last(p)) |>
    group_by(var) |> sevcheck::get_intervals(d),
  fixDepth_post |> mutate(var="fixDepth") |> 
    drop_na() |>
    group_by(var, iter) |> summarise(d=first(p)-last(p)) |>
    group_by(var) |> sevcheck::get_intervals(d)
)







# p post sLL --------------------------------------------------------------

# There is no universal set of proportions, so instead create a grid and 
# summarise the proportions within 10km of each site
ensFull_LatLon <- read_csv("out/valid_df.csv") |>
  select(rowNum, date, sepaSite, sepaSiteNum, licePerFish_rtrt, starts_with("sim"), starts_with("c_sim")) |>
  select(-contains("avg")) |>
  mutate(across(starts_with("sim_"), ~.x - mean(.x), .names="c_{.col}")) |>
  left_join(site_i) |>
  select(-sepaSite) |>
  arrange(rowNum)
ensBlend_rec <- readRDS("out/ensembles/recipe_sLonLatD4_all.rds")
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_bbox <- st_bbox(mesh_fp)
mesh_land <- st_convex_hull(mesh_fp) |>
  st_difference(mesh_fp) |>
  st_crop(site_i |> st_as_sf(coords=c("easting", "northing"), crs=27700) |> st_buffer(10e3))


out_ensBlend <- readRDS(glue("out/ensembles/ensMix_all_sLonLatD4_FULL_stanfit.rds"))
dat_ensBlend <- readRDS(glue("out/ensembles/ensMix_all_sLonLatD4_FULL_standata.rds"))
map_df <- expand_grid(easting=seq(min(site_i$easting)-10e3, max(site_i$easting)+10e3, by=2e3),
                      northing=seq(min(site_i$northing)-10e3, max(site_i$northing)+10e3, by=2e3)) |>
  st_as_sf(coords=c("easting", "northing"), crs=27700, remove=F) |>
  st_intersection(site_i |> st_as_sf(coords=c("easting", "northing"), crs=27700) |> st_buffer(10e3)) |>
  st_intersection(mesh_fp) |>
  st_drop_geometry() |>
  mutate(sepaSiteNum=row_number()) |>
  bind_cols(ensFull_LatLon |> summarise(across(c(licePerFish_rtrt, contains("sim")), mean))) 

b_p_ls <- make_predictions_ensMix_sLonLat(out_ensBlend, 
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
gc()
  
p_a <- ggplot(b_p_post, aes(p, Simulation, colour=lab_short, fill=lab_short)) + 
  stat_slab(normalize="xy", scale=0.7, colour=NA, fill="dodgerblue4",
            aes(slab_alpha=after_stat(-pmax(abs(1-2*cdf), 0.5))), fill_type="gradient") +
  stat_pointinterval(.width=c(0.5, 0.8, 0.95), colour="dodgerblue4") +
  scale_y_discrete(limits=rev(sort(unique(b_p_post$Simulation)))) + 
  scale_slab_alpha_continuous(range=c(0.01, 0.8)) +
  xlim(0, 1) +
  labs(x=expression(paste("Ensemble blending weight (", italic(pi[~~k]), ")")),
       y="Simulation variant") +
  theme(legend.position="none",
        axis.title=element_text(size=9),
        axis.text=element_text(size=7))

sim_key <- read_csv("out/sim_2019-2023/sim_i.csv") |> 
  mutate(sim=paste0("sim_", i)) |>
  select(sim, salinityMort, eggTemp, swimSpeed, salinityThresh, fixDepth) |> 
  inner_join(sim_i) |>
  mutate(swimSpeed=if_else(lab_short=="2D", NA, swimSpeed),
         salinityThresh=if_else(lab_short=="2D", NA, salinityThresh)) |>
  rename(Simulation=lab)
swim_post <- b_p_post |>
  inner_join(sim_key |> select(Simulation, swimSpeed)) |>
  group_by(swimSpeed, iter, rowNum) |>
  summarise(p=sum(p)) |>
  ungroup() |>
  arrange(swimSpeed) |>
  mutate(swimSpeed=factor(swimSpeed, levels=c(0.001, 0.0001), 
                          labels=paste(c("0.10", "0.01"), "cm/s")))
gc()
eggTemp_post <- b_p_post |>
  inner_join(sim_key |> select(Simulation, eggTemp)) |>
  group_by(eggTemp, iter, rowNum) |>
  summarise(p=sum(p)) |>
  ungroup() |>
  arrange(eggTemp) |>
  mutate(eggTemp=factor(eggTemp, levels=c(TRUE, FALSE), 
                        labels=c("f(Temp.)", "Constant")))
gc()
salMort_post <- b_p_post |>
  inner_join(sim_key |> select(Simulation, salinityMort)) |>
  group_by(salinityMort, iter, rowNum) |>
  summarise(p=sum(p)) |>
  ungroup() |>
  arrange(salinityMort) |>
  mutate(salinityMort=factor(salinityMort, levels=c(TRUE, FALSE), 
                             labels=c("f(Sal.)", "Constant")))
gc()
salThresh_post <- b_p_post |>
  inner_join(sim_key |> select(Simulation, salinityThresh)) |>
  group_by(salinityThresh, iter, rowNum) |>
  summarise(p=sum(p)) |>
  ungroup() |>
  arrange(salinityThresh) |>
  mutate(salinityThresh=factor(salinityThresh, levels=c("high", "low"), 
                               labels=c("23-31 psu", "20-23 psu")))
gc()
fixDepth_post <- b_p_post |>
  inner_join(sim_key |> select(Simulation, fixDepth)) |>
  group_by(fixDepth, iter, rowNum) |>
  summarise(p=sum(p)) |>
  ungroup() |>
  arrange(fixDepth) |>
  mutate(fixDepth=factor(fixDepth, levels=c(TRUE, FALSE), 
                         labels=c("2D", "3D")))

p_b <- swim_post |> mutate(var="Vertical swim speed") |> rename(name=swimSpeed) |>
  bind_rows(eggTemp_post |> mutate(var="Egg production") |> rename(name=eggTemp)) |>
  bind_rows(salMort_post |> mutate(var="Larval mortality") |> rename(name=salinityMort)) |>
  bind_rows(salThresh_post |> mutate(var="Sinking threshold") |> rename(name=salinityThresh)) |>
  bind_rows(fixDepth_post |> mutate(var="Dimensionality") |> rename(name=fixDepth)) |>
  drop_na() |>
  mutate(var=factor(var, levels=c("Dimensionality", "Egg production", "Larval mortality",
                                  "Vertical swim speed", "Sinking threshold")),
         name=lvls_reorder(name, c(1, 2, 3, 5, 6, 7, 4, 9, 8))) |>
  ggplot(aes(p, name)) + 
  stat_slab(normalize="xy", scale=0.7, colour=NA, fill="dodgerblue4",
            aes(slab_alpha=after_stat(-pmax(abs(1-2*cdf), 0.5))), fill_type="gradient") +
  stat_pointinterval(.width=c(0.5, 0.8, 0.95), shape=1, colour="dodgerblue4") +
  scale_slab_alpha_continuous(range=c(0.01, 0.8)) +
  xlim(0, 1) +
  labs(x=expression(paste("Total blending weight (", italic(Sigma~~pi[~~k]), ")"))) +
  facet_grid(var~., scales="free_y", switch="y") +
  scale_y_discrete(position="right") +
  theme(legend.position="none",
        axis.title=element_text(size=9),
        axis.text=element_text(size=7),
        axis.title.y=element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(size=7))

p_full <- cowplot::plot_grid(p_a, p_b, ncol=2, labels=c("a", "b"))
ggsave(glue("figs/pub/ens_Blend_p_sLonLatD4.png"), p_full, width=7, height=6)

# Posterior summaries
bind_rows(
  swim_post |> mutate(var="swimSpeed") |> rename(name=swimSpeed) |> 
    drop_na() |>
    group_by(var, name) |> sevcheck::get_intervals(p),
  eggTemp_post |> mutate(var="eggTemp") |> rename(name=eggTemp) |> 
    drop_na() |>
    group_by(var, name) |> sevcheck::get_intervals(p),
  salMort_post |> mutate(var="salinityMort") |> rename(name=salinityMort) |> 
    drop_na() |>
    group_by(var, name) |> sevcheck::get_intervals(p),
  salThresh_post |> mutate(var="salinityThresh") |> rename(name=salinityThresh) |> 
    drop_na() |>
    group_by(var, name) |> sevcheck::get_intervals(p),
  fixDepth_post |> mutate(var="fixDepth") |> rename(name=fixDepth) |> 
    drop_na() |>
    group_by(var, name) |> sevcheck::get_intervals(p)
) |>
  mutate(var=factor(var, levels=c("fixDepth", "eggTemp", "salinityMort",
                                  "swimSpeed", "salinityThresh")),
         name=lvls_reorder(name, c(1, 2, 3, 5, 6, 7, 4, 9, 8))) |>
  arrange(var, desc(name))


# Comparisons
bind_rows(
  eggTemp_post |> mutate(var="eggTemp") |> 
    drop_na() |>
    group_by(var, iter, rowNum) |> summarise(d=first(p)-last(p)) |>
    group_by(var) |> sevcheck::get_intervals(d),
  salMort_post |> mutate(var="salinityMort") |>
    drop_na() |>
    group_by(var, iter, rowNum) |> summarise(d=first(p)-last(p)) |>
    group_by(var) |> sevcheck::get_intervals(d),
  swim_post |> mutate(var="swimSpeed") |>
    drop_na() |>
    group_by(var, iter, rowNum) |> summarise(d=first(p)-last(p)) |>
    group_by(var) |> sevcheck::get_intervals(d),
  salThresh_post |> mutate(var="salinityThresh") |> 
    drop_na() |>
    group_by(var, iter, rowNum) |> summarise(d=first(p)-last(p)) |>
    group_by(var) |> sevcheck::get_intervals(d),
  fixDepth_post |> mutate(var="fixDepth") |> 
    drop_na() |>
    group_by(var, iter, rowNum) |> summarise(d=first(p)-last(p)) |>
    group_by(var) |> sevcheck::get_intervals(d)
)






# p maps sLL --------------------------------------------------------------


out_ensBlend <- readRDS(glue("out/ensembles/ensMix_all_sLonLatD4_FULL_stanfit.rds"))
dat_ensBlend <- readRDS(glue("out/ensembles/ensMix_all_sLonLatD4_FULL_standata.rds"))
ensBlend_rec <- readRDS("out/ensembles/recipe_sLonLatD4_all.rds")

ensFull_LatLon <- read_csv("out/valid_df.csv") |>
  select(rowNum, date, sepaSite, sepaSiteNum, licePerFish_rtrt, starts_with("sim"), starts_with("c_sim")) |>
  select(-contains("avg")) |>
  mutate(across(starts_with("sim_"), ~.x - mean(.x), .names="c_{.col}")) |>
  left_join(site_i) |>
  select(-sepaSite) |>
  arrange(rowNum)
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_bbox <- st_bbox(mesh_fp)
mesh_land <- st_convex_hull(mesh_fp) |>
  st_difference(mesh_fp) |>
  st_crop(site_i |> st_as_sf(coords=c("easting", "northing"), crs=27700) |> st_buffer(20e3))

map_df <- expand_grid(easting=seq(min(site_i$easting)-30e3, max(site_i$easting)+30e3, by=2e3),
                      northing=seq(min(site_i$northing)-30e3, max(site_i$northing)+30e3, by=2e3)) |>
  st_as_sf(coords=c("easting", "northing"), crs=27700, remove=F) |>
  # st_intersection(site_i |> st_as_sf(coords=c("easting", "northing"), crs=27700) |> st_buffer(10e3)) |>
  st_intersection(st_buffer(mesh_fp, 5e3)) |>
  st_drop_geometry() |>
  mutate(sepaSiteNum=row_number()) |>
  bind_cols(ensFull_LatLon |> summarise(across(c(licePerFish_rtrt, contains("sim")), mean))) 
p_ls <- make_predictions_ensMix_sLonLat(out_ensBlend, newdata=bake(ensBlend_rec, map_df), iter=1000, mode="b_p")
sim_details <- read_csv("out/sim_2019-2023/sim_i.csv") |>
  mutate(salinityThresh=if_else(fixDepth, "NA", salinityThresh),
         swimSpeed=if_else(fixDepth, 100, swimSpeed))
# sim_details <- sim_details |> filter(!fixDepth)
p_summary <- list(p_dim__2D=apply(p_ls[,,sim_details$fixDepth], 1:2, sum),
                  p_dim__3D=apply(p_ls[,,!sim_details$fixDepth], 1:2, sum),
                  p_salThresh__High=apply(p_ls[,,sim_details$salinityThresh=="high"], 1:2, sum),
                  p_salThresh__Low=apply(p_ls[,,sim_details$salinityThresh=="low"], 1:2, sum),
                  p_mort__Constant=apply(p_ls[,,!sim_details$salinityMort], 1:2, sum),
                  p_mort__Salinity=apply(p_ls[,,sim_details$salinityMort], 1:2, sum),
                  p_egg__Constant=apply(p_ls[,,!sim_details$eggTemp], 1:2, sum),
                  p_egg__Temperature=apply(p_ls[,,sim_details$eggTemp], 1:2, sum),
                  p_swim__Slow=apply(p_ls[,,sim_details$swimSpeed==0.0001], 1:2, sum),
                  p_swim__Fast=apply(p_ls[,,sim_details$swimSpeed==0.001], 1:2, sum))
# p_summary$p_salThresh__High <- p_summary$p_salThresh__High/p_summary$p_dim__3D
# p_summary$p_salThresh__Low <- p_summary$p_salThresh__Low/p_summary$p_dim__3D
# p_summary$p_swim__Slow <- p_summary$p_swim__Slow/p_summary$p_dim__3D
# p_summary$p_swim__Fast <- p_summary$p_swim__Fast/p_summary$p_dim__3D

p_labs <- tibble(var=c("dim", "mort", "swim", "egg", "salThresh"),
                 lab=c("Dimensions", "Larval mortality", "Vertical swim speed", 
                       "Egg production", "Sinking threshold"))

p_means <- map_df |> select(easting, northing) |>
  bind_cols(map_dfr(p_summary, rowMeans)) |>
  pivot_longer(starts_with("p_")) |>
  mutate(var=str_remove(str_split_fixed(name, "__", 2)[,1], "p_"),
         level=str_split_fixed(name, "__", 2)[,2]) |>
  mutate(level=factor(level, levels=c("2D", "3D", 
                                      "Constant", "Salinity", "Temperature",
                                      "Slow", "Fast", "Low", "High"))) |>
  select(-name) |>
  left_join(p_labs)
mean_plotlist <- p_means |>
  group_by(var) |>
  group_split() |>
  map(~.x |>
        ggplot() +
        geom_raster(aes(easting, northing, fill=value)) + 
        stat_contour(aes(easting, northing, z=value), breaks=seq(0, 1, 0.1), 
                     colour="white", linewidth=0.1) +
        geom_sf(data=mesh_land, fill="grey90", colour="grey40", linewidth=0.1) +
        geom_point(data=site_i, aes(easting, northing), shape=1, colour="black", size=0.4) +
        scale_fill_viridis_c("Posterior mean\ntotal blending weight",
                             limits=c(0, 1), begin=0.1, 
                             breaks=seq(0, 1, by=0.1), 
                             labels=c("0", "", "", "", "", "0.5", "", "", "", "", "1")) +
        scale_x_continuous(limits=range(site_i$easting)*c(0.96, 1), oob=scales::oob_keep) +
        scale_y_continuous(limits=range(site_i$northing), oob=scales::oob_keep) +
        labs(subtitle=first(.x$lab)) +
        facet_wrap(~level, nrow=1) +
        theme(axis.title=element_blank(),
              axis.text=element_blank(),
              axis.ticks=element_blank(),
              plot.subtitle=element_text(hjust=0.05)))
mean_legend <- get_legend(mean_plotlist[[1]])
mean_plotlist <- map(mean_plotlist, ~.x + theme(legend.position="none"))

mean_plot <- cowplot::plot_grid(mean_plotlist[[1]], mean_plotlist[[2]],
                                mean_plotlist[[3]], mean_plotlist[[4]],
                                mean_plotlist[[5]], mean_legend,
                                ncol=2, nrow=3, byrow=F, label_x=-0.02,
                                labels=c(LETTERS[1:5], ""), align="hv", axis="trbl")
ggsave("figs/pub/ens_Blend_p_map-mn_sLonLatD4.png", mean_plot, height=10, width=6, dpi=500)

mean_plot <- cowplot::plot_grid(mean_plotlist[[1]], mean_plotlist[[2]],
                                mean_plotlist[[3]], mean_plotlist[[4]],
                                mean_plotlist[[5]], mean_legend,
                                ncol=3, nrow=2, byrow=F, label_x=-0.02,
                                labels=c(LETTERS[1:5], ""), align="hv", axis="trbl")
ggsave("figs/talk/ens_Blend_p_map-mn_sLonLatD4.png", mean_plot, height=6.75, width=9, dpi=500)








# ensFcst explainers ------------------------------------------------------

library(DALEX); library(DALEXtra); library(tidymodels)

wf_fitted <- readRDS("out/ensembles/licePerFish_best_fitted_5wk_rmse.rds")
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
  
library(terra)
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg") 
site_10k <- site_i |>
  st_as_sf(coords=c("easting", "northing"), crs=27700) |>
  st_buffer(dist=10e3)

map_df <- ensFull_df |> 
  filter(date > ymd("2021-05-01")) |>
  summarise(across(starts_with("sim"), median)) |>
  mutate(coords=list(
    expand_grid(easting=seq(min(site_i$easting)-10e3, max(site_i$easting)+10e3, length.out=100),
                northing=seq(min(site_i$northing)-10e3, max(site_i$northing)+10e3, length.out=100),
                # date=ymd(paste0("2023-", 1:12, "-01"))))) |>
                date=seq(ymd("2023-01-01"), ymd("2023-12-31"), by="1 week")))) |>
  unnest(coords) %>%
  mutate(pred=predict(wf_fitted, new_data=.)$.pred)
map_df |>
  group_by(date) |>
  group_split() |>
  map_dfr(~.x |> select(easting, northing, pred) |>
            rast(crs="epsg:27700") |>
            mask(mesh_fp) |> 
            mask(site_10k) |>
            as.data.frame(xy=T) |>
            mutate(date=.x$date[1])) |>
  ggplot() + 
  geom_raster(aes(x, y, fill=pred)) + 
  scale_fill_viridis_c(option="turbo") +
  facet_grid(.~date) + 
  theme(legend.position="bottom")

map_df <- ensFull_df |> 
  filter(date > ymd("2021-05-01")) |>
  mutate(date=ymd("2023-01-01") + ((yday(date)-1) %/% 7)*7) |>
  group_by(date) |>
  summarise(across(starts_with("sim"), median)) |>
  ungroup() |>
  full_join(
    expand_grid(easting=seq(min(site_i$easting)-10e3, max(site_i$easting)+10e3, length.out=100),
                northing=seq(min(site_i$northing)-10e3, max(site_i$northing)+10e3, length.out=100),
                date=seq(ymd("2023-01-01"), ymd("2023-12-31"), by="1 week"))) |>
  drop_na() %>%
  mutate(pred=predict(wf_fitted, new_data=.)$.pred)
map_df |>
  group_by(date) |>
  group_split() |>
  map_dfr(~.x |> select(easting, northing, pred) |>
            rast(crs="epsg:27700") |>
            mask(mesh_fp) |> 
            mask(site_10k) |>
            as.data.frame(xy=T) |>
            mutate(date=.x$date[1])) |>
  ggplot() + 
  geom_raster(aes(x, y, fill=pred)) + 
  geom_sf(data=mesh_land, colour=NA, fill="grey40") +
  scale_fill_viridis_c(option="plasma", begin=0, end=0.95, guide="none") +
  facet_wrap(~date, nrow=4) + 
  theme(legend.position="bottom")
map_df |>
  group_by(date) |>
  group_split() |>
  map_dfr(~.x |> select(easting, northing, pred) |>
            rast(crs="epsg:27700") |>
            mask(mesh_fp) |> 
            mask(site_10k) |>
            as.data.frame(xy=T) |>
            mutate(date=.x$date[1])) |>
  group_by(x, y) |>
  mutate(pred_rel=(pred-min(pred))/(max(pred)-min(pred))) |>
  ungroup() |>
  ggplot() + 
  geom_raster(aes(x, y, fill=pred_rel)) + 
  geom_sf(data=mesh_land, colour=NA, fill="grey40") +
  scale_fill_viridis_c(option="plasma", begin=0, end=0.95, guide="none") +
  facet_wrap(~date, nrow=4) + 
  theme(legend.position="bottom")





map_df |> 
  ggplot() +
  geom_raster(aes(easting, northing, fill=pred)) + 
  geom_sf(data=mesh_fp, fill=NA) +
  scale_fill_viridis_c(option="turbo") +
  facet_wrap(~date, nrow=1)

ydayDep_df <- model_profile(wf_explain, variables=c("northing", "easting"), type="conditional", center=F)
ydayDep_df$agr_profiles |>
  ggplot(aes(`_x_`, `_yhat_`, colour=`_vname_`)) + 
  geom_point()

ydayDep_df <- model_profile(wf_explain, variables="date", center=F)
ydayDep_df$agr_profiles |>
  mutate(yday=if_else(`_vname_`=="ydayCos", acos(`_x_`), asin(`_x_`))*366) |>
  ggplot(aes(yday, `_yhat_`^4, colour=`_vname_`)) + 
  geom_line()


plot(map(str_pad(1:20, 2, "left", "0"), 
         ~model_profile(wf_explain, variables=paste0("c_sim_", .x), type="accumulated")))
pdep_ls <- map(str_pad(1:20, 2, "left", "0"), 
               ~model_profile(wf_explain, variables=paste0("c_sim_", .x)))
pdep_ls[[1]]$cp_profiles |>
  # reduce(bind_rows) |>
  ggplot(aes(x=c_sim_01, y=`_yhat_`, group=`_ids_`)) + 
  geom_line()

# scatterplots ------------------------------------------------------------

ensCV_df <- read_csv("out/ensemble_CV.csv")

preds_df <- ensCV_df |>
  select(rowNum, licePerFish_rtrt, lice_g05, IP_sim_avg2D, IP_sim_avg3D,
         IP_predBlend, IP_predF1) |>
  pivot_longer(starts_with("IP")) |>
  mutate(name=factor(name, 
                     levels=paste0("IP_", c("null", "sim_avg2D", "sim_avg3D", 
                                            "predF1", "predBlend")),
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
ggsave("figs/pub/predictions_CV_scatterplot.png", p, width=5, height=7)


map(unique(ensCV_df$sepaSite), 
    ~(ensCV_df |>
      filter(sepaSite==.x) |>
      pivot_longer(starts_with("IP_")) |> 
      ggplot(aes(value, licePerFish_rtrt)) + 
      geom_abline() + 
      geom_hline(yintercept=0.5^0.25, linetype=2) + 
      geom_vline(xintercept=0.5^0.25, linetype=2) + 
      geom_point(alpha=0.5, shape=1) + 
      facet_wrap(~name, nrow=4) + 
        xlim(0, 2) + ylim(0, 1.75) + coord_equal()) |>
      ggsave(glue("figs/siteScatter/{.x}.png"), plot=_, width=10, height=10))

map(unique(ensCV_df$date), 
    ~(ensCV_df |>
        filter(date==.x) |>
        pivot_longer(starts_with("IP_")) |> 
        ggplot(aes(value, licePerFish_rtrt)) + 
        geom_abline() + 
        geom_hline(yintercept=0.5^0.25, linetype=2) + 
        geom_vline(xintercept=0.5^0.25, linetype=2) + 
        geom_point(alpha=0.5, shape=1) + 
        facet_wrap(~name, nrow=4) + 
        xlim(0, 2) + ylim(0, 1.75) + coord_equal()) |>
      ggsave(glue("figs/dateScatter/{.x}.png"), plot=_, width=10, height=10))




# IP sLL ------------------------------------------------------------------


set.seed(1003)
out_ensBlend <- readRDS(glue("out/ensembles/ensMix_all_sLonLatD4_FULL_stanfit.rds"))
dat_ensBlend <- readRDS(glue("out/ensembles/ensMix_all_sLonLatD4_FULL_standata.rds"))

IP_LatLon <- read_csv("out/valid_df.csv") |>
  select(rowNum, date, sepaSite, sepaSiteNum, licePerFish_rtrt, starts_with("sim"), starts_with("c_sim")) |>
  select(-contains("avg")) |>
  mutate(across(starts_with("sim_"), ~.x - mean(.x), .names="c_{.col}")) |>
  group_by(sepaSite, sepaSiteNum) |>
  summarise(rowNum=first(rowNum),
            date=first(date),
            licePerFish_rtrt=mean(licePerFish_rtrt),
            across(contains("sim"), ~mean(.x, na.rm=T))) |>
  ungroup() |>
  left_join(site_i) |>
  arrange(sepaSiteNum)
ensBlend_rec <- readRDS("out/ensembles/recipe_sLonLatD4_all.rds")
b_p_ls <- make_predictions_ensMix_sLonLat(out_ensBlend, 
                                          newdata=bake(ensBlend_rec, IP_LatLon), 
                                          iter=3000, mode="b_p") 
site_p_post <- map_dfr(1:dim(b_p_ls)[2], 
                    ~as_tibble(b_p_ls[,.x,]) |>
                      set_names(dat_ensBlend$sim_names) |>
                      mutate(sepaSite=IP_LatLon$sepaSite,
                             iter=.x)) |> 
  pivot_longer(starts_with("sim"), names_to="sim", values_to="p") |>
  nest(p=c(iter, p))
rm(b_p_ls); rm(out_ensBlend); gc()

influx_df <- readRDS("out/sim_2019-2023/processed/connectivity_day.rds") |>
  select(sepaSite, sepaSite, date, sim, influx_m2) |>
  mutate(sim=paste0("sim_", sim),
         influx_m3_4rt=(replace_na(influx_m2, 0)/2)^0.25)

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
  saveRDS("out/sim_2019-2023/processed/influx_ens.rds")


influx_ens <- readRDS("out/sim_2019-2023/processed/influx_ens.rds")

thresholds <- c(0, 1e-4, 1e-3, 1e-2, 1e-1, 1)
thresh_cols <- c("white", viridis::turbo(length(thresholds)+1))

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
ggsave("figs/pub/ens_influx_daily_sLonLatD4.png", fig_influx, width=7, height=4, dpi=400)

thresholds <- c(0, 1e-3, 1e-2, 1e-1, 1)
thresh_cols <- c("white", viridis::turbo(length(thresholds)+1))

fig_influx <- influx_ens |>
  group_by(date) |>
  rename(lice=lice_mn) |>
  summarise(lt_t1=mean(lice==0),
            lt_t2=mean(lice > thresholds[1] & lice < thresholds[2]),
            lt_t3=mean(between(lice, thresholds[2], thresholds[3])),
            lt_t4=mean(between(lice, thresholds[3], thresholds[4])),
            lt_t5=mean(between(lice, thresholds[4], thresholds[5])),
            lt_t6=mean(lice > thresholds[5])) |>
  ungroup() |>
  pivot_longer(starts_with("lt_"), names_to="threshold", values_to="propSites") |>
  mutate(threshold=factor(threshold, 
                          labels=c("0", 
                                   paste(thresholds[1:4], "-", thresholds[2:5]),
                                   paste(">", thresholds[5]))),
         threshold_num=as.numeric(threshold)) |>
  filter(threshold_num != 1) |>
  ggplot(aes(date, propSites, fill=threshold_num, group=threshold_num)) +
  geom_area(colour="grey30", linewidth=0.05, outline.type="both") +
  scale_x_date(date_breaks="1 year", date_minor_breaks="3 months", date_labels="%Y") +
  scale_y_continuous("Proportion of active farms", limits=c(0,1)) +
  scale_fill_viridis_b(expression(paste("Copepodids" %.% "m"^"-2" %.% "h"^"-1")),
                       option="turbo", begin=0.05,
                       breaks=c(2.5, 3.5, 4.5, 5.5),
                       labels=c("0.001", "0.01", "0.1", "1")) +
  theme(panel.grid.major.x=element_line(colour="grey", linewidth=0.6),
        panel.grid.minor.x=element_line(colour="grey", linewidth=0.2),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=14),
        legend.position="bottom", 
        legend.title=element_text(size=12),
        legend.key.width=unit(1.5, "cm"), 
        legend.key.height=unit(0.2, "cm"),
        strip.text=element_text(size=14))
ggsave("figs/talk/ens_influx_daily_sLonLatD4.png", fig_influx, width=6, height=4, dpi=400)




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
ggsave("figs/pub/ens_influx_SpringNeap_14d.png", fig_influx, width=10, height=4, dpi=400)






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






# density sLL -------------------------------------------------------------

out_ensBlend <- readRDS(glue("out/ensembles/ensMix_all_sLonLatD4_FULL_stanfit.rds"))
dat_ensBlend <- readRDS(glue("out/ensembles/ensMix_all_sLonLatD4_FULL_standata.rds"))
ensBlend_rec <- readRDS("out/ensembles/recipe_sLonLatD4_all.rds")

# from 03a*.R -- (weekly copepodid IP)^0.25 in each grid cell
f <- dirf("out/sim_2019-2023/processed/weekly", "Mature")
ens_ls <- vector("list", length(f))
ps_lims <- tibble(ens_mn=c(0,0),
                  ens_CL005=c(0,0),
                  ens_CL025=c(0,0),
                  ens_CL975=c(0,0),
                  ens_CL995=c(0,0),
                  ens_CI95width=c(0,0),
                  ens_CI99width=c(0,0),
                  sim_sd=c(0,0),
                  ens_mn_orig=c(0,0),
                  ens_CL005_orig=c(0,0),
                  ens_CL025_orig=c(0,0),
                  ens_CL975_orig=c(0,0),
                  ens_CL995_orig=c(0,0),
                  ens_CI95width_orig=c(0,0),
                  ens_CI99width_orig=c(0,0))
p_dir <- "out/ensembles/p_meshCentroids/"

library(doFuture)
plan(multisession, workers=70)
for(i in 1:length(f)) {
  timestep <- ymd("2019-04-01") + dhours(as.numeric(str_sub(str_split_fixed(f[i], "_t_", 2)[,2], 1, -5)))
  ps_i <- readRDS(f[i]) |>
    select(i, all_of(dat_ensBlend$sim_names))
  ps_mx <- as.matrix(ps_i |> select(-i))
  
  # Calculate ensIP in parallel -- all on 4th rt scale
  ensIP <- foreach(j=1:nrow(ps_i), .combine=rbind, .inorder=TRUE, 
                   .options.future=list(globals=structure(TRUE, add=c("ps_mx", "p_dir", "ps_i")))) %dofuture% {
    ens_j <- ps_mx[j,,drop=F] %*% readRDS(glue("{p_dir}/i_{as.integer(ps_i$i[j])}.rds"))
    c(mean(ens_j), quantile(ens_j, probs=c(0.005, 0.025, 0.975, 0.995)), sd(ens_j))
  }
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
                        sim_sd=ensIP[,6]) |>
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
saveRDS(ps_lims, "out/sim_2019-2023/processed/ps_lims.rds")

timesteps <- ymd("2019-04-01") + dhours(as.numeric(str_sub(str_split_fixed(f, "_t_", 2)[,2], 1, -5)))
ens_df <- map2_dfr(ens_ls, timesteps, ~.x |> mutate(date=.y))
saveRDS(ens_df, "out/sim_2019-2023/processed/ens_weekly.rds")

ens_avg <- ens_df |>
  filter(date >= "2021-05-01") |>
  group_by(i) |>
  summarise(across(where(is.numeric), .fn=list(mn=mean, md=median))) |>
  ungroup()
saveRDS(ens_avg, "out/sim_2019-2023/processed/ens_avg_sLonLatD4_all.rds")




# maps --------------------------------------------------------------------

# Left side: Ensemble mean(copepodid density)
# Right side: Ensemble mean(weekly CI width)
ens_df <- readRDS("out/sim_2019-2023/processed/ens_weekly.rds")
ens_avg <- readRDS("out/sim_2019-2023/processed/ens_avg_sLonLatD4_all.rds")

# WeStCOMS mesh
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_sf <- st_read("data/WeStCOMS2_mesh.gpkg") |> select(i, geom)
linnhe_mesh <- mesh_sf |> 
  st_crop(c(xmin=150000, xmax=220000, ymin=710000, ymax=785000))
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
        legend.key.width=unit(0.4, "cm"),
        legend.title=element_blank(),
        legend.text=element_text(size=7),
        axis.title=element_blank())
linnhe_panel <- ggplot() +
  geom_sf(data=mesh_fp, fill="grey", colour="grey30", linewidth=0.2) +
  guides(fill=guide_colourbar(title.position="top", direction="horizontal")) +
  scale_x_continuous(limits=c(160000, 216000), breaks=c(-5.8, -5.4, -5)) +
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

mn_lims <- c(0, 0.75)
mn_breaks <- c(0, 0.01, 0.1, 0.3)

ci_lims <- c(0, 0.15)
ci_breaks <- c(0, 1e-6, 1e-4)
ci_labs <- c("0", "1e-6", "1e-4")


# mn_lims <- c(0, 0.1)
# mn_breaks <- c(0, 0.01, 0.05, 0.1)

ens_avg <- ens_avg |>
  mutate(ens_mn_mn=pmin(ens_mn_mn/2, mn_lims[2]),
         ens_CI95width_mn=pmin(ens_CI95width_mn, ci_lims[2]))

ens_map <- vector("list", 6)
# Mean densities
ens_map[[1]] <- westcoms_panel + 
  geom_sf(data=ens_avg |> right_join(mesh_sf, y=_),
          aes(fill=ens_mn_mn), colour=NA) + 
  scale_fill_viridis_c("",
                       option="turbo", limits=mn_lims,
                       breaks=mn_breaks^0.25, labels=mn_breaks) +
  annotate("text", x=79000, y=547000, label="Ensemble mean", size=3) +
  annotate("text", x=79000, y=525000, 
           label=expression("cop." %.% "m"^"-3" %.% "h"^"-1"), parse=T, size=3)
ens_map[[2]] <- linnhe_panel + 
  geom_sf(data=ens_avg |> right_join(linnhe_mesh, y=_),
          aes(fill=ens_mn_mn), colour=NA) + 
  scale_fill_viridis_c(option="turbo", limits=mn_lims,
                       breaks=mn_breaks^0.25, labels=mn_breaks) 
ens_map[[3]] <- skye_panel + 
  geom_sf(data=ens_avg |> right_join(skye_mesh, y=_),
          aes(fill=ens_mn_mn), colour=NA) + 
  scale_fill_viridis_c(option="turbo", limits=mn_lims,
                       breaks=mn_breaks^0.25, labels=mn_breaks) 
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
          rel_heights=c(2.1, 0.98, 1.23), rel_widths=c(1, 1)) |>
  ggsave("figs/pub/ens_map_sLonLatD4_all.png", plot=_, width=4.75, height=9.55, dpi=300)

ggsave("figs/talk/ens_map_WeStCOMS.png", ens_map[[1]], width=3.25, height=7, dpi=300)





# fig overview inset ------------------------------------------------------

ens_avg <- readRDS("out/sim_2019-2023/processed/ens_avg_sLonLatD4_all.rds")

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
ggsave("figs/pub/fig_overview_example_map.png", width=5.5, height=3)

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
ggsave("figs/pub/fig_overview_example_map2.png", width=5.5, height=7)




# r stormcloud ------------------------------------------------------------

# ensCV_df <- read_csv("out/ensemble_CV.csv")

# sim_09 would be selected as 'optimal'... probably
metric_ranks |> 
  filter(metric %in% c("rmse", "ROC_AUC")) |> 
  group_by(type, sim) |> 
  summarise(mnRank=mean(rank)) |> 
  group_by(sim) |> 
  summarise(mn=mean(mnRank)) |> 
  arrange(mn)

farm_r.df <- metrics_by_farm |>
  filter(N >= 30) |>
  filter(sim %in% c("predF1", "predBlend", "sim_09", "sim_17")) |>
  # filter(grepl("null|avg|pred", sim)) |>
  # filter(sim != "nullFarm") |>
  left_join(sim_i) |>
  droplevels() |>
  select(sepaSite, sim, r, rmse, ROC_AUC, PR_AUC, lab, lab_short) |>
  pivot_longer(all_of(c("r", "rmse", "ROC_AUC", "PR_AUC")), names_to="metric") |>
  mutate(type="By farm") |>
  mutate(metric=factor(metric, levels=c("ROC_AUC", "PR_AUC", "rsq", "r", "rmse"),
                     labels=c("'AUC'['ROC']", "'AUC'['PR']", "R^2", "rho", "RMSE"))) |>
  arrange(lab) |>
  mutate(type=factor(type, 
                     levels=c("global", "By farm", "By week"),
                     labels=c("Global", "'By farm'", "'By week'"))) |>
  filter(!is.na(value)) |>
  mutate(lab=lvls_revalue(lab, c("Ens['Fcst']", "Ens['Blend']", "Opt['3D']", "Opt['2D']")))
week_r.df <- metrics_by_week |>
  filter(N >= 30) |>
  filter(sim %in% c("predF1", "predBlend", "sim_09", "sim_17")) |>
  # filter(grepl("null|avg|pred", sim)) |>
  # filter(sim != "nullTime") |>
  left_join(sim_i) |>
  droplevels() |>
  select(date, sim, r, rmse, ROC_AUC, PR_AUC, lab, lab_short) |>
  pivot_longer(all_of(c("r", "rmse", "ROC_AUC", "PR_AUC")), names_to="metric") |>
  filter(!(sim=="nullTime" & metric=="ROC_AUC")) |>
  mutate(type="By week") |>
  mutate(metric=factor(metric, levels=c("ROC_AUC", "PR_AUC", "rsq", "r", "rmse"),
                       labels=c("'AUC'['ROC']", "'AUC'['PR']", "R^2", "rho", "RMSE"))) |>
  arrange(lab) |>
  mutate(type=factor(type, 
                     levels=c("global", "By farm", "By week"),
                     labels=c("Global", "'By farm'", "'By week'"))) |>
  filter(!is.na(value)) |>
  mutate(lab=lvls_revalue(lab, c("Ens['Fcst']", "Ens['Blend']", "Opt['3D']", "Opt['2D']")))
dummy_lims <- expand_grid(sim="null", 
                          lab="Null", 
                          lab_short="Null",
                          metric=unique(farm_r.df$metric),
                          lim=c("min", "max")) |>
  mutate(value=c(-1, 1, 0, 0.9, 0.4, 1, 0, 1))

mn_ci <- bind_rows(farm_r.df, week_r.df) |>
  filter(metric %in% c("'AUC'['ROC']", "RMSE")) |>
  filter(lab %in% c("Ens['Fcst']", "Ens['Blend']", "Opt['3D']", "Opt['2D']")) |>
  filter(!lab %in% c("Mean3D", "Mean2D")) |>
  droplevels() |>
  group_by(sim, lab, lab_short, metric, type) |>
  summarise(mn=mean(value, na.rm=T),
            md=median(value, na.rm=T),
            N=sum(!is.na(value)),
            se=sd(value, na.rm=T)/sqrt(N),
            ci_lo=mn - qt(0.975, N-1)*se,
            ci_hi=mn + qt(0.975, N-1)*se)

all_metrics_medians <- all_metrics_df |>
  filter(metric %in% c("'AUC'['ROC']", "RMSE")) |>
  filter(grepl("(Ens|3D|2D)", lab)) |>
  droplevels() |>
  mutate(type=paste0("'", type, "'")) |>
  mutate(lab_short=case_when(sim=="sim_avg2D" ~ "2D",
                             sim=="sim_avg3D" ~ "3D",
                             .default=lab_short),
         lab_short=factor(lab_short, levels=levels(farm_r.df$lab_short)))

# Maps among weeks are much more stable, generally better
# More variability among farms in predicting time series
p <- bind_rows(farm_r.df, week_r.df) |>
  filter(metric %in% c("'AUC'['ROC']", "RMSE")) |>
  filter(lab %in% c("Ens['Fcst']", "Ens['Blend']", "Opt['3D']", "Opt['2D']")) |>
  droplevels() |>
  mutate(metric=lvls_reorder(metric, 2:1)) |>
  ggplot(aes(value, lab, fill=lab_short, colour=lab_short)) + 
  geom_dots(side="bottom", scale=0.5) + 
  stat_slab(normalize="xy", scale=0.5, colour=NA, fill_type="gradient",
            aes(slab_alpha=after_stat(-pmax(abs(1-2*cdf), 0.25)))) +
  stat_pointinterval(.width=c(0.5, 0.8), colour="black", fatten_point=1.25) +
  geom_rug(data=all_metrics_medians |> filter(sim %in% mn_ci$sim), 
           aes(x=value), sides="b", length=unit(0.075, "npc"), linewidth=0.5, alpha=0.75) + 
  geom_rug(data=all_metrics_medians |> filter(! sim %in% mn_ci$sim), 
           aes(x=value), sides="b", length=unit(0.035, "npc"), linewidth=0.2) + 
  scale_slab_alpha_continuous(range=c(0.01, 0.75), guide="none") +
  scale_fill_manual(values=c("grey40", "red",
                             scico(2, begin=0.2, end=0.7, palette="broc", direction=1),
                             scico(2, begin=0.2, end=0.7, palette="broc", direction=1))) +
  scale_colour_manual(values=c("grey40", "red",
                               scico(2, begin=0.2, end=0.7, palette="broc", direction=1),
                               scico(2, begin=0.2, end=0.7, palette="broc", direction=1))) +
  scale_y_discrete(breaks=levels(mn_ci$lab), labels=parse(text=levels(mn_ci$lab)), limits=levels(mn_ci$lab)) +
  labs(x="Cross validation score") +
  facet_grid(type~metric, scales="free_x", labeller="label_parsed") +
  theme(panel.grid.major=element_line(colour="grey90", linewidth=0.2),
        strip.text=element_text(size=11),
        axis.title.x=element_text(size=9),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=9),
        legend.position="none")
ggsave("figs/pub/metric_stormclouds_sLonLatD4_all_n30.png", p, width=10, height=5.25, dpi=300)


p <- bind_rows(farm_r.df, week_r.df) |>
  filter(metric %in% c("'AUC'['ROC']", "RMSE")) |>
  filter(lab %in% c("Ens['Fcst']", "Ens['Blend']", "Mean3D", "Mean2D")) |>
  droplevels() |>
  mutate(metric=lvls_reorder(metric, 2:1)) |>
  ggplot(aes(value, lab, fill=lab_short, colour=lab_short)) + 
  geom_dots(side="bottom", scale=0.5) + 
  stat_slab(normalize="xy", scale=0.5, colour=NA, fill_type="gradient",
            aes(slab_alpha=after_stat(-pmax(abs(1-2*cdf), 0.5)))) +
  stat_pointinterval(.width=c(0.5, 0.8), colour="black", fatten_point=1.25) +
  geom_rug(data=mn_ci, aes(x=md), linewidth=1, sides="b", length=unit(0.05, "npc")) + 
  scale_slab_alpha_continuous(range=c(0.1, 0.5), guide="none") +
  scale_fill_manual(values=c("grey40", "red",
                               scico(2, begin=0.2, end=0.7, palette="broc", direction=1))) +
  scale_colour_manual(values=c("grey40", "red",
                               scico(2, begin=0.2, end=0.7, palette="broc", direction=1))) +
  scale_y_discrete(breaks=levels(mn_ci$lab), labels=parse(text=levels(mn_ci$lab))) +
  labs(x="Cross validation score") +
  facet_grid(type~metric, scales="free_x", labeller="label_parsed") +
  theme(panel.grid.major=element_line(colour="grey90", linewidth=0.2),
        strip.text=element_text(size=11),
        axis.title.x=element_text(size=11),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        legend.position="none")
ggsave("figs/talk/metric_stormclouds_RMSE_AUC_n30.png", p, width=6.75, height=4.5, dpi=300)



bind_rows(farm_r.df, week_r.df) |>
  filter(metric %in% c("'AUC'['ROC']", "RMSE")) |>
  filter(lab %in% c("Ens['Fcst']", "Ens['Blend']", "Mean3D", "Mean2D")) |>
  filter(lab != "Mean2D") |>
  droplevels() |>
  mutate(metric=lvls_reorder(metric, 2:1)) |>
  arrange(metric, sepaSite, date, lab) |>
  group_by(metric, sepaSite, date) |>
  mutate(diff_2D=(value - last(value))/last(value)) |>
  ungroup() |>
  ggplot(aes(diff_2D, lab, fill=lab_short, colour=lab_short)) + 
  geom_dots(side="bottom", scale=0.5) + 
  stat_slab(normalize="xy", scale=0.5, colour=NA, fill_type="gradient",
            aes(slab_alpha=after_stat(-pmax(abs(1-2*cdf), 0.5)))) +
  stat_pointinterval(.width=c(0.5, 0.8), colour="black", fatten_point=1.25) +
  # geom_rug(data=mn_ci, aes(x=md), linewidth=1, sides="b", length=unit(0.05, "npc")) + 
  scale_slab_alpha_continuous(range=c(0.1, 0.5), guide="none") +
  scale_fill_manual(values=c("grey40", "red",
                             scico(2, begin=0.2, end=0.7, palette="broc", direction=1))) +
  scale_colour_manual(values=c("grey40", "red",
                               scico(2, begin=0.2, end=0.7, palette="broc", direction=1))) +
  scale_y_discrete(breaks=levels(mn_ci$lab), labels=parse(text=levels(mn_ci$lab))) +
  labs(x="Cross validation score") +
  facet_grid(type~metric, scales="free_x", labeller="label_parsed") +
  theme(panel.grid.major=element_line(colour="grey90", linewidth=0.2),
        strip.text=element_text(size=11),
        axis.title.x=element_text(size=11),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        legend.position="none")
bind_rows(farm_r.df, week_r.df) |>
  filter(metric %in% c("'AUC'['ROC']", "RMSE")) |>
  filter(lab %in% c("Ens['Fcst']", "Ens['Blend']", "Mean3D", "Mean2D")) |>
  filter(lab != "Mean2D") |>
  droplevels() |>
  mutate(metric=lvls_reorder(metric, 2:1)) |>
  arrange(metric, sepaSite, date, lab) |>
  group_by(metric, sepaSite, date) |>
  mutate(diff_2D=(value - last(value))) |>
  ungroup() |>
  group_by(metric, lab) |>
  summarise(mnPctDiff=mean(diff_2D))


# rank stormcloud ---------------------------------------------------------

# ensCV_df <- read_csv("out/ensemble_CV.csv")

# dummy_lims <- expand_grid(sim="null", 
#                           lab="Null", 
#                           lab_short="Null",
#                           metric=unique(metric_ranks$metric),
#                           lim=c("min", "max")) |>
#   mutate(value=c(-1, 1, 0, 0.9, 0.4, 1, 0, 1))
p <- metric_ranks |>
  filter(N > 10) |>
  # filter(grepl("null|avg|pred", sim)) |>
  # filter(lab %in% c("Null['farm']", "Null['time']",
  #                   "'2D.3'",
  #                   paste0("'3D.", c(1, 3, 5, 9, 11, 13), "'"),
  #                   "Mean2D", "Mean3D", 
  #                   paste0("Ens['", c("Blend", "Fcst"), "']"))) |>
  mutate(lab=factor(lab, 
                    levels=c("Null['farm']", "Null['time']", 
                             paste0("'2D.", c(4, 2, 3, 1), "'"),
                             "Mean2D",
                             paste0("'3D.", c(16, 8, 12, 4, 14, 6, 10, 2, 15, 7, 11, 3, 13, 5, 1, 9), "'"),
                             "Mean3D",
                             paste0("Ens['", c("Blend", "Fcst"), "']")) |>
                      rev())) |>
  # filter(!grepl("5", sim)) |>
  mutate(metric=factor(metric, levels=c("ROC_AUC", "PR_AUC", "rsq", "r", "rmse"),
                       labels=c("'AUC'['ROC']", "'AUC'['PR']", "R^2", "rho", "RMSE"))) |>
  arrange(lab) |>
  mutate(type=factor(type, 
                     levels=c("global", "byFarm", "byWeek"),
                     labels=c("Global", "'By farm'", "'By week'"))) |>
  ggplot(aes(rank, lab, fill=lab_short, colour=lab_short)) + 
  # geom_point(data=dummy_lims, colour="white", size=0.05, position=position_nudge(y=0.2)) +
  geom_dots(side="bottom", scale=0.5) + 
  stat_histinterval(normalize="xy", scale=0.5, colour=NA, alpha=0.5, 
                    breaks=breaks_fixed(width=2)) +
  stat_pointinterval(.width=c(0.5, 0.95), colour="black", fatten_point=1.2) +
  scale_fill_scico_d(palette="glasgow", guide="none", end=0.8) +
  scale_colour_scico_d(palette="glasgow", guide="none", end=0.8) +
  # scale_y_discrete(breaks=levels(farm_r.df$lab), labels=parse(text=levels(farm_r.df$lab))) +
  labs(x="Rank") +
  facet_grid(type~metric, scales="free_x", labeller="label_parsed") +
  theme(panel.grid.major=element_line(colour="grey90", linewidth=0.2),
        strip.text=element_text(size=11),
        axis.title.x=element_text(size=9),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=9))
ggsave("figs/pub/metric_stormclouds_ranks_sLonLatD4_all.png", p, width=10, height=10, dpi=300)
  

metric_ranks |>
  filter(N > 10) |>
  filter(lab %in% c("Null['farm']", "Null['time']",
                    "'2D.3'",
                    paste0("'3D.", c(1, 3, 5, 9, 11, 13), "'"),
                    "Mean2D", "Mean3D",
                    paste0("Ens['", c("Blend", "Fcst"), "']"))) |>
  mutate(metric=factor(metric, levels=c("ROC_AUC", "PR_AUC", "rsq", "r", "rmse"),
                       labels=c("'AUC'['ROC']", "'AUC'['PR']", "R^2", "rho", "RMSE"))) |>
  arrange(lab) |>
  mutate(type=factor(type, 
                     levels=c("global", "byFarm", "byWeek"),
                     labels=c("Global", "'By farm'", "'By week'"))) |>
  # group_by(type, metric, date, sepaSite) |>
  # mutate(rank=min_rank(value_lowGood)) |>
  group_by(lab, type) |>
  summarise(mnRank=mean(rank),
            sdRank=sd(rank),
            q25=quantile(rank, probs=0.25),
            q75=quantile(rank, probs=0.75),
            mdRank=median(rank)) |>
  arrange(type, mnRank) |>
  print(n=28)

metric_ranks |>
  filter(N > 10) |>
  # filter(lab %in% c("Null['farm']", "Null['time']",
  #                   "'2D.3'",
  #                   paste0("'3D.", c(1, 3, 5, 9, 11, 13), "'"),
  #                   "Mean2D", "Mean3D",
  #                   paste0("Ens['", c("Blend", "Fcst"), "']"))) |>
  mutate(metric=factor(metric, levels=c("ROC_AUC", "PR_AUC", "rsq", "r", "rmse"),
                       labels=c("'AUC'['ROC']", "'AUC'['PR']", "R^2", "rho", "RMSE"))) |>
  arrange(lab) |>
  mutate(type=factor(type, 
                     levels=c("global", "byFarm", "byWeek"),
                     labels=c("Global", "'By farm'", "'By week'"))) |>
  # group_by(type, metric, date, sepaSite) |>
  # mutate(rank=min_rank(value_lowGood)) |>
  group_by(lab) |>
  summarise(mnRank=mean(rank),
            sdRank=sd(rank),
            q25=quantile(rank, probs=0.25),
            q75=quantile(rank, probs=0.75),
            mdRank=median(rank)) |>
  arrange(mnRank) |>
  print(n=28)
  


# # maps of r ---------------------------------------------------------------
# 
# mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
# ensCV_df <- read_csv("out/ensemble_CV.csv")
# 
# r_info <- tibble(breaks=seq(-1, 1, by=0.25)) |>
#   mutate(break_labs=as.character(round(breaks, 1)),
#          break_labs=if_else(row_number() %% 2 == 0, "", break_labs),
#          letter=letters[row_number()],
#          mdpt=(breaks + (lead(breaks)-breaks)/2)) 
# 
# p_ls <- vector("list", 3)
# mods <- c("IP_predF1", "IP_predF5", "IP_predBlend")
# for(i in seq_along(mods)) {
#   mod_lab <- as.character(filter(sim_i, grepl(str_sub(mods[i], 4, -1), sim))$lab) |>
#     str_remove("Ens\\['") |> str_remove("']")
#   farm_r.df <- ensCV_df |> 
#     rename_with(~"predColumn", .cols=matches(mods[i])) |>
#     group_by(sepaSite) |> 
#     summarise(r=cor(licePerFish_rtrt, predColumn, use="pairwise", method="spearman")) |> 
#     filter(!is.na(r)) |> 
#     mutate(letter=cut(r, breaks=r_info$breaks, labels=letters[1:(length(r_info$breaks)-1)])) |>
#     left_join(site_i)
#   farm_r_count <- farm_r.df |>
#     count(letter) |>
#     mutate(scaled=n/max(n)) |>
#     full_join(r_info |> select(letter, mdpt) |> drop_na()) |>
#     mutate(n=replace_na(n, 0),
#            scaled=replace_na(scaled, 0)) |>
#     arrange(letter)
#   low_polygon <- tibble(x=c(81000, 96000, 96000, 81000)+4000,
#                         y=rep(c(0, 54800/nrow(farm_r_count)), each=2) + 652000)
#   x_rng <- diff(range(low_polygon$x))
#   y_rng <- diff(range(low_polygon$y))
#   farm_r_bar.df <- map_dfr(1:nrow(farm_r_count), 
#                            ~low_polygon |> 
#                              mutate(mdpt=farm_r_count$mdpt[.x],
#                                     x=if_else(x==max(x), 
#                                               x, 
#                                               max(x)-x_rng*farm_r_count$scaled[.x]),
#                                     y=y + (y_rng*(.x-1)))
#   )
#   farm_r_count_labs <- farm_r_bar.df |>
#     group_by(mdpt) |>
#     summarise(x=min(x), y=mean(y)) |>
#     ungroup() |>
#     left_join(farm_r_count) |>
#     mutate(prop=paste0(round(n/sum(n)*100), "%"))
#   
#   col_lab <- expr(rho~":"~Ens[!!mod_lab])
#   p_ls[[i]] <- farm_r.df |> 
#     ggplot() + 
#     geom_sf(data=mesh_fp, fill="grey90", colour="grey", size=0.1) + 
#     geom_point(aes(easting, northing, fill=r), shape=21, size=2, 
#                position=position_jitter(width=2e3, height=2e3, seed=2)) +
#     geom_polygon(data=farm_r_bar.df, aes(x, y, fill=mdpt, group=mdpt), colour="grey10", linewidth=0.15) +
#     geom_text(data=farm_r_count_labs, aes(x, y, label=prop), 
#               size=2, hjust=1, vjust=0.5, nudge_x=-1000) +
#     colorspace::scale_fill_binned_diverging(
#       name=col_lab, palette="Blue-Red 3", rev=T, 
#       limits=c(-1,1), breaks=r_info$breaks, labels=r_info$break_labs,
#       l1=20, l2=90, p2=2) +
#     scale_y_continuous(limits=c(630000, 955000), oob=scales::oob_keep,
#                        breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
#     scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
#                        limits=c(75000, 235000), oob=scales::oob_keep) +
#     theme(legend.position="inside",
#           legend.position.inside=c(ifelse(i==3, 0.195, 0.185), 0.203),
#           legend.background=element_blank(),
#           legend.key.height=unit(0.395, "cm"),
#           legend.key.width=unit(0.0, "cm"),
#           legend.text=element_text(size=6),
#           legend.title=element_text(size=8, vjust=1, hjust=1),
#           legend.ticks=element_line(colour="grey10", linewidth=0.25),
#           legend.ticks.length=unit(0.04, "cm"),
#           axis.title=element_blank()) 
# }
# ggarrange(p_ls[[1]], p_ls[[2]], p_ls[[3]], nrow=1, common.legend=F, labels="auto") |> 
#   ggsave("figs/pub/ens_farm-r_map.png", plot=_ , width=9, height=5.5, dpi=300)
# 
# 
# 
# 
# 
# # maps of r + IDW ---------------------------------------------------------
# 
# library(terra)
# #ensCV_df <- read_csv("out/ensemble_CV.csv")
# mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
# mesh_rast <- st_read("data/WeStCOMS2_meshFootprint.gpkg") |>
#   rast(resolution=500)
# 
# r_info <- tibble(breaks=seq(-1, 1, by=0.25)) |>
#   mutate(break_labs=as.character(round(breaks, 1)),
#          break_labs=if_else(row_number() %% 2 == 0, "", break_labs),
#          letter=letters[row_number()],
#          mdpt=(breaks + (lead(breaks)-breaks)/2)) 
# farm_r.df <- ensCV_df |>
#   group_by(sepaSite) |>
#   summarise(across(starts_with("IP"), 
#                    ~cor(.x, licePerFish_rtrt, use="pairwise", method="spearman"))) |>
#   inner_join(site_i)
# 
# p_ls <- vector("list", 3)
# mods <- c("IP_predF1", "IP_predBlend", "IP_sim_avg3D", "IP_sim_avg2D")
# for(i in seq_along(mods)) {
#   mod_lab <- as.character(filter(sim_i, grepl(str_sub(mods[i], 4, -1), sim))$lab) |>
#     str_remove("Ens\\['") |> str_remove("']")
#   # if(grepl("Blend", mods[i])) mod_lab <- paste0(mod_lab, "  ")
#   if(grepl("Mean", mod_lab)) {
#     col_lab <- paste0("Mean", str_sub(mod_lab, -2, -1))
#   } else {
#     col_lab <- expr(Ens[!!mod_lab])
#   }
#   map_interp <- interpIDW(mesh_rast, 
#                           farm_r.df |> 
#                             rename_with(~"predColumn", .cols=matches(mods[i])) |>
#                             select(easting, northing, predColumn) |> 
#                             drop_na() |>
#                             as.matrix(),
#                           radius=1000e3) |>
#     mask(mesh_fp)
#   farm_r.df_i <- farm_r.df |> 
#     rename_with(~"predColumn", .cols=matches(mods[i])) |>
#     filter(!is.na(predColumn)) |> 
#     select(sepaSite, predColumn, easting, northing) |>
#     mutate(letter=cut(predColumn, 
#                       breaks=r_info$breaks, 
#                       labels=letters[1:(length(r_info$breaks)-1)]))
#   farm_r_count <- farm_r.df_i |>
#     count(letter) |>
#     mutate(scaled=n/max(n)) |>
#     full_join(r_info |> select(letter, mdpt) |> drop_na()) |>
#     mutate(n=replace_na(n, 0),
#            scaled=replace_na(scaled, 0)) |>
#     arrange(letter)
#   low_polygon <- tibble(x=c(81000, 96000, 96000, 81000)+4000,
#                         y=rep(c(0, 54800/nrow(farm_r_count)), each=2) + 652000)
#   x_rng <- diff(range(low_polygon$x))
#   y_rng <- diff(range(low_polygon$y))
#   farm_r_bar.df <- map_dfr(1:nrow(farm_r_count), 
#                            ~low_polygon |> 
#                              mutate(mdpt=farm_r_count$mdpt[.x],
#                                     x=if_else(x==max(x), 
#                                               x, 
#                                               max(x)-x_rng*farm_r_count$scaled[.x]),
#                                     y=y + (y_rng*(.x-1)))
#   )
#   farm_r_count_labs <- farm_r_bar.df |>
#     group_by(mdpt) |>
#     summarise(x=min(x), y=mean(y)) |>
#     ungroup() |>
#     left_join(farm_r_count) |>
#     mutate(prop=paste0(round(n/sum(n)*100), "%"))
#   
#   p_ls[[i]] <- as_tibble(map_interp) |>
#     bind_cols(crds(map_interp)) |>
#     ggplot() + 
#     geom_sf(data=mesh_fp, fill="grey90", colour="grey", size=0.1) + 
#     geom_raster(aes(x, y, fill=lyr.1)) + 
#     colorspace::scale_fill_continuous_diverging(
#       name=col_lab, palette="Blue-Red 3", rev=T, l1=20, l2=90, p2=2,
#       limits=c(-1,1), breaks=c(-1, 0, 1), guide="none") +
#     new_scale_fill() +
#     geom_point(data=farm_r.df_i, aes(easting, northing, fill=predColumn), 
#                shape=21, size=1, stroke=0.25, colour="grey10") +
#     geom_polygon(data=farm_r_bar.df, aes(x, y, fill=mdpt, group=mdpt),
#                  colour="grey10", linewidth=0.15) +
#     geom_text(data=farm_r_count_labs, aes(x, y, label=prop), 
#               size=2, hjust=1, vjust=0.5, nudge_x=-1000) +
#     colorspace::scale_fill_binned_diverging(
#       name=col_lab, palette="Blue-Red 3", rev=T, 
#       limits=c(-1,1), breaks=r_info$breaks, labels=r_info$break_labs,
#       l1=20, l2=90, p2=2) +
#     scale_y_continuous(limits=c(630000, 955000), oob=scales::oob_keep,
#                        breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
#     scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
#                        limits=c(75000, 235000), oob=scales::oob_keep) +
#     theme(legend.position="inside",
#           legend.position.inside=c(ifelse(i==3, 0.195, 0.185), 0.203),
#           legend.background=element_blank(),
#           legend.key.height=unit(0.395, "cm"),
#           legend.key.width=unit(0.0, "cm"),
#           legend.text=element_text(size=6),
#           legend.title=element_text(size=8, vjust=1, hjust=1),
#           legend.ticks=element_line(colour="grey10", linewidth=0.25),
#           legend.ticks.length=unit(0.04, "cm"),
#           axis.title=element_blank()) 
# }
# 
# ggarrange(p_ls[[1]], p_ls[[2]], p_ls[[3]], p_ls[[4]], nrow=1, common.legend=F, labels="auto") |> 
#   ggsave("figs/pub/ens_farm-r+IDW_map_sLonLatD4_all.png", plot=_ , width=12, height=6, dpi=300)
# 
# 
# ggsave("figs/talk/r+IDW_ensFcst.png", p_ls[[1]], width=3, height=5.5)
# ggsave("figs/talk/r+IDW_ensBlend.png", p_ls[[2]], width=3, height=5.5)
# ggsave("figs/talk/r+IDW_ens3D.png", p_ls[[3]], width=3, height=5.5)
# ggsave("figs/talk/r+IDW_ens2D.png", p_ls[[4]], width=3, height=5.5)
# 
# 
# 
# 
# # maps of rmse + IDW ---------------------------------------------------------
# 
# library(terra)
# #ensCV_df <- read_csv("out/ensemble_CV.csv")
# mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
# mesh_rast <- st_read("data/WeStCOMS2_meshFootprint.gpkg") |>
#   rast(resolution=500)
# 
# rmse_info <- tibble(breaks=seq(0.1, 0.65, by=0.05)) |>
#   mutate(break_labs=as.character(round(breaks, 1)),
#          break_labs=if_else(row_number() %% 2 == 0, "", break_labs),
#          letter=letters[row_number()],
#          mdpt=(breaks + (lead(breaks)-breaks)/2)) 
# farm_rmse.df <- ensCV_df |>
#   group_by(sepaSite) |>
#   summarise(across(starts_with("IP"), 
#                    ~yardstick::rmse_vec(.x, truth=licePerFish_rtrt))) |>
#   inner_join(site_i)
# 
# p_ls <- vector("list", 4)
# mods <- c("IP_predF1", "IP_predBlend", "IP_sim_avg3D", "IP_sim_avg2D")
# for(i in seq_along(mods)) {
#   mod_lab <- as.character(filter(sim_i, grepl(str_sub(mods[i], 4, -1), sim))$lab) |>
#     str_remove("Ens\\['") |> str_remove("']")
#   # if(grepl("Blend", mods[i])) mod_lab <- paste0(mod_lab, "  ")
#   if(grepl("Mean", mod_lab)) {
#     col_lab <- paste0("Mean", str_sub(mod_lab, -2, -1))
#   } else {
#     col_lab <- expr(Ens[!!mod_lab])
#   }
#   map_interp <- interpIDW(mesh_rast, 
#                           farm_rmse.df |> 
#                             rename_with(~"predColumn", .cols=matches(mods[i])) |>
#                             select(easting, northing, predColumn) |> 
#                             drop_na() |>
#                             as.matrix(),
#                           radius=1000e3) |>
#     mask(mesh_fp)
#   farm_rmse.df_i <- farm_rmse.df |> 
#     rename_with(~"predColumn", .cols=matches(mods[i])) |>
#     filter(!is.na(predColumn)) |> 
#     select(sepaSite, predColumn, easting, northing) |>
#     mutate(letter=cut(predColumn, 
#                       breaks=rmse_info$breaks, 
#                       labels=letters[1:(length(rmse_info$breaks)-1)]))
#   farm_rmse_count <- farm_rmse.df_i |>
#     count(letter) |>
#     mutate(scaled=n/max(n)) |>
#     full_join(rmse_info |> select(letter, mdpt) |> drop_na()) |>
#     mutate(n=replace_na(n, 0),
#            scaled=replace_na(scaled, 0)) |>
#     arrange(letter)
#   low_polygon <- tibble(x=c(81000, 96000, 96000, 81000)+4000,
#                         y=rep(c(0, 54800/nrow(farm_rmse_count)), each=2) + 652000)
#   x_rng <- diff(range(low_polygon$x))
#   y_rng <- diff(range(low_polygon$y))
#   farm_rmse_bar.df <- map_dfr(1:nrow(farm_rmse_count), 
#                            ~low_polygon |> 
#                              mutate(mdpt=farm_rmse_count$mdpt[.x],
#                                     x=if_else(x==max(x), 
#                                               x, 
#                                               max(x)-x_rng*farm_rmse_count$scaled[.x]),
#                                     y=y + (y_rng*(.x-1)))
#   )
#   farm_rmse_count_labs <- farm_rmse_bar.df |>
#     group_by(mdpt) |>
#     summarise(x=min(x), y=mean(y)) |>
#     ungroup() |>
#     left_join(farm_rmse_count) |>
#     mutate(prop=paste0(round(n/sum(n)*100), "%"))
#   
#   p_ls[[i]] <- as_tibble(map_interp) |>
#     bind_cols(crds(map_interp)) |>
#     ggplot() + 
#     geom_sf(data=mesh_fp, fill="grey90", colour="grey", size=0.1) + 
#     geom_raster(aes(x, y, fill=lyr.1)) + 
#     scale_fill_viridis_c(name=col_lab, option="turbo", limits=c(0.1, 0.65), guide="none") +
#     # colorspace::scale_fill_continuous_diverging(
#     #   name=col_lab, palette="Blue-Red 3", rev=F, l1=20, l2=90, p2=2, mid=0.5,
#     #   limits=c(0.1,0.6), guide="none") +
#     new_scale_fill() +
#     geom_point(data=farm_rmse.df_i, aes(easting, northing, fill=predColumn), 
#                shape=21, size=1, stroke=0.25, colour="grey10") +
#     # geom_polygon(data=farm_rmse_bar.df, aes(x, y, fill=mdpt, group=mdpt),
#     #              colour="grey10", linewidth=0.15) +
#     # geom_text(data=farm_rmse_count_labs, aes(x, y, label=prop), 
#     #           size=2, hjust=1, vjust=0.5, nudge_x=-1000) +
#     scale_fill_viridis_b(name=col_lab, option="turbo", limits=c(0.1, 0.65),
#                          breaks=rmse_info$breaks, labels=rmse_info$break_labs) +
#     # colorspace::scale_fill_binned_diverging(
#     #   name=col_lab, palette="Blue-Red 3", rev=F, mid=0.5,
#     #   limits=c(0.1,0.6), breaks=rmse_info$breaks, labels=rmse_info$break_labs,
#     #   l1=20, l2=90, p2=2) +
#     scale_y_continuous(limits=c(630000, 955000), oob=scales::oob_keep,
#                        breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
#     scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
#                        limits=c(75000, 235000), oob=scales::oob_keep) +
#     theme(legend.position="inside",
#           legend.position.inside=c(ifelse(i==3, 0.195, 0.185), 0.203),
#           legend.background=element_blank(),
#           legend.key.height=unit(0.395, "cm"),
#           legend.key.width=unit(0.0, "cm"),
#           legend.text=element_text(size=6),
#           legend.title=element_text(size=8, vjust=1, hjust=1),
#           legend.ticks=element_line(colour="grey10", linewidth=0.25),
#           legend.ticks.length=unit(0.04, "cm"),
#           axis.title=element_blank()) 
#   if(i > 1) p_ls[[i]] <- p_ls[[i]] + theme(axis.text.y=element_blank())
# }
# 
# ggarrange(p_ls[[1]], p_ls[[2]], p_ls[[3]], p_ls[[4]], nrow=1, common.legend=F, labels="auto") |> 
#   ggsave("figs/pub/ens_farm-rmse+IDW_map_sLonLatD4_all.png", plot=_ , width=12, height=6, dpi=300)
# 
# 
# ggsave("figs/talk/rmse+IDW_ensFcst.png", p_ls[[1]], width=3, height=5.5)
# ggsave("figs/talk/rmse+IDW_ensBlend.png", p_ls[[2]], width=3, height=5.5)
# ggsave("figs/talk/rmse+IDW_ens3D.png", p_ls[[3]], width=3, height=5.5)
# ggsave("figs/talk/rmse+IDW_ens2D.png", p_ls[[4]], width=3, height=5.5)
# 
# 
# 
# # maps of mae + IDW ---------------------------------------------------------
# 
# library(terra)
# #ensCV_df <- read_csv("out/ensemble_CV.csv")
# mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
# mesh_rast <- st_read("data/WeStCOMS2_meshFootprint.gpkg") |>
#   rast(resolution=500)
# 
# mae_info <- tibble(breaks=seq(0.05, 0.55, by=0.1)) |>
#   mutate(break_labs=as.character(round(breaks, 1)),
#          break_labs=if_else(row_number() %% 2 == 0, "", break_labs),
#          letter=letters[row_number()],
#          mdpt=(breaks + (lead(breaks)-breaks)/2)) 
# farm_mae.df <- ensCV_df |>
#   group_by(sepaSite) |>
#   summarise(across(starts_with("IP"), 
#                    ~yardstick::mae_vec(.x, truth=licePerFish_rtrt))) |>
#   inner_join(site_i)
# 
# p_ls <- vector("list", 4)
# mods <- c("IP_predF1", "IP_predBlend", "IP_sim_avg3D", "IP_sim_avg2D")
# for(i in seq_along(mods)) {
#   mod_lab <- as.character(filter(sim_i, grepl(str_sub(mods[i], 4, -1), sim))$lab) |>
#     str_remove("Ens\\['") |> str_remove("']")
#   # if(grepl("Blend", mods[i])) mod_lab <- paste0(mod_lab, "  ")
#   if(grepl("Mean", mod_lab)) {
#     col_lab <- paste0("Mean", str_sub(mod_lab, -2, -1))
#   } else {
#     col_lab <- expr(Ens[!!mod_lab])
#   }
#   map_interp <- interpIDW(mesh_rast, 
#                           farm_mae.df |> 
#                             rename_with(~"predColumn", .cols=matches(mods[i])) |>
#                             select(easting, northing, predColumn) |> 
#                             drop_na() |>
#                             as.matrix(),
#                           radius=1000e3) |>
#     mask(mesh_fp)
#   farm_mae.df_i <- farm_mae.df |> 
#     rename_with(~"predColumn", .cols=matches(mods[i])) |>
#     filter(!is.na(predColumn)) |> 
#     select(sepaSite, predColumn, easting, northing) |>
#     mutate(letter=cut(predColumn, 
#                       breaks=mae_info$breaks, 
#                       labels=letters[1:(length(mae_info$breaks)-1)]))
#   farm_mae_count <- farm_mae.df_i |>
#     count(letter) |>
#     mutate(scaled=n/max(n)) |>
#     full_join(mae_info |> select(letter, mdpt) |> drop_na()) |>
#     mutate(n=replace_na(n, 0),
#            scaled=replace_na(scaled, 0)) |>
#     arrange(letter)
#   low_polygon <- tibble(x=c(81000, 96000, 96000, 81000)+4000,
#                         y=rep(c(0, 54800/nrow(farm_mae_count)), each=2) + 652000)
#   x_rng <- diff(range(low_polygon$x))
#   y_rng <- diff(range(low_polygon$y))
#   farm_mae_bar.df <- map_dfr(1:nrow(farm_mae_count), 
#                               ~low_polygon |> 
#                                 mutate(mdpt=farm_mae_count$mdpt[.x],
#                                        x=if_else(x==max(x), 
#                                                  x, 
#                                                  max(x)-x_rng*farm_mae_count$scaled[.x]),
#                                        y=y + (y_rng*(.x-1)))
#   )
#   farm_mae_count_labs <- farm_mae_bar.df |>
#     group_by(mdpt) |>
#     summarise(x=min(x), y=mean(y)) |>
#     ungroup() |>
#     left_join(farm_mae_count) |>
#     mutate(prop=paste0(round(n/sum(n)*100), "%"))
#   
#   p_ls[[i]] <- as_tibble(map_interp) |>
#     bind_cols(crds(map_interp)) |>
#     ggplot() + 
#     geom_sf(data=mesh_fp, fill="grey90", colour="grey", size=0.1) + 
#     geom_raster(aes(x, y, fill=lyr.1)) + 
#     # colorspace::scale_fill_continuous_sequential(
#     #   name=col_lab, palette="Reds", rev=T,
#     #   limits=c(0.1, 0.9), guide="none") +
#     colorspace::scale_fill_continuous_diverging(
#       name=col_lab, palette="Blue-Red 3", rev=F, l1=20, l2=90, p2=2, mid=0.5,
#       limits=c(0.05,0.55), guide="none") +
#     new_scale_fill() +
#     geom_point(data=farm_mae.df_i, aes(easting, northing, fill=predColumn), 
#                shape=21, size=1, stroke=0.25, colour="grey10") +
#     geom_polygon(data=farm_mae_bar.df, aes(x, y, fill=mdpt, group=mdpt),
#                  colour="grey10", linewidth=0.15) +
#     geom_text(data=farm_mae_count_labs, aes(x, y, label=prop), 
#               size=2, hjust=1, vjust=0.5, nudge_x=-1000) +
#     colorspace::scale_fill_binned_diverging(
#       name=col_lab, palette="Blue-Red 3", rev=F, mid=0.5,
#       limits=c(0.05,0.55), breaks=mae_info$breaks, labels=mae_info$break_labs,
#       l1=20, l2=90, p2=2) +
#     scale_y_continuous(limits=c(630000, 955000), oob=scales::oob_keep,
#                        breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
#     scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
#                        limits=c(75000, 235000), oob=scales::oob_keep) +
#     theme(legend.position="inside",
#           legend.position.inside=c(ifelse(i==3, 0.195, 0.185), 0.203),
#           legend.background=element_blank(),
#           legend.key.height=unit(0.395, "cm"),
#           legend.key.width=unit(0.0, "cm"),
#           legend.text=element_text(size=6),
#           legend.title=element_text(size=8, vjust=1, hjust=1),
#           legend.ticks=element_line(colour="grey10", linewidth=0.25),
#           legend.ticks.length=unit(0.04, "cm"),
#           axis.title=element_blank()) 
# }
# 
# ggarrange(p_ls[[1]], p_ls[[2]], p_ls[[3]], p_ls[[4]], nrow=1, common.legend=F, labels="auto") |> 
#   ggsave("figs/pub/ens_farm-mae+IDW_map_sLonLatD4_all.png", plot=_ , width=12, height=6, dpi=300)
# 
# 
# ggsave("figs/talk/mae+IDW_ensFcst.png", p_ls[[1]], width=3, height=5.5)
# ggsave("figs/talk/mae+IDW_ensBlend.png", p_ls[[2]], width=3, height=5.5)
# ggsave("figs/talk/mae+IDW_ens3D.png", p_ls[[3]], width=3, height=5.5)
# ggsave("figs/talk/mae+IDW_ens2D.png", p_ls[[4]], width=3, height=5.5)
# 




# maps of AUC + IDW ---------------------------------------------------------

library(terra)

mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_rast <- st_read("data/WeStCOMS2_meshFootprint.gpkg") |>
  rast(resolution=500)

farm_AUC.df <- metrics_by_farm |>
  filter(N >= 30) |>
  # filter(grepl("avg|pred", sim)) |>
  filter(sim %in% c("predF1", "predBlend", "sim_09", "sim_17")) |>
  left_join(sim_i) |>
  droplevels() |>
  inner_join(site_i) |>
  drop_na(ROC_AUC) |>
  mutate(lab=lvls_revalue(lab, c("Ens['Fcst']", "Ens['Blend']", "Opt['3D']", "Opt['2D']")))

AUC_map_interp_df <- map(unique(farm_AUC.df$lab), 
                     ~interpIDW(mesh_rast, 
                                farm_AUC.df |>
                                  filter(lab==.x) |>
                                  select(easting, northing, ROC_AUC) |>
                                  drop_na() |>
                                  as.matrix(),
                                radius=1000e3) |>
                       mask(mesh_fp)) |>
  map2_dfr(unique(farm_AUC.df$lab), 
           ~as_tibble(.x) |>
             bind_cols(crds(.x)) |>
             mutate(lab=.y)) |>
  mutate(lab=factor(lab, levels=levels(farm_AUC.df$lab)))

p_AUC <- AUC_map_interp_df |>
  ggplot() + 
  geom_sf(data=mesh_fp, fill="grey90", colour="grey", size=0.1) + 
  geom_raster(aes(x, y, fill=lyr.1)) + 
  stat_contour(aes(x, y, z=lyr.1), breaks=seq(0, 1, by=0.1), 
               colour="grey40", linewidth=0.1) +
  geom_point(data=farm_AUC.df, aes(easting, northing, fill=ROC_AUC), 
             shape=21, size=1, stroke=0.25, colour="grey10") +
  colorspace::scale_fill_binned_diverging(
    name="AUC", palette="Blue-Red 3", rev=T, mid=0.5,
    limits=c(0,1), breaks=seq(0, 1, by=0.1), labels=seq(0, 1, by=0.1),
    l1=20, l2=90, p2=2) +
  scale_y_continuous(limits=c(630000, 955000), oob=scales::oob_keep,
                     breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
  scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
                     limits=c(75000, 235000), oob=scales::oob_keep) +
  facet_grid(.~lab, labeller=label_parsed) +
  theme(axis.title=element_blank(),
        legend.text=element_text(size=6),
        legend.title=element_text(size=8),
        legend.ticks=element_line(colour="grey40", linewidth=0.1),
        legend.key.width=unit(0.2, "cm"),
        legend.key.height=unit(0.6, "cm")) 



farm_RMSE.df <- metrics_by_farm |>
  filter(N >= 30) |>
  # filter(grepl("avg|pred", sim)) |>
  filter(sim %in% c("predF1", "predBlend", "sim_09", "sim_17")) |>
  left_join(sim_i) |>
  droplevels() |>
  inner_join(site_i) |>
  drop_na(rmse) |>
  mutate(lab=lvls_revalue(lab, c("Ens['Fcst']", "Ens['Blend']", "Opt['3D']", "Opt['2D']")))

RMSE_map_interp_df <- map(unique(farm_RMSE.df$lab), 
                     ~interpIDW(mesh_rast, 
                                farm_RMSE.df |>
                                  filter(lab==.x) |>
                                  select(easting, northing, rmse) |>
                                  drop_na() |>
                                  as.matrix(),
                                radius=1000e3) |>
                       mask(mesh_fp)) |>
  map2_dfr(unique(farm_RMSE.df$lab), 
           ~as_tibble(.x) |>
             bind_cols(crds(.x)) |>
             mutate(lab=.y)) |>
  mutate(lab=factor(lab, levels=levels(farm_RMSE.df$lab)))

p_RMSE <- RMSE_map_interp_df |>
  ggplot() + 
  geom_sf(data=mesh_fp, fill="grey90", colour="grey", size=0.1) + 
  geom_raster(aes(x, y, fill=lyr.1)) + 
  stat_contour(aes(x, y, z=lyr.1), breaks=seq(0.1, 0.6, by=0.05), 
               colour="grey40", linewidth=0.1) +
  geom_point(data=farm_RMSE.df, aes(easting, northing, fill=rmse), 
             shape=21, size=1, stroke=0.25, colour="grey10") +
  colorspace::scale_fill_binned_diverging(
    name="RMSE", palette="Tropic", mid=0.35, 
    limits=c(0.1, 0.6), breaks=seq(0.1, 0.6, by=0.05),
    l1=20, l2=90, p2=2) +
  scale_y_continuous(limits=c(630000, 955000), oob=scales::oob_keep,
                     breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
  scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
                     limits=c(75000, 235000), oob=scales::oob_keep) +
  facet_grid(.~lab, labeller=label_parsed) +
  theme(axis.title=element_blank(),
        legend.text=element_text(size=6),
        legend.title=element_text(size=8),
        legend.ticks=element_line(colour="grey40", linewidth=0.1),
        legend.key.width=unit(0.2, "cm"),
        legend.key.height=unit(0.6, "cm")) 



ggarrange(p_AUC, p_RMSE, nrow=2, common.legend=F, labels="auto") |> 
  ggsave("figs/pub/ens_AUC-RMSE_map_sLonLatD4_all_continuous.png", plot=_ , width=10, height=10, dpi=300)




# date performance --------------------------------------------------------

mods <- c("predF1", "predF5", "predBlend")

metric_date_df <- metrics_by_week |>
  filter(date > "2021-06-01") |>
  filter(sim %in% mods) |>
  pivot_longer(all_of(c("ROC_AUC", "PR_AUC", "r", "rmse")), names_to="metric") |>
  mutate(metric=factor(metric, levels=c("ROC_AUC", "PR_AUC", "rsq", "r", "rmse"),
                       labels=c("'AUC'['ROC']", "'AUC'['PR']", "R^2", "rho", "RMSE")))
metric_date_df |>
  ggplot(aes(date, value, colour=sim)) +
  geom_point(size=0.5, shape=1) + 
  geom_line(stat="smooth", method="loess", formula=y~x, se=F, span=0.2) +
  geom_text(data=metric_date_df |> slice_max(date), 
            aes(date, value, label=sim),
            hjust=0, nudge_x=20, vjust=0.5, size=2.5, parse=T) +
  scale_x_date(date_breaks="1 year", date_minor_breaks="3 months", 
               date_labels="%Y", expand=expansion(mult=c(0.05, 0.1))) +
  scale_y_continuous("Among farm value") +
  # scale_colour_manual(values=colorspace::divergingx_hcl("Earth", n=6)[c(1,2,5,6)]) +
  theme_bw() + 
  facet_wrap(~metric, ncol=1, strip.position="right", scales="free_y", labeller=label_parsed) +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=9),
        panel.grid.major.y=element_line(colour="grey85", linewidth=0.4),
        panel.grid.minor.y=element_line(colour="grey90", linewidth=0.2),
        axis.text=element_text(size=7))







# 2D vs 3D IDW r ----------------------------------------------------------

library(terra)
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_rast <- st_read("data/WeStCOMS2_meshFootprint.gpkg") |>
  rast(resolution=500)

farm_r.df <- ensCV_df |>
  group_by(sepaSite) |>
  summarise(across(starts_with("IP"), ~cor(.x, licePerFish_rtrt, method="spearman"))) |>
  inner_join(site_i)

interp_3D <- interpIDW(mesh_rast, 
                        farm_r.df |> 
                          select(easting, northing, IP_sim_avg3D) |> 
                          drop_na() |>
                          as.matrix(),
                        radius=1000e3) |>
  mask(mesh_fp)
interp_2D <- interpIDW(mesh_rast, 
                       farm_r.df |> 
                         select(easting, northing, IP_sim_avg2D) |> 
                         drop_na() |>
                         as.matrix(),
                       radius=1000e3) |>
  mask(mesh_fp)
interp_comp_df <- full_join(
  as_tibble(interp_3D) |> bind_cols(crds(interp_3D)) |> rename(r_3D=lyr.1),
  as_tibble(interp_2D) |> bind_cols(crds(interp_2D)) |> rename(r_2D=lyr.1)
) |>
  mutate(d3m2=r_3D-r_2D)
farm_comp_df <- 
  farm_r.df |>
  mutate(d3m2=IP_sim_avg3D - IP_sim_avg2D) |>
  select(easting, northing, d3m2) 
p <- interp_comp_df |>
  ggplot() + 
  geom_sf(data=mesh_fp, fill="grey90", colour="grey", size=0.1) + 
  geom_raster(aes(x, y, fill=d3m2)) + 
  annotate("text", x=224000, y=862500, label="r :", size=2.5) +
  # scale_fill_gradient2() +
  colorspace::scale_fill_binned_diverging(
    palette="Blue-Red 3", rev=T, l1=20, l2=90, p2=2,
    limits=c(-max(abs(interp_comp_df$d3m2)), max(abs(interp_comp_df$d3m2))), 
    breaks=seq(-2, 2, by=0.5)) +
  # colorspace::scale_fill_continuous_diverging(
  #   palette="Blue-Red 3", rev=T, l1=20, l2=90, p2=2,
  #   limits=c(-max(abs(interp_comp_df$d3m2)), max(abs(interp_comp_df$d3m2))), 
  #   breaks=c(-1, 0, 1)) +
  scale_y_continuous(limits=c(620000, 975000), oob=scales::oob_keep,
                     breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
  scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
                     limits=c(63000, 243000), oob=scales::oob_keep) +
  theme(legend.position="inside",
        legend.position.inside=c(0.9, 0.57),
        legend.text=element_text(size=6),
        legend.title=element_text(size=7, hjust=0.1),
        legend.background=element_blank(),
        legend.key.height=unit(0.3, "cm"),
        legend.key.width=unit(0.15, "cm"),
        axis.title=element_blank()) 
ggsave("figs/pub/ens_farm-r_map_3D-2D.png", plot=p, width=3, height=5.25, dpi=300)




# 2D vs 3D IDW rmse ----------------------------------------------------------

library(terra)
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_rast <- st_read("data/WeStCOMS2_meshFootprint.gpkg") |>
  rast(resolution=500)

farm_rmse.df <- ensCV_df |>
  group_by(sepaSite) |>
  summarise(across(starts_with("IP"), ~yardstick::rmse_vec(.x, truth=licePerFish_rtrt))) |>
  inner_join(site_i)

interp_3D <- interpIDW(mesh_rast, 
                       farm_rmse.df |> 
                         select(easting, northing, IP_sim_avg3D) |> 
                         drop_na() |>
                         as.matrix(),
                       radius=1000e3) |>
  mask(mesh_fp)
interp_2D <- interpIDW(mesh_rast, 
                       farm_rmse.df |> 
                         select(easting, northing, IP_sim_avg2D) |> 
                         drop_na() |>
                         as.matrix(),
                       radius=1000e3) |>
  mask(mesh_fp)
interp_comp_df <- full_join(
  as_tibble(interp_3D) |> bind_cols(crds(interp_3D)) |> rename(rmse_3D=lyr.1),
  as_tibble(interp_2D) |> bind_cols(crds(interp_2D)) |> rename(rmse_2D=lyr.1)
) |>
  mutate(d3m2=rmse_3D-rmse_2D)
farm_comp_df <- 
  farm_rmse.df |>
  mutate(d3m2=IP_sim_avg3D - IP_sim_avg2D) |>
  select(easting, northing, d3m2) 
p <- interp_comp_df |>
  ggplot() + 
  geom_sf(data=mesh_fp, fill="grey90", colour="grey", size=0.1) + 
  geom_raster(aes(x, y, fill=d3m2)) + 
  annotate("text", x=224000, y=862500, label="rmse :", size=2.5) +
  # scale_fill_gradient2() +
  colorspace::scale_fill_binned_diverging(
    palette="Blue-Red 3", rev=F, l1=20, l2=90, p2=2,
    limits=c(-max(abs(interp_comp_df$d3m2)), max(abs(interp_comp_df$d3m2))), 
    breaks=seq(-0.25, 0.25, by=0.05)) +
  # colorspace::scale_fill_continuous_diverging(
  #   palette="Blue-Red 3", rev=T, l1=20, l2=90, p2=2,
  #   limits=c(-max(abs(interp_comp_df$d3m2)), max(abs(interp_comp_df$d3m2))), 
  #   breaks=c(-1, 0, 1)) +
  scale_y_continuous(limits=c(620000, 975000), oob=scales::oob_keep,
                     breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
  scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
                     limits=c(63000, 243000), oob=scales::oob_keep) +
  theme(legend.position="inside",
        legend.position.inside=c(0.865, 0.57),
        legend.text=element_text(size=6),
        legend.title=element_text(size=7, hjust=0.1),
        legend.background=element_blank(),
        legend.key.height=unit(0.375, "cm"),
        legend.key.width=unit(0.15, "cm"),
        axis.title=element_blank()) 
ggsave("figs/pub/ens_farm-rmse_map_3D-2D.png", plot=p, width=3, height=5.25, dpi=300)



# 2D vs 3D IDW AUC ----------------------------------------------------------

library(terra)
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_rast <- st_read("data/WeStCOMS2_meshFootprint.gpkg") |>
  rast(resolution=500)

farm_AUC.df <- ensCV_df |>
  group_by(sepaSite) |>
  summarise(across(starts_with("IP"), 
                   ~yardstick::roc_auc_vec(.x, truth=lice_g05, event_level="second"))) |>
  inner_join(site_i)

interp_3D <- interpIDW(mesh_rast, 
                       farm_AUC.df |> 
                         select(easting, northing, IP_sim_avg3D) |> 
                         drop_na() |>
                         as.matrix(),
                       radius=1000e3) |>
  mask(mesh_fp)
interp_2D <- interpIDW(mesh_rast, 
                       farm_AUC.df |> 
                         select(easting, northing, IP_sim_avg2D) |> 
                         drop_na() |>
                         as.matrix(),
                       radius=1000e3) |>
  mask(mesh_fp)
interp_comp_df <- full_join(
  as_tibble(interp_3D) |> bind_cols(crds(interp_3D)) |> rename(AUC_3D=lyr.1),
  as_tibble(interp_2D) |> bind_cols(crds(interp_2D)) |> rename(AUC_2D=lyr.1)
) |>
  mutate(d3m2=AUC_3D-AUC_2D)
farm_comp_df <- 
  farm_AUC.df |>
  mutate(d3m2=IP_sim_avg3D - IP_sim_avg2D) |>
  select(easting, northing, d3m2) 
p <- interp_comp_df |>
  ggplot() + 
  geom_sf(data=mesh_fp, fill="grey90", colour="grey", size=0.1) + 
  geom_raster(aes(x, y, fill=d3m2)) + 
  annotate("text", x=224000, y=862500, label="AUC :", size=2.5) +
  # scale_fill_gradient2() +
  colorspace::scale_fill_binned_diverging(
    palette="Blue-Red 3", rev=T, l1=20, l2=90, p2=2, 
    limits=round(c(-max(abs(interp_comp_df$d3m2)), max(abs(interp_comp_df$d3m2))), 2), 
    breaks=seq(-1, 1, by=0.1)) +
  # colorspace::scale_fill_continuous_diverging(
  #   palette="Blue-Red 3", rev=T, l1=20, l2=90, p2=2,
  #   limits=c(-max(abs(interp_comp_df$d3m2)), max(abs(interp_comp_df$d3m2))), 
  #   breaks=c(-1, 0, 1)) +
  scale_y_continuous(limits=c(620000, 975000), oob=scales::oob_keep,
                     breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
  scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
                     limits=c(63000, 243000), oob=scales::oob_keep) +
  theme(legend.position="inside",
        legend.position.inside=c(0.865, 0.57),
        legend.text=element_text(size=6),
        legend.title=element_text(size=7, hjust=0.1),
        legend.background=element_blank(),
        legend.key.height=unit(0.375, "cm"),
        legend.key.width=unit(0.15, "cm"),
        axis.title=element_blank()) 
ggsave("figs/pub/ens_farm-AUC_map_3D-2D.png", plot=p, width=3, height=5.25, dpi=300)










# 2D vs EnsBlend IDW r ----------------------------------------------------------

library(terra)
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_rast <- st_read("data/WeStCOMS2_meshFootprint.gpkg") |>
  rast(resolution=500)

farm_r.df <- ensCV_df |>
  group_by(sepaSite) |>
  summarise(across(starts_with("IP"), ~cor(.x, licePerFish_rtrt, method="spearman"))) |>
  inner_join(site_i)

interp_3D <- interpIDW(mesh_rast, 
                       farm_r.df |> 
                         select(easting, northing, IP_predBlend) |> 
                         drop_na() |>
                         as.matrix(),
                       radius=1000e3) |>
  mask(mesh_fp)
interp_2D <- interpIDW(mesh_rast, 
                       farm_r.df |> 
                         select(easting, northing, IP_sim_avg2D) |> 
                         drop_na() |>
                         as.matrix(),
                       radius=1000e3) |>
  mask(mesh_fp)
interp_comp_df <- full_join(
  as_tibble(interp_3D) |> bind_cols(crds(interp_3D)) |> rename(r_3D=lyr.1),
  as_tibble(interp_2D) |> bind_cols(crds(interp_2D)) |> rename(r_2D=lyr.1)
) |>
  mutate(d3m2=r_3D-r_2D)
farm_comp_df <- 
  farm_r.df |>
  mutate(d3m2=IP_predBlend - IP_sim_avg2D) |>
  select(easting, northing, d3m2) 
p <- interp_comp_df |>
  ggplot() + 
  geom_sf(data=mesh_fp, fill="grey90", colour="grey", size=0.1) + 
  geom_raster(aes(x, y, fill=d3m2)) + 
  annotate("text", x=224000, y=862500, label="r :", size=2.5) +
  # scale_fill_gradient2() +
  colorspace::scale_fill_binned_diverging(
    palette="Blue-Red 3", rev=T, l1=20, l2=90, p2=2,
    limits=c(-max(abs(interp_comp_df$d3m2)), max(abs(interp_comp_df$d3m2))), 
    breaks=seq(-2, 2, by=0.5)) +
  # colorspace::scale_fill_continuous_diverging(
  #   palette="Blue-Red 3", rev=T, l1=20, l2=90, p2=2,
  #   limits=c(-max(abs(interp_comp_df$d3m2)), max(abs(interp_comp_df$d3m2))), 
  #   breaks=c(-1, 0, 1)) +
  scale_y_continuous(limits=c(620000, 975000), oob=scales::oob_keep,
                     breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
  scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
                     limits=c(63000, 243000), oob=scales::oob_keep) +
  theme(legend.position="inside",
        legend.position.inside=c(0.9, 0.57),
        legend.text=element_text(size=6),
        legend.title=element_text(size=7, hjust=0.1),
        legend.background=element_blank(),
        legend.key.height=unit(0.3, "cm"),
        legend.key.width=unit(0.15, "cm"),
        axis.title=element_blank()) 
ggsave("figs/pub/ens_farm-r_map_ensBlend-2D.png", plot=p, width=3, height=5.25, dpi=300)




# 2D vs EnsBlend IDW rmse ----------------------------------------------------------

library(terra)
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_rast <- st_read("data/WeStCOMS2_meshFootprint.gpkg") |>
  rast(resolution=500)

farm_rmse.df <- ensCV_df |>
  group_by(sepaSite) |>
  summarise(across(starts_with("IP"), ~yardstick::rmse_vec(.x, truth=licePerFish_rtrt))) |>
  inner_join(site_i)

interp_3D <- interpIDW(mesh_rast, 
                       farm_rmse.df |> 
                         select(easting, northing, IP_predBlend) |> 
                         drop_na() |>
                         as.matrix(),
                       radius=1000e3) |>
  mask(mesh_fp)
interp_2D <- interpIDW(mesh_rast, 
                       farm_rmse.df |> 
                         select(easting, northing, IP_sim_avg2D) |> 
                         drop_na() |>
                         as.matrix(),
                       radius=1000e3) |>
  mask(mesh_fp)
interp_comp_df <- full_join(
  as_tibble(interp_3D) |> bind_cols(crds(interp_3D)) |> rename(rmse_3D=lyr.1),
  as_tibble(interp_2D) |> bind_cols(crds(interp_2D)) |> rename(rmse_2D=lyr.1)
) |>
  mutate(d3m2=rmse_3D-rmse_2D)
farm_comp_df <- 
  farm_rmse.df |>
  mutate(d3m2=IP_predBlend - IP_sim_avg2D) |>
  select(easting, northing, d3m2) 
p <- interp_comp_df |>
  ggplot() + 
  geom_sf(data=mesh_fp, fill="grey90", colour="grey", size=0.1) + 
  geom_raster(aes(x, y, fill=d3m2)) + 
  annotate("text", x=224000, y=862500, label="rmse :", size=2.5) +
  # scale_fill_gradient2() +
  # colorspace::scale_fill_binned_diverging(
  #   palette="Blue-Red 3", rev=F, l1=20, l2=90, p2=2,
  #   limits=round(c(-max(abs(interp_comp_df$d3m2)), max(abs(interp_comp_df$d3m2))), 3), 
  #   breaks=seq(-0.3, 0.3, by=0.05)) +
  colorspace::scale_fill_continuous_diverging(
    palette="Blue-Red 3", rev=F, l1=20, l2=90, p2=2,
    limits=c(-max(abs(interp_comp_df$d3m2)), max(abs(interp_comp_df$d3m2)))) +
  scale_y_continuous(limits=c(620000, 975000), oob=scales::oob_keep,
                     breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
  scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
                     limits=c(63000, 243000), oob=scales::oob_keep) +
  theme(legend.position="inside",
        legend.position.inside=c(0.865, 0.57),
        legend.text=element_text(size=6),
        legend.title=element_text(size=7, hjust=0.1),
        legend.background=element_blank(),
        legend.key.height=unit(0.375, "cm"),
        legend.key.width=unit(0.15, "cm"),
        axis.title=element_blank()) 
ggsave("figs/pub/ens_farm-rmse_map_ensBlend-2D.png", plot=p, width=3, height=5.25, dpi=300)



# 2D vs EnsBlend IDW AUC ----------------------------------------------------------

library(terra)
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_rast <- st_read("data/WeStCOMS2_meshFootprint.gpkg") |>
  rast(resolution=500)

farm_AUC.df <- ensCV_df |>
  group_by(sepaSite) |>
  summarise(across(starts_with("IP"), 
                   ~yardstick::roc_auc_vec(.x, truth=lice_g05, event_level="second"))) |>
  inner_join(site_i)

interp_3D <- interpIDW(mesh_rast, 
                       farm_AUC.df |> 
                         select(easting, northing, IP_predBlend) |> 
                         drop_na() |>
                         as.matrix(),
                       radius=1000e3) |>
  mask(mesh_fp)
interp_2D <- interpIDW(mesh_rast, 
                       farm_AUC.df |> 
                         select(easting, northing, IP_sim_avg2D) |> 
                         drop_na() |>
                         as.matrix(),
                       radius=1000e3) |>
  mask(mesh_fp)
interp_comp_df <- full_join(
  as_tibble(interp_3D) |> bind_cols(crds(interp_3D)) |> rename(AUC_3D=lyr.1),
  as_tibble(interp_2D) |> bind_cols(crds(interp_2D)) |> rename(AUC_2D=lyr.1)
) |>
  mutate(d3m2=AUC_3D-AUC_2D)
farm_comp_df <- 
  farm_AUC.df |>
  mutate(d3m2=IP_predBlend - IP_sim_avg2D) |>
  select(easting, northing, d3m2) 
p <- interp_comp_df |>
  ggplot() + 
  geom_sf(data=mesh_fp, fill="grey90", colour="grey", size=0.1) + 
  geom_raster(aes(x, y, fill=d3m2)) + 
  annotate("text", x=224000, y=862500, label="AUC :", size=2.5) +
  # scale_fill_gradient2() +
  colorspace::scale_fill_binned_diverging(
    palette="Blue-Red 3", rev=T, l1=20, l2=90, p2=2, 
    limits=round(c(-max(abs(interp_comp_df$d3m2)), max(abs(interp_comp_df$d3m2))), 2), 
    breaks=seq(-1, 1, by=0.1)) +
  # colorspace::scale_fill_continuous_diverging(
  #   palette="Blue-Red 3", rev=T, l1=20, l2=90, p2=2,
  #   limits=c(-max(abs(interp_comp_df$d3m2)), max(abs(interp_comp_df$d3m2))), 
  #   breaks=c(-1, 0, 1)) +
  scale_y_continuous(limits=c(620000, 975000), oob=scales::oob_keep,
                     breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
  scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
                     limits=c(63000, 243000), oob=scales::oob_keep) +
  theme(legend.position="inside",
        legend.position.inside=c(0.865, 0.57),
        legend.text=element_text(size=6),
        legend.title=element_text(size=7, hjust=0.1),
        legend.background=element_blank(),
        legend.key.height=unit(0.375, "cm"),
        legend.key.width=unit(0.15, "cm"),
        axis.title=element_blank()) 
ggsave("figs/pub/ens_farm-AUC_map_ensBlend-2D.png", plot=p, width=3, height=5.25, dpi=300)












# 2D vs EnsFcst IDW r ----------------------------------------------------------

library(terra)
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_rast <- st_read("data/WeStCOMS2_meshFootprint.gpkg") |>
  rast(resolution=500)

farm_r.df <- ensCV_df |>
  group_by(sepaSite) |>
  summarise(across(starts_with("IP"), ~cor(.x, licePerFish_rtrt, method="spearman"))) |>
  inner_join(site_i)

interp_3D <- interpIDW(mesh_rast, 
                       farm_r.df |> 
                         select(easting, northing, IP_predF1) |> 
                         drop_na() |>
                         as.matrix(),
                       radius=1000e3) |>
  mask(mesh_fp)
interp_2D <- interpIDW(mesh_rast, 
                       farm_r.df |> 
                         select(easting, northing, IP_sim_avg2D) |> 
                         drop_na() |>
                         as.matrix(),
                       radius=1000e3) |>
  mask(mesh_fp)
interp_comp_df <- full_join(
  as_tibble(interp_3D) |> bind_cols(crds(interp_3D)) |> rename(r_3D=lyr.1),
  as_tibble(interp_2D) |> bind_cols(crds(interp_2D)) |> rename(r_2D=lyr.1)
) |>
  mutate(d3m2=r_3D-r_2D)
farm_comp_df <- 
  farm_r.df |>
  mutate(d3m2=IP_predF1 - IP_sim_avg2D) |>
  select(easting, northing, d3m2) 
p <- interp_comp_df |>
  ggplot() + 
  geom_sf(data=mesh_fp, fill="grey90", colour="grey", size=0.1) + 
  geom_raster(aes(x, y, fill=d3m2)) + 
  annotate("text", x=224000, y=862500, label="r :", size=2.5) +
  # scale_fill_gradient2() +
  colorspace::scale_fill_binned_diverging(
    palette="Blue-Red 3", rev=T, l1=20, l2=90, p2=2,
    limits=c(-max(abs(interp_comp_df$d3m2)), max(abs(interp_comp_df$d3m2))), 
    breaks=seq(-2, 2, by=0.5)) +
  # colorspace::scale_fill_continuous_diverging(
  #   palette="Blue-Red 3", rev=T, l1=20, l2=90, p2=2,
  #   limits=c(-max(abs(interp_comp_df$d3m2)), max(abs(interp_comp_df$d3m2))), 
  #   breaks=c(-1, 0, 1)) +
  scale_y_continuous(limits=c(620000, 975000), oob=scales::oob_keep,
                     breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
  scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
                     limits=c(63000, 243000), oob=scales::oob_keep) +
  theme(legend.position="inside",
        legend.position.inside=c(0.9, 0.57),
        legend.text=element_text(size=6),
        legend.title=element_text(size=7, hjust=0.1),
        legend.background=element_blank(),
        legend.key.height=unit(0.3, "cm"),
        legend.key.width=unit(0.15, "cm"),
        axis.title=element_blank()) 
ggsave("figs/pub/ens_farm-r_map_ensFcst-2D.png", plot=p, width=3, height=5.25, dpi=300)




# 2D vs EnsFcst IDW rmse ----------------------------------------------------------

library(terra)
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_rast <- st_read("data/WeStCOMS2_meshFootprint.gpkg") |>
  rast(resolution=500)

farm_rmse.df <- ensCV_df |>
  group_by(sepaSite) |>
  summarise(across(starts_with("IP"), ~yardstick::rmse_vec(.x, truth=licePerFish_rtrt))) |>
  inner_join(site_i)

interp_3D <- interpIDW(mesh_rast, 
                       farm_rmse.df |> 
                         select(easting, northing, IP_predF1) |> 
                         drop_na() |>
                         as.matrix(),
                       radius=1000e3) |>
  mask(mesh_fp)
interp_2D <- interpIDW(mesh_rast, 
                       farm_rmse.df |> 
                         select(easting, northing, IP_sim_avg2D) |> 
                         drop_na() |>
                         as.matrix(),
                       radius=1000e3) |>
  mask(mesh_fp)
interp_comp_df <- full_join(
  as_tibble(interp_3D) |> bind_cols(crds(interp_3D)) |> rename(rmse_3D=lyr.1),
  as_tibble(interp_2D) |> bind_cols(crds(interp_2D)) |> rename(rmse_2D=lyr.1)
) |>
  mutate(d3m2=rmse_3D-rmse_2D)
farm_comp_df <- 
  farm_rmse.df |>
  mutate(d3m2=IP_predF1 - IP_sim_avg2D) |>
  select(easting, northing, d3m2) 
p <- interp_comp_df |>
  ggplot() + 
  geom_sf(data=mesh_fp, fill="grey90", colour="grey", size=0.1) + 
  geom_raster(aes(x, y, fill=d3m2)) + 
  annotate("text", x=224000, y=862500, label="rmse :", size=2.5) +
  # scale_fill_gradient2() +
  colorspace::scale_fill_binned_diverging(
    palette="Blue-Red 3", rev=F, l1=20, l2=90, p2=2,
    limits=c(-max(abs(interp_comp_df$d3m2)), max(abs(interp_comp_df$d3m2))), 
    breaks=seq(-0.3, 0.3, by=0.05)) +
  # colorspace::scale_fill_continuous_diverging(
  #   palette="Blue-Red 3", rev=T, l1=20, l2=90, p2=2,
  #   limits=c(-max(abs(interp_comp_df$d3m2)), max(abs(interp_comp_df$d3m2))), 
  #   breaks=c(-1, 0, 1)) +
  scale_y_continuous(limits=c(620000, 975000), oob=scales::oob_keep,
                     breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
  scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
                     limits=c(63000, 243000), oob=scales::oob_keep) +
  theme(legend.position="inside",
        legend.position.inside=c(0.865, 0.57),
        legend.text=element_text(size=6),
        legend.title=element_text(size=7, hjust=0.1),
        legend.background=element_blank(),
        legend.key.height=unit(0.375, "cm"),
        legend.key.width=unit(0.15, "cm"),
        axis.title=element_blank()) 
ggsave("figs/pub/ens_farm-rmse_map_ensFcst-2D.png", plot=p, width=3, height=5.25, dpi=300)



# 2D vs EnsFcst IDW AUC ----------------------------------------------------------

library(terra)
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_rast <- st_read("data/WeStCOMS2_meshFootprint.gpkg") |>
  rast(resolution=500)

farm_AUC.df <- ensCV_df |>
  group_by(sepaSite) |>
  summarise(across(c(starts_with("IP"), starts_with("pr")), 
                   ~yardstick::roc_auc_vec(.x, truth=lice_g05, event_level="second"))) |>
  inner_join(site_i)

interp_3D <- interpIDW(mesh_rast, 
                       farm_AUC.df |> 
                         select(easting, northing, pr_predF1) |> 
                         drop_na() |>
                         as.matrix(),
                       radius=1000e3) |>
  mask(mesh_fp)
interp_2D <- interpIDW(mesh_rast, 
                       farm_AUC.df |> 
                         select(easting, northing, IP_sim_avg2D) |> 
                         drop_na() |>
                         as.matrix(),
                       radius=1000e3) |>
  mask(mesh_fp)
interp_comp_df <- full_join(
  as_tibble(interp_3D) |> bind_cols(crds(interp_3D)) |> rename(AUC_3D=lyr.1),
  as_tibble(interp_2D) |> bind_cols(crds(interp_2D)) |> rename(AUC_2D=lyr.1)
) |>
  mutate(d3m2=AUC_3D-AUC_2D)
farm_comp_df <- 
  farm_AUC.df |>
  mutate(d3m2=pr_predF1 - IP_sim_avg2D) |>
  select(easting, northing, d3m2) 
p <- interp_comp_df |>
  ggplot() + 
  geom_sf(data=mesh_fp, fill="grey90", colour="grey", size=0.1) + 
  geom_raster(aes(x, y, fill=d3m2)) + 
  annotate("text", x=224000, y=862500, label="AUC :", size=2.5) +
  # scale_fill_gradient2() +
  colorspace::scale_fill_binned_diverging(
    palette="Blue-Red 3", rev=T, l1=20, l2=90, p2=2, 
    limits=round(c(-max(abs(interp_comp_df$d3m2)), max(abs(interp_comp_df$d3m2))), 2), 
    breaks=seq(-1, 1, by=0.1)) +
  # colorspace::scale_fill_continuous_diverging(
  #   palette="Blue-Red 3", rev=T, l1=20, l2=90, p2=2,
  #   limits=c(-max(abs(interp_comp_df$d3m2)), max(abs(interp_comp_df$d3m2))), 
  #   breaks=c(-1, 0, 1)) +
  scale_y_continuous(limits=c(620000, 975000), oob=scales::oob_keep,
                     breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
  scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
                     limits=c(63000, 243000), oob=scales::oob_keep) +
  theme(legend.position="inside",
        legend.position.inside=c(0.865, 0.57),
        legend.text=element_text(size=6),
        legend.title=element_text(size=7, hjust=0.1),
        legend.background=element_blank(),
        legend.key.height=unit(0.375, "cm"),
        legend.key.width=unit(0.15, "cm"),
        axis.title=element_blank()) 
ggsave("figs/pub/ens_farm-AUC_map_ensFcst-2D.png", plot=p, width=3, height=5.25, dpi=300)

















# bad performers ----------------------------------------------------------

metrics_by_farm |> 
  filter(N >= 30) |>
  filter(grepl("avg|pred", sim)) |>
  arrange(desc(rmse))
metrics_by_farm |>
  ggplot(aes(ROC_AUC, sepaSite)) +
  geom_boxplot()





# variability among simulations -------------------------------------------

c_df <- readRDS("out/sim_2019-2023/processed/connectivity_wk.rds")
p <- vector("list", 3)
p[[1]] <- c_df |> 
  group_by(sepaSite, week) |> 
  summarise(influx_rng=diff(range(influx_m2, na.rm=T))) |> 
  group_by(sepaSite) |> 
  summarise(md_rng=median(influx_rng, na.rm=T)) |> 
  left_join(site_i) |> 
  ggplot() + 
  geom_sf(data=mesh_fp) + 
  geom_point(aes(easting, northing, fill=md_rng), shape=21, size=2) +
  scale_fill_viridis_c(option="turbo", begin=0.1) + 
  scale_y_continuous(limits=c(630000, 955000), oob=scales::oob_keep,
                     breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
  scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
                     limits=c(75000, 235000), oob=scales::oob_keep)
p[[2]] <- c_df |> 
  group_by(sepaSite, week) |> 
  summarise(influx_CV=sd(influx_m2, na.rm=T)/mean(influx_m2, na.rm=T)) |> 
  group_by(sepaSite) |> 
  summarise(mn_CV=mean(influx_CV, na.rm=T)) |> 
  left_join(site_i) |> 
  ggplot() + 
  geom_sf(data=mesh_fp) + 
  geom_point(aes(easting, northing, fill=mn_CV), shape=21, size=2) +
  scale_fill_viridis_c(option="turbo", begin=0.1) + 
  scale_y_continuous(limits=c(630000, 955000), oob=scales::oob_keep,
                     breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
  scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
                     limits=c(75000, 235000), oob=scales::oob_keep)
p[[3]] <- c_df |> 
  group_by(sepaSite, week) |> 
  summarise(influx_sd=sd(influx_m2, na.rm=T)) |> 
  group_by(sepaSite) |> 
  summarise(mn_sd=mean(influx_sd^0.25, na.rm=T)) |> 
  left_join(site_i) |> 
  ggplot() + 
  geom_sf(data=mesh_fp) + 
  geom_point(aes(easting, northing, fill=mn_sd), shape=21, size=2) +
  scale_fill_viridis_c(option="turbo", begin=0.1) + 
  scale_y_continuous(limits=c(630000, 955000), oob=scales::oob_keep,
                     breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
  scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
                     limits=c(75000, 235000), oob=scales::oob_keep)
ggpubr::ggarrange(plotlist=p, nrow=1, common.legend=F)

c_df |> 
  group_by(sepaSite, week) |> 
  summarise(influx_sd=sd(influx_m2^0.25, na.rm=T)) |>
  ggplot(aes(week, influx_sd)) + geom_line() + facet_wrap(~sepaSite)
c_df |> 
  group_by(sepaSite, week) |> 
  summarise(influx_sd=sd(influx_m2^0.25, na.rm=T)) |>
  ggplot(aes(week, influx_sd, group=sepaSite)) + geom_line(alpha=0.25)
c_df |> 
  mutate(influx_m2=influx_m2^0.25) |>
  ggplot(aes(week, influx_m2, colour=sim)) + geom_line() + facet_wrap(~sepaSite)




# vertical distributions --------------------------------------------------

z_ens <- map_dfr(c("WeStCOMS", "Linnhe", "Skye"),
                 ~readRDS(glue("{out_dir}/processed/summary_daily_z_{.x}.rds")) |>
                   mutate(region=.x)) |>
  mutate(region=factor(region, levels=c("WeStCOMS", "Linnhe", "Skye"),
                       labels=c("Full domain", "Loch Linnhe", "Skye")))
gc()

z_ens |>
  ggplot(aes(day, prop, fill=z, colour=z, group=z)) +
  geom_area(outline.type="upper", linewidth=0.2) +
  scale_y_continuous("Ensemble proportion of copepodids (daily)") +
  scale_x_date(date_breaks="1 month", date_labels="%b") +
  scale_fill_viridis_b("Depth bin (m)", direction=-1,
                       breaks=c(1, 2, 5, 10, 15, 20, 25, 30)-0.01,
                       labels=c(1, 2, 5, 10, 15, 20, 25, 30)) +
  scale_colour_viridis_b("Depth bin (m)", direction=-1,
                         breaks=c(1, 2, 5, 10, 15, 20, 25, 30)-0.01,
                         labels=c(1, 2, 5, 10, 15, 20, 25, 30)) +
# scale_fill_viridis_b("Depth bin (m)", direction=-1,
#                      breaks=seq(0, 30, by=5)[-1]-0.01,
#                      labels=seq(0, 30, by=5)[-1]) +
# scale_colour_viridis_b("Depth bin (m)", direction=-1,
#                        breaks=seq(0, 30, by=5)[-1]-0.01,
#                        labels=seq(0, 30, by=5)[-1]) +
  # scale_fill_viridis_c(direction=-1) +
  # scale_colour_viridis_c(direction=-1) +
  facet_grid(region~.) +
  theme(axis.title.x=element_blank(),
        panel.grid.major.y=element_line(colour="grey90", linewidth=0.2),
        legend.position="bottom", 
        legend.key.height=unit(0.2, "cm"), 
        legend.key.width=unit(1.5, "cm"))
ggsave("figs/pub/ens_z_distribution.png", width=4.5, height=8)


z_ens |>
  filter(day > "2023-01-14") |>
  ggplot(aes(day, mean_ens_sd/N, fill=z, colour=z, group=z)) +
  geom_area(outline.type="upper", linewidth=0.2) +
  scale_y_continuous("Ensemble proportion of copepodids (daily)") +
  scale_x_date(date_breaks="1 month", date_labels="%b") +
  scale_fill_viridis_b("Depth bin (m)", direction=-1,
                       breaks=c(1, 2, 5, 10, 15, 20, 25, 30)-0.01,
                       labels=c(1, 2, 5, 10, 15, 20, 25, 30)) +
  scale_colour_viridis_b("Depth bin (m)", direction=-1,
                         breaks=c(1, 2, 5, 10, 15, 20, 25, 30)-0.01,
                         labels=c(1, 2, 5, 10, 15, 20, 25, 30)) +
  facet_grid(region~.) +
  theme(axis.title.x=element_blank(),
        panel.grid.major.y=element_line(colour="grey90", linewidth=0.2),
        legend.position="bottom", 
        legend.key.height=unit(0.2, "cm"), 
        legend.key.width=unit(1.5, "cm"))










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
