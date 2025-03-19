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
    tibble(sim=c("predFcst", "predBlend", 
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
  select(rowNum, sepaSite, CV_k, year, date, licePerFish_rtrt, lice_g05) |>
  left_join(read_csv("out/candidates/CV_candidate_predictions.csv")) |>
  left_join(read_csv("out/ensembles/CV_avg_predictions.csv")) |>
  left_join(read_csv("out/ensembles/CV_ensBlend_predictions.csv") |>
              select(rowNum, IP_sLonLatD3_n20) |> rename(IP_predBlend=IP_sLonLatD3_n20)) |>
  left_join(read_csv("out/ensembles/CV_ensFc-1_rmse.csv") |>
              select(rowNum, .pred) |> rename(IP_predFwk_RMSE=.pred)) |>
  left_join(read_csv("out/ensembles/CV_ensFc-1_roc_auc.csv") |>
              select(rowNum, .pred_TRUE) |> rename(IP_predFwk_ROC=.pred_TRUE)) |>
  left_join(read_csv("out/ensembles/CV_ensFc-1_average_precision.csv") |>
              select(rowNum, .pred_TRUE) |> rename(IP_predFwk_PR=.pred_TRUE)) |>
  left_join(read_csv("out/ensembles/CV_ensFc-CV-farms_rmse.csv") |>
              select(rowNum, .pred) |> rename(IP_predFf_RMSE=.pred)) |>
  left_join(read_csv("out/ensembles/CV_ensFc-CV-farms_roc_auc.csv") |>
              select(rowNum, .pred_TRUE) |> rename(IP_predFf_ROC=.pred_TRUE)) |>
  left_join(read_csv("out/ensembles/CV_ensFc-CV-farms_average_precision.csv") |>
              select(rowNum, .pred_TRUE) |> rename(IP_predFf_PR=.pred_TRUE))

folds <- unique(ensFull_df$CV_k)
ensNull_time <- ensNull_farm <- vector("list", length(folds))
for(k in seq_along(folds)) {
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
    mutate(CV_k=folds[k])
}
ensCV_df <- ensCV_df |>
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
              N=n(),
              prop_g05=mean(lice_g05=="TRUE"),
              prop_0=mean(licePerFish_rtrt==0)) |>
    ungroup()
ensF_metrics_by_farm <- metrics_by_farm |>
  filter(grepl("predF", sim)) |>
  group_by(sepaSite) |>
  summarise(rmse=sum(rmse * (sim=="predFf_RMSE")),
            sim="predFcst", 
            across(any_of(c("N", "prop_g05", "prop_0", "minPRAUC")), first))
metrics_by_farm <- metrics_by_farm |>
  filter(!grepl("predF", sim)) |>
  bind_rows(ensF_metrics_by_farm)

# Mean among site
metrics_by_week <- ensCV_df |>
    pivot_longer(starts_with("IP_"), names_to="sim") |>
    mutate(sim=str_remove(sim, "IP_")) |>
    group_by(date, sim) |>
    summarise(rmse=rmse_vec(value, truth=licePerFish_rtrt),
              N=n(),
              prop_g05=mean(lice_g05=="TRUE"),
              prop_0=mean(licePerFish_rtrt==0)) |>
    ungroup()
ensF_metrics_by_week <- metrics_by_week |>
  filter(grepl("predF", sim)) |>
  group_by(date) |>
  summarise(rmse=sum(rmse * (sim=="predFwk_RMSE")),
            sim="predFcst", 
            across(any_of(c("N", "prop_g05", "prop_0", "minPRAUC")), first))
metrics_by_week <- metrics_by_week |>
  filter(!grepl("predF", sim)) |>
  bind_rows(ensF_metrics_by_week)

# Medians
metrics_by_farm_md <- metrics_by_farm |>
  filter(N >= 30) |>
  group_by(sim) |>
  summarise(rmse=median(rmse, na.rm=T),
            N=mean(N, na.rm=T),
            prop_g05=mean(prop_g05),
            prop_0=mean(prop_0, na.rm=T)) |>
  ungroup()
metrics_by_week_md <- metrics_by_week |>
  filter(N >= 30) |>
  group_by(sim) |>
  summarise(rmse=median(rmse, na.rm=T),
            N=mean(N, na.rm=T),
            prop_g05=mean(prop_g05),
            prop_0=mean(prop_0, na.rm=T)) |>
  ungroup()
metrics_by_farm_mn <- metrics_by_farm |>
  filter(N >= 30) |>
  group_by(sim) |>
  summarise(rmse=mean(rmse, na.rm=T),
            N=mean(N, na.rm=T),
            prop_g05=mean(prop_g05),
            prop_0=mean(prop_0, na.rm=T)) |>
  ungroup()
metrics_by_week_mn <- metrics_by_week |>
  filter(N >= 30) |>
  group_by(sim) |>
  summarise(rmse=mean(rmse, na.rm=T),
            N=mean(N, na.rm=T),
            prop_g05=mean(prop_g05),
            prop_0=mean(prop_0, na.rm=T)) |>
  ungroup()

plot_metric_ordered <- function(df, m, highlight="sim_04") {
  df_highlight <- df |>
    filter(sim %in% highlight)
  df |>
    arrange({{m}}) |>
    mutate(sim=factor(sim, levels=unique(sim))) |>
    ggplot(aes({{m}}, sim)) + 
    geom_vline(data=df_highlight, aes(xintercept={{m}}, colour=sim), linetype=2) +
    geom_hline(data=df_highlight, aes(yintercept=sim, colour=sim), linetype=2) +
    geom_point() + 
    theme(legend.position="none")
}

library(cowplot)
highlight_sims <- c("sim_04", #paste0("sLonLatD3_n", c(5, 10, 20)), 
                    "predBlend", "predFcst",
                    paste0("predFwk_", c("RMSE", "ROC", "PR")),
                    paste0("predFf_", c("RMSE", "ROC", "PR")))
map(list(metrics_by_farm_md, metrics_by_week_md),
    ~plot_metric_ordered(.x, rmse, highlight_sims)) |>
  plot_grid(plotlist=_, ncol=2)
map(list(metrics_by_farm_mn, metrics_by_week_mn),
    ~plot_metric_ordered(.x, rmse, highlight_sims)) |>
  plot_grid(plotlist=_, ncol=2)


metric_ranks <- bind_rows(
  metrics_by_farm |>
    filter(sim %in% c("predFcst", "predBlend",
                      paste0("sim_0", 1:9), paste0("sim_", 10:20),
                      paste0("sLonLatD", 3:4, "_n5"),
                      paste0("sLonLatD", 3:4, "_n10"),
                      paste0("sLonLatD", 3:4, "_n20"))) |>
    select(sepaSite, sim, N, rmse) |>
    mutate(value_lowGood=rmse,
           type="byFarm") |>
    drop_na() |>
    group_by(sepaSite) |>
    mutate(rank=min_rank(value_lowGood)) |>
    ungroup(),
  metrics_by_week |>
    filter(sim %in% c("predFcst", "predBlend",
                      paste0("sim_0", 1:9), paste0("sim_", 10:20),
                      paste0("sLonLatD", 3, "_n5"),
                      paste0("sLonLatD", 3, "_n10"),
                      paste0("sLonLatD", 3, "_n20"))) |>
    select(date, sim, N, rmse) |>
    mutate(value_lowGood=rmse,
           type="byWeek") |>
    drop_na() |>
    group_by(date) |>
    mutate(rank=min_rank(value_lowGood)) |>
    ungroup()
  ) |>
  left_join(sim_i)

metric_ranks |>
  filter(N >= 30) |>
  group_by(sim, lab, lab_short, type) |>
  summarise(mn=median(rank, na.rm=T)) |>
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
  filter(N >= 30) |>
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
  metrics_by_farm_md |> mutate(type="byFarm"),
  metrics_by_week_md |> mutate(type="byWeek")
) |>
  pivot_longer(any_of(c("rmse", "rsq", "rho", "r", "ROC_AUC", "PR_AUC")), names_to="metric") |>
  filter(metric %in% c("rmse")) |>
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


ms_rmse <- all_metrics_df |> filter(metric=="RMSE") |>
  metric_plot_base(theme="ms") + 
  scale_y_continuous("Median cross-validation RMSE", limits=c(0.22, 0.4), oob=scales::oob_keep, 
                     breaks=seq(0, 1, by=0.05), minor_breaks=seq(0, 1, by=0.01))

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
  scale_size_manual(values=c(rep(2.5, 4), 1.5, 1.5, rep(1, 2))) +
  scale_alpha_manual(values=c(1, 1, 1, 1, 1, 1, 0.5, 0.5)) +
  ylim(0.575, 1.175) +
  theme(legend.position="none",
        plot.margin=margin(t=0, b=0, l=0, r=0),
        panel.border=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())
p <- plot_grid(ms_rmse, ms_legend, 
               align="h", axis="tb", nrow=1, rel_widths=c(1.12, 0.4))
ggsave("figs/pub/validation_metrics_CV_medians.png", p, width=3, height=4)






# weekly performance ------------------------------------------------------

metric_date_df <- ensCV_df |>
  filter(date >= "2021-05-01") |>
  group_by(date) |>
  summarise(rmse=rmse_vec(IP_predBlend, truth=licePerFish_rtrt)) |>
  ungroup() |>
  mutate(sim="predBlend") |>
  bind_rows(
    ensCV_df |>
      filter(date >= "2021-05-01") |>
      group_by(date) |>
      summarise(rmse=rmse_vec(IP_predFwk_RMSE, truth=licePerFish_rtrt)) |>
      ungroup() |>
      mutate(sim="predFcst")
  ) |>
  bind_rows(
    ensCV_df |>
      filter(date >= "2021-05-01") |>
      group_by(date) |>
      summarise(rmse=rmse_vec(IP_sim_20, truth=licePerFish_rtrt)) |>
      ungroup() |>
      mutate(sim="sim_20")
  ) |>
  bind_rows(
    ensCV_df |>
      filter(date >= "2021-05-01") |>
      group_by(date) |>
      summarise(rmse=rmse_vec(IP_sim_04, truth=licePerFish_rtrt)) |>
      ungroup() |>
      mutate(sim="sim_04")
  ) |>
  pivot_longer(2, names_to="metric") |>
  mutate(metric=factor(metric, levels=c("ROC_AUC", "PR_AUC", "rsq", "r", "rmse"),
                       labels=c("'AUC'['ROC']", "'AUC'['PR']", "R^2", "rho", "RMSE"))) |>
  left_join(sim_i) |>
  droplevels()

metric_date_labs <- expand_grid(date=max(metric_date_df$date),
                                metric=levels(metric_date_df$metric),
                                lab=levels(metric_date_df$lab)) |>
  mutate(value=c(0.21, 0.32, 0.305, 0.33)) |>
  filter(metric %in% c("RMSE", "'AUC'['ROC']"))

p <- metric_date_df |>
  filter(metric %in% c("RMSE", "'AUC'['ROC']")) |>
  ggplot(aes(date, value, colour=lab)) + 
  geom_point(size=0.5, shape=1) + 
  geom_line(stat="smooth", method="gam", formula=y~s(x), se=F) +
  scale_x_date(date_breaks="1 year", #date_minor_breaks="3 months", 
               date_labels="%Y", expand=expansion(mult=c(0.05, 0.1))) +
  ylab("Cross validation RMSE by week") +
  scale_colour_manual("Model", 
                      values=c("black", "red",
                               scico(2, begin=0.2, end=0.7, palette="broc", direction=1)), 
                      labels=c(bquote(Ens['Fcst']), 
                               bquote(Ens['Blend']), 
                               "Opt: 3D",
                               "Opt: 2D")) +
  theme_bw() + 
  guides(colour=guide_legend(override.aes=list(size=1))) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size=9),
        panel.grid.major.y=element_line(colour="grey85", linewidth=0.4),
        panel.grid.minor.y=element_line(colour="grey90", linewidth=0.2),
        axis.text=element_text(size=8))
ggsave("figs/pub/validation_metrics_byWeek.png", p, width=6.5, height=4)



# Performance summaries: values
all_metrics_df |> 
  filter(sim %in% c("predBlend", "predFcst", "sim_04", "sim_20")) |> 
  arrange(type, value)








# farm performance --------------------------------------------------------

metric_farm_df <- ensCV_df |>
  filter(date >= "2021-05-01") |>
  group_by(sepaSite) |>
  summarise(rmse=rmse_vec(IP_predBlend, truth=licePerFish_rtrt)) |>
  ungroup() |>
  mutate(sim="predBlend") |>
  bind_rows(
    ensCV_df |>
      filter(date > "2021-05-01") |>
      group_by(sepaSite) |>
      summarise(rmse=rmse_vec(IP_predFf_RMSE, truth=licePerFish_rtrt)) |>
      ungroup() |>
      mutate(sim="predFcst")
  ) |>
  bind_rows(
    ensCV_df |>
      filter(date >= "2021-05-01") |>
      group_by(sepaSite) |>
      summarise(rmse=rmse_vec(IP_sim_20, truth=licePerFish_rtrt)) |>
      ungroup() |>
      mutate(sim="sim_20")
  ) |>
  bind_rows(
    ensCV_df |>
      filter(date >= "2021-05-01") |>
      group_by(sepaSite) |>
      summarise(rmse=rmse_vec(IP_sim_04, truth=licePerFish_rtrt)) |>
      ungroup() |>
      mutate(sim="sim_04")
  ) |>
  pivot_longer(2, names_to="metric") |>
  mutate(metric=factor(metric, levels=c("ROC_AUC", "PR_AUC", "rsq", "r", "rmse"),
                       labels=c("'AUC'['ROC']", "'AUC'['PR']", "R^2", "rho", "RMSE"))) |>
  left_join(sim_i) |>
  droplevels()




# Blending proportions ------------------------------------------------------

mod <- "n20_sLonLatD3"
out_ensBlend <- readRDS(glue("out/ensembles/ensBlend_{mod}_FULL_stanfit.rds"))
dat_ensBlend <- readRDS(glue("out/ensembles/ensBlend_{mod}_FULL_standata.rds"))
ensFull_LatLon <- read_csv("out/valid_df_2021-2024.csv") |>
  select(rowNum, date, CV_k, sepaSite, sepaSiteNum, licePerFish_rtrt, starts_with("sim")) |>
  select(-contains("avg")) |>
  mutate(across(starts_with("sim_"), ~.x - mean(.x), .names="c_{.col}")) |>
  left_join(site_i) |>
  select(-sepaSite) |>
  arrange(rowNum)
ensBlend_rec <- make_spline_recipe(ensFull_LatLon, 3, sim_i$sim[1:20])
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_bbox <- st_bbox(mesh_fp)
mesh_land <- st_convex_hull(mesh_fp) |>
  st_difference(mesh_fp) |>
  st_crop(site_i |> st_as_sf(coords=c("easting", "northing"), crs=27700) |> st_buffer(10e3))


# blending proportions by site
map_df <- site_i |>
  mutate(sepaSiteNum=row_number()) |>
  bind_cols(ensFull_LatLon |> summarise(across(c(licePerFish_rtrt, contains("sim")), mean))) 
b_p_ls <- make_predictions_ensBlend_sLonLat(out_ensBlend, 
                                            newdata=bake(ensBlend_rec, map_df), 
                                            iter=3000, mode="b_p") 
b_p_post <- map_dfr(1:dim(b_p_ls)[2], 
                      ~as_tibble(b_p_ls[,.x,]) |>
                        set_names(dat_ensBlend$sim_names) |>
                        mutate(rowNum=row_number(),
                               iter=.x)) |> 
  pivot_longer(starts_with("sim"), names_to="sim", values_to="p") |>
  left_join(sim_i) |>
  rename(Simulation=lab)
gc()
  
p_a <- b_p_post |> 
  ggplot(aes(p, group=rowNum)) + 
  geom_line(alpha=0.2, linewidth=0.25, stat="density", adjust=2) +
  labs(x=expression(paste("Ensemble blending weight (", italic(pi[~~k]), ")")),
       y="log density") +
  scale_y_continuous(transform="log1p") +
  facet_wrap(~Simulation, labeller=label_parsed, scales="free_y", ncol=1, strip.position="right") +
  theme(axis.title=element_text(size=9),
        axis.text=element_text(size=7),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave("figs/pub/ensBlend_p_sitePosterior.png", p_a, width=4, height=16)



# blending proportions by parameter: maps
map_df <- expand_grid(easting=seq(min(site_i$easting)-30e3, max(site_i$easting)+30e3, by=4e3),
                      northing=seq(min(site_i$northing)-30e3, max(site_i$northing)+30e3, by=4e3)) |>
  st_as_sf(coords=c("easting", "northing"), crs=27700, remove=F) |>
  st_intersection(st_buffer(mesh_fp, 5e3)) |>
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
sim_key <- read_csv("out/sim_2021-2024/sim_i.csv") |> 
  mutate(sim=paste0("sim_", i)) |>
  select(-outDir) |> 
  inner_join(sim_i) |>
  mutate(across(matches("swim|Thresh"), ~if_else(lab_short=="2D", NA, .x)),
         across(matches("swim"), ~abs(.x)),
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


param_post_sum |>
  pivot_longer(-(1:2), names_to="var", values_to="value") |>
  ggplot(aes(value)) + 
  geom_density() + 
  facet_wrap(~var, scales="free")

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
                                   "lightThreshNauplius", "lightThreshCopepodid")))

param_map_ls <- param_post_sum |>
  pivot_longer(-(1:2), names_to="var", values_to="value") |>
  group_by(rowNum, var) |>
  summarise(post_mn=mean(value)) |>
  ungroup() |>
  inner_join(var_pretty, by=join_by(var)) |>
  arrange(var_order) |>
  full_join(map_df |> select(easting, northing) |> mutate(rowNum=row_number()),
            by=join_by(rowNum)) |>
  group_split(var_order, .keep=TRUE)




param_map_plot_ls <- map(param_map_ls, ~make_param_map_plot(.x, site_i, mesh_land))
p <- plot_grid(plotlist=param_map_plot_ls, align="hv", axis="tblr", nrow=2)
ggsave("figs/pub/ensBlend_p_map-mn.png", p, height=10.2, width=12, dpi=200)




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

sim_params |> 
  ggplot(aes(param_val, i, colour=i=="04")) + 
  geom_point() + 
  facet_wrap(~param_name, scales="free")

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
mod <- "n20_sLonLatD3"
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
ensBlend_rec <- make_spline_recipe(IP_LatLon, 3, sim_i$sim[1:20])
b_p_ls <- make_predictions_ensBlend_sLonLat(out_ensBlend, 
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

influx_df <- readRDS("out/sim_2021-2024/processed/connectivity_day.rds") |>
  select(sepaSite, sepaSite, date, sim, influx_m2) |>
  mutate(sim=paste0("sim_", sim),
         influx_m3_4rt=(replace_na(influx_m2, 0)/20)^0.25)

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
ggsave("figs/pub/ensBlend_influx_daily_20.png", fig_influx, width=7, height=4, dpi=400)

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

set.seed(1003)
mod <- "n20_sLonLatD3"
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
                  ens_mn_orig=c(0,0),
                  ens_CL005_orig=c(0,0),
                  ens_CL025_orig=c(0,0),
                  ens_CL975_orig=c(0,0),
                  ens_CL995_orig=c(0,0),
                  ens_CI95width_orig=c(0,0),
                  ens_CI99width_orig=c(0,0))
p_dir <- "out/ensembles/p_meshCentroids/"

library(doFuture)
plan(multicore, workers=20)
for(i in 1:length(f)) {
  timestep <- ymd("2021-01-01") + dhours(as.numeric(str_sub(str_split_fixed(f[i], "_t_", 2)[,2], 1, -5)))
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
saveRDS(ps_lims, "out/sim_2021-2024/processed/ps_lims.rds")

timesteps <- ymd("2019-04-01") + dhours(as.numeric(str_sub(str_split_fixed(f, "_t_", 2)[,2], 1, -5)))
ens_df <- map2_dfr(ens_ls, timesteps, ~.x |> mutate(date=.y))
saveRDS(ens_df, "out/sim_2021-2024/processed/ens_weekly.rds")

ens_avg <- ens_df |>
  filter(date >= "2021-05-01") |>
  group_by(i) |>
  summarise(across(where(is.numeric), .fn=list(mn=mean, md=median))) |>
  ungroup()
saveRDS(ens_avg, "out/sim_2021-2024/processed/ens_avg_sLonLatD3_n20.rds")




# maps --------------------------------------------------------------------

# Left side: Ensemble mean(copepodid density)
# Right side: Ensemble mean(weekly CI width)
# ens_df <- readRDS("out/sim_2021-2024/processed/ens_weekly.rds")
ens_avg <- readRDS("out/sim_2021-2024/processed/ens_avg_sLonLatD3_n20.rds")

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
  ggsave("figs/pub/ens_map_sLonLatD3_n20.png", plot=_, width=4.75, height=9.1, dpi=600)

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




# rmse stormcloud ---------------------------------------------------------

# ensCV_df <- read_csv("out/ensemble_CV.csv")

# sim_04 would be selected as 'optimal'
metric_ranks |> 
  group_by(type, sim) |> 
  summarise(mnRank=mean(rank)) |> 
  group_by(sim) |> 
  summarise(mn=mean(mnRank)) |> 
  arrange(mn)

farm_r.df <- metrics_by_farm |>
  filter(N >= 30) |>
  filter(sim %in% c("predFcst", "predBlend", "sim_04", "sim_17", "sim_avg3D", "sim_avg2D")) |>
  # filter(grepl("null|avg|pred", sim)) |>
  # filter(sim != "nullFarm") |>
  left_join(sim_i) |>
  droplevels() |>
  select(sepaSite, sim, rmse, lab, lab_short) |>
  pivot_longer("rmse", names_to="metric") |>
  mutate(type="By farm") |>
  mutate(metric=factor(metric, levels=c("ROC_AUC", "PR_AUC", "rsq", "r", "rmse"),
                     labels=c("'AUC'['ROC']", "'AUC'['PR']", "R^2", "rho", "RMSE"))) |>
  arrange(lab) |>
  mutate(type=factor(type, 
                     levels=c("global", "By farm", "By week"),
                     labels=c("Global", "'By farm'", "'By week'"))) |>
  filter(!is.na(value)) |>
  mutate(lab=lvls_revalue(lab, c("Ens['Fcst']", "Ens['Blend']", "Mean['3D']", "Mean['2D']", "Opt['3D']", "Opt['2D']")),
         lab=lvls_reorder(lab, c(1,2,5,3,6,4)))
week_r.df <- metrics_by_week |>
  filter(N >= 30) |>
  filter(sim %in% c("predFcst", "predBlend", "sim_04", "sim_17", "sim_avg2D", "sim_avg3D")) |>
  # filter(grepl("null|avg|pred", sim)) |>
  # filter(sim != "nullTime") |>
  left_join(sim_i) |>
  droplevels() |>
  select(date, sim, rmse, lab, lab_short) |>
  pivot_longer("rmse", names_to="metric") |>
  filter(!(sim=="nullTime" & metric=="ROC_AUC")) |>
  mutate(type="By week") |>
  mutate(metric=factor(metric, levels=c("ROC_AUC", "PR_AUC", "rsq", "r", "rmse"),
                       labels=c("'AUC'['ROC']", "'AUC'['PR']", "R^2", "rho", "RMSE"))) |>
  arrange(lab) |>
  mutate(type=factor(type, 
                     levels=c("global", "By farm", "By week"),
                     labels=c("Global", "'By farm'", "'By week'"))) |>
  filter(!is.na(value)) |>
  mutate(lab=lvls_revalue(lab, c("Ens['Fcst']", "Ens['Blend']", "Mean['3D']", "Mean['2D']", "Opt['3D']", "Opt['2D']")),
         lab=lvls_reorder(lab, c(1,2,5,3,6,4)))
mn_ci <- bind_rows(farm_r.df, week_r.df) |>
  # filter(lab %in% c("Ens['Fcst']", "Ens['Blend']", "Opt['3D']", "Opt['2D']")) |>
  # filter(!lab %in% c("Mean3D", "Mean2D")) |>
  droplevels() |>
  group_by(sim, lab, lab_short, metric, type) |>
  summarise(mn=mean(value, na.rm=T),
            md=median(value, na.rm=T),
            N=sum(!is.na(value)),
            se=sd(value, na.rm=T)/sqrt(N),
            ci_lo=mn - qt(0.975, N-1)*se,
            ci_hi=mn + qt(0.975, N-1)*se)

all_metrics_medians <- all_metrics_df |>
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
  droplevels() |>
  ggplot(aes(value, lab, fill=lab_short, colour=lab_short)) + 
  geom_dots(side="bottom", scale=0.5) + 
  stat_slab(normalize="xy", scale=0.5, colour=NA, fill_type="gradient",
            aes(slab_alpha=after_stat(-pmax(abs(1-2*cdf), 0.25)))) +
  stat_pointinterval(.width=c(0.5, 0.8), colour="black", fatten_point=1.25) +
  # geom_rug(data=all_metrics_medians |> filter(sim %in% mn_ci$sim), 
  #          aes(x=value), sides="b", length=unit(0.075, "npc"), linewidth=0.5, alpha=0.75) + 
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
  labs(x="Cross validation RMSE") +
  facet_grid(type~., scales="free_x", labeller="label_parsed") +
  theme(panel.grid.major=element_line(colour="grey90", linewidth=0.2),
        strip.text=element_text(size=11),
        axis.title.x=element_text(size=9),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=9),
        legend.position="none")
ggsave("figs/pub/RMSE_stormclouds.png", p, width=10, height=5.25, dpi=300)





# rank stormcloud ---------------------------------------------------------

# ensCV_df <- read_csv("out/ensemble_CV.csv")

p <- metric_ranks |>
  filter(N > 30) |>
  mutate(lab=factor(lab, 
                    levels=c("Null['farm']", "Null['time']", 
                             paste0("'2D.", c(2, 3, 1, 4), "'"),
                             "Mean2D",
                             paste0("'3D.", c(3, 7, 6, 12, 2, 13, 11, 15, 
                                              5, 16, 14, 1, 8, 9, 10, 4), "'"),
                             "Mean3D",
                             paste0("Ens['", c("Blend", "Fcst"), "']")) |>
                      rev())) |>
  mutate(metric=factor(metric, levels=c("ROC_AUC", "PR_AUC", "rsq", "r", "rho", "rmse"),
                       labels=c("'AUC'['ROC']", "'AUC'['PR']", "R^2", "r", "rho", "RMSE"))) |>
  arrange(lab) |>
  mutate(type=factor(type, 
                     levels=c("global", "byFarm", "byWeek"),
                     labels=c("Global", "'By farm'", "'By week'"))) |>
  ggplot(aes(rank, lab, fill=lab_short, colour=lab_short)) + 
  geom_dots(side="bottom", scale=0.5) + 
  stat_histinterval(normalize="xy", scale=0.5, colour=NA, alpha=0.5, 
                    breaks=breaks_fixed(width=2)) +
  stat_pointinterval(.width=c(0.5, 0.95), colour="black", fatten_point=1.2) +
  scale_fill_scico_d(palette="glasgow", guide="none", end=0.8) +
  scale_colour_scico_d(palette="glasgow", guide="none", end=0.8) +
  labs(x="Rank") +
  facet_grid(type~., scales="free_x", labeller="label_parsed") +
  theme(panel.grid.major=element_line(colour="grey90", linewidth=0.2),
        strip.text=element_text(size=11),
        axis.title.x=element_text(size=9),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=9))
ggsave("figs/pub/RMSE_stormclouds_ranks_all.png", p, width=10, height=10, dpi=300)

metric_ranks |>
  filter(N >= 30) |>
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
  filter(N > 30) |>
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

rmse_info <- tibble(breaks=seq(0.1, 0.65, by=0.05)) |>
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
mods <- c("IP_predFf_RMSE", "IP_predBlend", "IP_sim_avg3D", "IP_sim_avg2D")
col_labs <- c(expression(Ens['Fcst']), expression(Ens['Blend']), 
              expression(Mean['3D']), expression(Mean['2D']))
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
    geom_polygon(data=farm_rmse_bar.df, aes(x, y, fill=mdpt, group=mdpt),
                 colour="grey10", linewidth=0.15) +
    geom_text(data=farm_rmse_count_labs, aes(x, y, label=prop),
              size=2, hjust=1, vjust=0.5, nudge_x=-1000) +
    colorspace::scale_fill_binned_diverging(
      name=col_lab, palette="Blue-Red 3", l1=20, l2=90, p2=2, mid=0.375,
      limits=c(0.1, 0.65), breaks=rmse_info$breaks, labels=rmse_info$break_labs) +
    scale_y_continuous(limits=c(630000, 955000), oob=scales::oob_keep,
                       breaks=c(56, 58), labels=paste0(c(56, 58), "\u00B0N")) +
    scale_x_continuous(breaks=c(-7, -5), labels=paste0(c(7, 5), "\u00B0W"),
                       limits=c(75000, 235000), oob=scales::oob_keep) +
    theme(legend.position="inside",
          legend.position.inside=c(c(0.198,0.19,0.195,0.195)[i], 0.202),
          legend.background=element_blank(),
          legend.key.height=unit(0.43, "cm"),
          legend.key.width=unit(0.0, "cm"),
          legend.text=element_text(size=6),
          legend.title=element_text(size=8, vjust=1, hjust=1),
          legend.ticks=element_line(colour="grey10", linewidth=0.25),
          legend.ticks.length=unit(0.04, "cm"),
          axis.title=element_blank())
  if(i > 1) p_ls[[i]] <- p_ls[[i]] + theme(axis.text.y=element_blank())
}

cowplot::plot_grid(plotlist=p_ls, labels="auto", align="hv", axis="tblr", nrow=1) |>
  ggsave("figs/pub/ens_farm-rmse+IDW_map.png", plot=_ , width=13, height=6, dpi=300)


ggsave("figs/talk/rmse+IDW_ensFcst.png", p_ls[[1]], width=3, height=5.5)
ggsave("figs/talk/rmse+IDW_ensBlend.png", p_ls[[2]], width=3, height=5.5)
ggsave("figs/talk/rmse+IDW_ens3D.png", p_ls[[3]], width=3, height=5.5)
ggsave("figs/talk/rmse+IDW_ens2D.png", p_ls[[4]], width=3, height=5.5)






# maps of AUC + IDW ---------------------------------------------------------

library(terra)

mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_rast <- st_read("data/WeStCOMS2_meshFootprint.gpkg") |>
  rast(resolution=500)


farm_RMSE.df <- metrics_by_farm |>
  filter(N >= 30) |>
  # filter(grepl("avg|pred", sim)) |>
  filter(sim %in% c("predFcst", "predBlend", "sim_04", "sim_20")) |>
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
  ggsave("figs/pub/ens_AUC-RMSE_map_D3n20_NEW.png", plot=_ , width=10, height=10, dpi=300)








# bad performers ----------------------------------------------------------

metrics_by_farm |> 
  filter(N >= 30) |>
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
