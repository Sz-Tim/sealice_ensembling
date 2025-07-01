# Project: Sealice IP Ensemble
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Blending ensemble



# setup -------------------------------------------------------------------

library(tidyverse)
library(glue)
library(rstan)
library(recipes)
rstan_options(auto_write=T)
source("code/00_fn.R")
theme_set(theme_bw())


# load datasets -----------------------------------------------------------

# Full dataset
ensFull_df <- read_csv("out/valid_df_2021-2024_FULL.csv") |>
  mutate(across(starts_with("sim_"), ~.x - mean(.x), .names="c_{.col}")) 
folds <- unique(ensFull_df$CV_k)

sim_i <- read_csv("out/sim_2021-2024/sim_i.csv", show_col_types=F) |>
  mutate(sim=paste0("sim_", i),
         lab_short=if_else(fixDepth, "2D", "3D")) |>
  group_by(lab_short) |>
  mutate(lab=paste0(lab_short, ".", row_number())) |>
  ungroup() |>
  select(sim, lab_short, lab) |>
  bind_rows(
    tibble(sim=c("predFcst", "predBlend", "sim_avg2D", "sim_avg3D", "null"),
           lab_short=c("Ens['Fcst']", "Ens['Blend']", "Mean2D", "Mean3D", "Null"),
           lab=c("Ens['Fcst']", "Ens['Blend']", "Mean2D", "Mean3D", "Null"))
  ) |>
  mutate(lab=factor(lab, 
                    levels=c("Ens['Fcst']", "Ens['Blend']", "Mean2D", "Mean3D", 
                             paste0("2D.", 1:20), paste0("3D.", 1:20), "Null")),
         lab_short=factor(lab_short, 
                          levels=c("Ens['Fcst']", "Ens['Blend']", "Mean2D", "Mean3D", 
                                   "2D", "3D", "Null")))

site_i <- read_csv("data/farm_sites.csv")
ensFull_LatLon <- read_csv("out/valid_df_2021-2024_FULL.csv") |>
  select(rowNum, date, CV_k, sepaSite, sepaSiteNum, licePerFish_rtrt, starts_with("sim")) |>
  select(-contains("avg")) |>
  mutate(across(starts_with("sim_"), ~.x - mean(.x), .names="c_{.col}")) |>
  left_join(site_i) |>
  select(-sepaSite) |>
  arrange(rowNum)



# cross-validation --------------------------------------------------------

set.seed(1001)
mods <- expand_grid(nSims=c(paste0("n", c(3, 5, 10, 15))),
            sample=1:100) |>
  rowwise() |>
  mutate(sim=list(sample(1:20, as.numeric(str_sub(nSims, 2, -1))))) |>
  ungroup() |>
  mutate(mod=list(paste0("sLonLat_", c("", "RE_"), "GQ_D3"))) |>
  unnest(mod) |>
  filter(mod=="sLonLat_RE_GQ_D3")
# saveRDS(mods, "out/ensembles/ensEval/ensBlend_EvalModSpecs.rds")

# mods <- mods |>
#   filter(nSims=="n3") |>
#   arrange(desc(mod), (sample))

CV_ensBlend <- CV_ensAvg <- vector("list", nrow(mods))

for(i in rev(1:nrow(mods))) {
  
  CV_k_ls <- CV_k_ensAvg_ls <- vector("list", length(folds))
  
  for(k in seq_along(folds)) {
    # EnsBlend setup
    recipe_i <- make_spline_recipe(ensFull_LatLon, 
                                   as.numeric(str_split_fixed(mods$mod[i], "D", 2)[2]), 
                                   sim_i$sim[mods$sim[[i]]])
    full_df <- bake(recipe_i, ensFull_LatLon)
    test_rows <- which(ensFull_LatLon$CV_k == folds[k])
    dat_rstan <- make_data_rstan_sLonLat_GQ(full_df, test_rows)
    pars <- c(#"b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu", #"b_p",
              "GQ_Ypred"#, "GQ_mu", "GQ_hu", "GQ_IP_ens",
              # paste0("b_s_", c("easting", "northing", "easting_x_northing"))
              )
    
    # fit EnsBlend
    fname <- glue("out/ensembles/ensEval/ensBlend_{mods$mod[i]}_{mods$nSims[i]}-{mods$sample[i]}_CV-{folds[k]}")
    if(file.exists(glue("{fname}_stanfit.rds"))) {
      cat("File exists:", fname, "\n")
      out_ensBlend <- readRDS(glue("{fname}_stanfit.rds"))
    } else {
      stanMod <- str_sub(str_split_fixed(mods$mod[i], "D", 2)[1], 1, -2)
      out_ensBlend <- stan(file=glue("code/stan/ensemble_mixture_model_{stanMod}.stan"),
                           model_name=mods$mod[i], data=dat_rstan,
                           chains=3, cores=3, iter=3000, warmup=2000,
                           control=list(adapt_delta=0.95, max_treedepth=20),
                           pars=pars)
      colMeans(rstan::extract(out_ensBlend, pars="GQ_Ypred")[[1]]) |>
        as_tibble() |>
        set_names(paste0("IP_", mods$mod[i], "_", mods$nSims[i], "_", mods$sample[i])) |>
        mutate(rowNum=test_rows) |>
        saveRDS(glue("{fname}_Ypred.rds"))
      # saveRDS(out_ensBlend, glue("{fname}_stanfit.rds"))
      saveRDS(dat_rstan, glue("{fname}_standata.rds")) 
    }
   
    # fit using avg AEIP
    dat_avg_df <- ensFull_df |>
      select(rowNum, sepaSite, sepaSiteNum, CV_k, date, licePerFish_rtrt,
             all_of(sim_i$sim[mods$sim[[i]]])) |>
      rowwise() |>
      mutate(sim_avg=mean(c_across(starts_with("sim")))) |>
      ungroup() |>
      mutate(c_sim_avg=c(scale(sim_avg))) |>
      select(-matches("sim_[0-9]"))
    test_rows <- which(dat_avg_df$CV_k == folds[k])
    fname_avg <- glue("out/ensembles/ensEval/ensAvg_{mods$nSims[i]}-{mods$sample[i]}_CV-{folds[k]}")
    if(file.exists(glue("{fname_avg}_stanfit.rds"))) {
      cat("File exists:", fname_avg, "\n")
      out_sim <- readRDS(glue("{fname_avg}_stanfit.rds"))
    } else {
      dat_rstan <- dat_avg_df |> make_data_rstan_GQ(test_rows)
      out_sim <- stan(file="code/stan/candidate_model_GQ.stan",
                      model_name=glue("avg-{folds[k]}"), data=dat_rstan,
                      chains=3, cores=3, iter=3000, warmup=2000,
                      pars=c("GQ_Ypred"))
      colMeans(rstan::extract(out_sim, pars="GQ_Ypred")[[1]]) |>
        as_tibble() |>
        set_names(paste0("IP_avg", "_", mods$nSims[i], "_", mods$sample[i])) |>
        mutate(rowNum=test_rows) |>
        saveRDS(glue("{fname_avg}_Ypred.rds"))
      # saveRDS(out_sim, glue("{fname_avg}_stanfit.rds"))
      saveRDS(dat_rstan, glue("{fname_avg}_standata.rds"))
    }
  }
}

# load predictions
pred_f <- tibble(pred_f=dir("out/ensembles/ensEval", "Ypred"),
                 dat_f=str_replace(pred_f, "Ypred", "standata"),
                 mod=str_split_fixed(pred_f, "_", 2)[,1],
                 nSims=str_split_fixed(str_split_fixed(pred_f, "_n", 2)[,2], 
                                       "-", 2)[,1],
                 sample=str_split_fixed(str_split_fixed(pred_f, "-", 2)[,2],
                                        "_", 2)[,1],
                 fold=str_split_fixed(str_split_fixed(pred_f, "CV-", 2)[,2],
                                      "_", 2)[,1]) |>
  group_by(mod, nSims, sample) |>
  mutate(finished=any(grepl(10, fold))) |>
  ungroup() |>
  filter(finished)

predBlend_df <- pred_f |>
  filter(mod=="ensBlend") |>
  mutate(preds=map(pred_f, ~readRDS(paste0("out/ensembles/ensEval/", .x)))) |>
  select(-pred_f, -dat_f, -fold) |>
  group_by(mod, nSims, sample) |>
  nest(dat=preds) |>
  rowwise() |>
  mutate(df=list(map_dfr(dat, ~.x))) |>
  ungroup() 
ensBlend_preds <- reduce(predBlend_df$df, full_join, by=join_by(rowNum))
ensBlend_preds |>
  write_csv("out/ensembles/CV_ensBlend_EvalCV_NEW_TEMP.csv")

predAvg_df <- pred_f |>
  filter(mod=="ensAvg") |>
  mutate(preds=map(pred_f, ~readRDS(paste0("out/ensembles/ensEval/", .x)))) |>
  select(-pred_f, -dat_f, -fold) |>
  group_by(mod, nSims, sample) |>
  nest(dat=preds) |>
  rowwise() |>
  mutate(df=list(map_dfr(dat, ~.x))) |>
  ungroup() 
ensAvg_preds <- reduce(predAvg_df$df, full_join, by=join_by(rowNum))
ensAvg_preds |>
  write_csv("out/ensembles/CV_ensAvg_EvalCV_NEW_TEMP.csv")






# compile predictions -----------------------------------------------------

mod_ids <- readRDS("out/ensembles/ensEval/ensBlend_EvalModSpecs.rds") |>
  inner_join(mods |> select(nSims, sample, mod)) |>
  mutate(ens_id=paste("IP", str_sub(mod, -2, -1), nSims, sample, sep="_"),
         sim_cols=map(sim,
                      ~tibble(s=paste0("IP_sim_", str_pad(.x, 2, "left", "0"))) |>
                        mutate(cand_id=paste0("cand", row_number())) |>
                        pivot_wider(names_from=cand_id, values_from=s))) |>
  bind_rows(tibble(ens_id=paste0("IP_D", 3:10, "_n20_1"),
                   mod=paste0("sLonLatD", 3:10),
                   nSims="n20",
                   sample=1,
                   sim=list(c(1:20)),
                   sim_cols=list(map_dfc(1:20,
                                         ~tibble(x=paste0("IP_sim_", str_pad(.x, 2, "left", "0"))) |>
                                           set_names(paste0("cand", .x)))))
  ) |>
  mutate(nSims=factor(nSims, levels=c("n3", "n5", "n10", "n15", "n20"))) |>
  arrange(nSims, sample, mod) |>
  mutate(ens_id_ord=factor(ens_id, levels=unique(ens_id)))
null_df <- read_csv("out/ensemble_CV.csv") |>
  select(rowNum, "IP_null0") |>
  pivot_longer(starts_with("IP_"), names_to="null_id", values_to="null_pred") |>
  inner_join(ensFull_LatLon |>
               select(rowNum, licePerFish_rtrt, date, sepaSiteNum))
candidate_df <- read_csv("out/candidates/CV_candidate_predictions.csv") |>
  pivot_longer(starts_with("IP_"), names_to="s", values_to="cand_pred") |>
  inner_join(ensFull_LatLon |>
               select(rowNum, licePerFish_rtrt, date, sepaSiteNum)) |>
  mutate(lice_g05=factor(licePerFish_rtrt^4 > 0.5))
avg_df <- read_csv("out/ensembles/CV_ensAvg_EvalCV_NEW_TEMP.csv") |>
  left_join(read_csv("out/ensembles/CV_avg_predictions.csv") |>
              select(rowNum, IP_sim_avg3D) |> rename(IP_avg_n20_1=IP_sim_avg3D)) |>
  pivot_longer(starts_with("IP_"), names_to="avg_id", values_to="avg_pred") |>
  separate_wider_delim(avg_id, "_", names=c("x", "avg", "nSim", "sample")) |>
  select(-x, -avg)
ens20_df <- read_csv("out/ensembles/CV_ensBlend_predictions.csv") |>
  select(rowNum, matches("IP_sLonLat.*n20")) |>
  inner_join(ensFull_LatLon |>
               select(rowNum, licePerFish_rtrt, date, sepaSiteNum)) |>
  mutate(lice_g05=factor(licePerFish_rtrt^4 > 0.5)) |>
  pivot_longer(starts_with("IP_"), names_to="ens_id", values_to="ens_pred") |>
  mutate(ens_id=paste0(str_remove(ens_id, "sLonLat"), "_1")) 
CV_df <- read_csv("out/ensembles/CV_ensBlend_EvalCV_NEW_TEMP.csv") |>
  rename_with(~str_remove(.x, "sLonLat_RE_GQ_")) |>
  inner_join(ensFull_LatLon |>
               select(rowNum, licePerFish_rtrt, date, sepaSiteNum)) |>
  mutate(lice_g05=factor(licePerFish_rtrt^4 > 0.5)) |>
  pivot_longer(starts_with("IP_"), names_to="ens_id", values_to="ens_pred") |>
  # bind_rows(ens20_df) |>
  separate_wider_delim(ens_id, "_", names=c("x", "D", "nSim", "sample"), cols_remove=F) |>
  select(-x) |>
  full_join(avg_df, by=join_by(rowNum, nSim, sample))

library(yardstick)
summary_fn <- list(means=function(x) {mean(x, na.rm=T)},
                   medians=function(x) {median(x, na.rm=T)})[2]
# RMSE: Site
RMSE_null_site <- null_df |>
  group_by(sepaSiteNum) |>
  mutate(N=n()) |>
  filter(N >= 30) |>
  summarise(across(contains("pred"), ~rmse_vec(.x, truth=licePerFish_rtrt))) |>
  ungroup()
RMSE_cand_site <- candidate_df |>
  group_by(sepaSiteNum, s) |>
  mutate(N=n()) |>
  filter(N >= 30) |>
  summarise(across(contains("pred"), ~rmse_vec(.x, truth=licePerFish_rtrt))) |>
  ungroup()
RMSE_ens_site <- CV_df |>
  group_by(sepaSiteNum, ens_id, D, nSim, sample) |>
  mutate(N=n()) |>
  filter(N >= 30) |>
  summarise(across(contains("pred"), ~rmse_vec(.x, truth=licePerFish_rtrt))) |>
  ungroup() |>
  left_join(mod_ids |>
              select(ens_id, sim_cols) |>
              unnest_wider(sim_cols))

# RMSE: Date
RMSE_null_date <- null_df |>
  group_by(date) |>
  mutate(N=n()) |>
  filter(N >= 30) |>
  summarise(across(contains("pred"), ~rmse_vec(.x, truth=licePerFish_rtrt))) |>
  ungroup()
RMSE_cand_date <- candidate_df |>
  group_by(date, s) |>
  mutate(N=n()) |>
  filter(N >= 30) |>
  summarise(across(contains("pred"), ~rmse_vec(.x, truth=licePerFish_rtrt))) |>
  ungroup()
RMSE_ens_date <- CV_df |>
  group_by(date, ens_id, D, nSim, sample) |>
  mutate(N=n()) |>
  filter(N >= 30) |>
  summarise(across(contains("pred"), ~rmse_vec(.x, truth=licePerFish_rtrt))) |>
  ungroup() |>
  left_join(mod_ids |>
              select(ens_id, sim_cols) |>
              unnest_wider(sim_cols))

# r: Site
r_null_site <- null_df |>
  group_by(sepaSiteNum) |>
  mutate(N=n()) |>
  filter(N >= 30) |>
  summarise(across(contains("pred"), ~cor(.x, licePerFish_rtrt, method="spearman", use="pairwise"))) |>
  ungroup()
r_cand_site <- candidate_df |>
  group_by(sepaSiteNum, s) |>
  mutate(N=n()) |>
  filter(N >= 30) |>
  summarise(across(contains("pred"), ~cor(.x, licePerFish_rtrt, method="spearman", use="pairwise"))) |>
  ungroup()
r_ens_site <- CV_df |>
  group_by(sepaSiteNum, ens_id, D, nSim, sample) |>
  mutate(N=n()) |>
  filter(N >= 30) |>
  summarise(across(contains("pred"), ~cor(.x, licePerFish_rtrt, method="spearman", use="pairwise"))) |>
  ungroup() |>
  left_join(mod_ids |>
              select(ens_id, sim_cols) |>
              unnest_wider(sim_cols))

# r: Date
r_null_date <- null_df |>
  group_by(date) |>
  mutate(N=n()) |>
  filter(N >= 30) |>
  summarise(across(contains("pred"), ~cor(.x, licePerFish_rtrt, method="spearman", use="pairwise"))) |>
  ungroup() |>
  mutate(null_pred=0)
r_cand_date <- candidate_df |>
  group_by(date, s) |>
  mutate(N=n()) |>
  filter(N >= 30) |>
  summarise(across(contains("pred"), ~cor(.x, licePerFish_rtrt, method="spearman", use="pairwise"))) |>
  ungroup()
r_ens_date <- CV_df |>
  group_by(date, ens_id, D, nSim, sample) |>
  mutate(N=n()) |>
  filter(N >= 30) |>
  summarise(across(contains("pred"), ~cor(.x, licePerFish_rtrt, method="spearman", use="pairwise"))) |>
  ungroup() |>
  left_join(mod_ids |>
              select(ens_id, sim_cols) |>
              unnest_wider(sim_cols))


# ROC-AUC: Site
ROC_AUC_null_site <- null_df |>
  mutate(lice_g05=factor(licePerFish_rtrt^4 >= 0.5)) |>
  group_by(sepaSiteNum) |>
  mutate(N=n()) |>
  filter(N >= 30 & n_distinct(lice_g05) > 1) |>
  summarise(across(contains("pred"), ~roc_auc_vec(.x, truth=lice_g05, event_level="second"))) |>
  ungroup()
ROC_AUC_cand_site <- candidate_df |>
  mutate(lice_g05=factor(licePerFish_rtrt^4 >= 0.5)) |>
  group_by(sepaSiteNum, s) |>
  mutate(N=n()) |>
  filter(N >= 30 & n_distinct(lice_g05) > 1) |>
  summarise(across(contains("pred"), ~roc_auc_vec(.x, truth=lice_g05, event_level="second"))) |>
  ungroup()
ROC_AUC_ens_site <- CV_df |>
  mutate(lice_g05=factor(licePerFish_rtrt^4 >= 0.5)) |>
  group_by(sepaSiteNum, ens_id, D, nSim, sample) |>
  mutate(N=n()) |>
  filter(N >= 30 & n_distinct(lice_g05) > 1) |>
  summarise(across(contains("pred"), ~roc_auc_vec(.x, truth=lice_g05, event_level="second"))) |>
  ungroup() |>
  left_join(mod_ids |>
              select(ens_id, sim_cols) |>
              unnest_wider(sim_cols))

# ROC-AUC: Date
ROC_AUC_null_date <- null_df |>
  mutate(lice_g05=factor(licePerFish_rtrt^4 >= 0.5)) |>
  group_by(date) |>
  mutate(N=n()) |>
  filter(N >= 30 & n_distinct(lice_g05) > 1) |>
  summarise(across(contains("pred"), ~roc_auc_vec(.x, truth=lice_g05, event_level="second"))) |>
  ungroup()
ROC_AUC_cand_date <- candidate_df |>
  mutate(lice_g05=factor(licePerFish_rtrt^4 >= 0.5)) |>
  group_by(date, s) |>
  mutate(N=n()) |>
  filter(N >= 30 & n_distinct(lice_g05) > 1) |>
  summarise(across(contains("pred"), ~roc_auc_vec(.x, truth=lice_g05, event_level="second"))) |>
  ungroup()
ROC_AUC_ens_date <- CV_df |>
  mutate(lice_g05=factor(licePerFish_rtrt^4 >= 0.5)) |>
  group_by(date, ens_id, D, nSim, sample) |>
  mutate(N=n()) |>
  filter(N >= 30 & n_distinct(lice_g05) > 1) |>
  summarise(across(contains("pred"), ~roc_auc_vec(.x, truth=lice_g05, event_level="second"))) |>
  ungroup() |>
  left_join(mod_ids |>
              select(ens_id, sim_cols) |>
              unnest_wider(sim_cols))


resampleEval_df <- bind_rows(
  RMSE_ens_date |>
    full_join(map(1:20,
                  ~RMSE_ens_date |> join_candidates_date(RMSE_cand_date, .x)) |>
                reduce(full_join, by=join_by(date, ens_id)),
              by=join_by(date, ens_id)) |>
    select(ens_id, D, nSim, sample, date, contains("pred")) |>
    pivot_longer(contains("pred")) |>
    drop_na() |>
    left_join(RMSE_null_date) |>
    mutate(metric="RMSE",
           type="date",
           group_num=as.numeric(as.factor(date))) |>
    select(-date),
  r_ens_date |>
    full_join(map(1:20,
                  ~r_ens_date |> join_candidates_date(r_cand_date, .x)) |>
                reduce(full_join, by=join_by(date, ens_id)),
              by=join_by(date, ens_id)) |>
    select(ens_id, D, nSim, sample, date, contains("pred")) |>
    pivot_longer(contains("pred")) |>
    drop_na() |>
    left_join(r_null_date) |>
    mutate(metric="r",
           type="date",
           group_num=as.numeric(as.factor(date))) |>
    select(-date),
  ROC_AUC_ens_date |>
    full_join(map(1:20,
                  ~ROC_AUC_ens_date |> join_candidates_date(ROC_AUC_cand_date, .x)) |>
                reduce(full_join, by=join_by(date, ens_id)),
              by=join_by(date, ens_id)) |>
    select(ens_id, D, nSim, sample, date, contains("pred")) |>
    pivot_longer(contains("pred")) |>
    drop_na() |>
    left_join(ROC_AUC_null_date) |>
    mutate(metric="ROC_AUC",
           type="date",
           group_num=as.numeric(as.factor(date))) |>
    select(-date),
  RMSE_ens_site |>
    full_join(map(1:20,
                  ~RMSE_ens_site |> join_candidates_siteNum(RMSE_cand_site, .x)) |>
                reduce(full_join, by=join_by(sepaSiteNum, ens_id)),
              by=join_by(sepaSiteNum, ens_id)) |>
    select(ens_id, D, nSim, sample, sepaSiteNum, contains("pred")) |>
    pivot_longer(contains("pred")) |>
    drop_na() |>
    left_join(RMSE_null_site) |>
    mutate(metric="RMSE",
           type="site",
           group_num=as.numeric(as.factor(sepaSiteNum))) |>
    select(-sepaSiteNum),
  r_ens_site |>
    full_join(map(1:20,
                  ~r_ens_site |> join_candidates_siteNum(r_cand_site, .x)) |>
                reduce(full_join, by=join_by(sepaSiteNum, ens_id)),
              by=join_by(sepaSiteNum, ens_id)) |>
    select(ens_id, D, nSim, sample, sepaSiteNum, contains("pred")) |>
    pivot_longer(contains("pred")) |>
    drop_na() |>
    left_join(r_null_site) |>
    mutate(metric="r",
           type="site",
           group_num=as.numeric(as.factor(sepaSiteNum))) |>
    select(-sepaSiteNum),
  ROC_AUC_ens_site |>
    full_join(map(1:20,
                  ~ROC_AUC_ens_site |> join_candidates_siteNum(ROC_AUC_cand_site, .x)) |>
                reduce(full_join, by=join_by(sepaSiteNum, ens_id)),
              by=join_by(sepaSiteNum, ens_id)) |>
    select(ens_id, D, nSim, sample, sepaSiteNum, contains("pred")) |>
    pivot_longer(contains("pred")) |>
    drop_na() |>
    left_join(ROC_AUC_null_site) |>
    mutate(metric="ROC_AUC",
           type="site",
           group_num=as.numeric(as.factor(sepaSiteNum))) |>
    select(-sepaSiteNum),
) |>
  # filter(name != "avg_pred") |>
  mutate(skill_range=if_else(metric %in% c("RMSE", "MAE"), 
                             0 - null_pred,
                             1 - null_pred),
         mod_m_null=value - null_pred,
         skill=mod_m_null / skill_range) |>
  arrange(ens_id, metric, type, group_num, name) |>
  group_by(ens_id, metric, type, group_num) |>
  mutate(ens_m_cand=last(value)-value,
         ens_m_cand_pct=ens_m_cand/value * 100,
         rank=if_else(metric %in% c("RMSE", "MAE"), 
                      min_rank(value), 
                      min_rank(desc(value)))) |>
  ungroup() |>
  mutate(nSim=factor(nSim, levels=levels(mod_ids$nSims)),
         ens_id=factor(ens_id, levels=levels(mod_ids$ens_id_ord))) |>
  mutate(candID=str_remove(name, "_pred")) |>
  left_join(mod_ids |> 
              select(ens_id, sim_cols) |>
              unnest(sim_cols) |>
              pivot_longer(starts_with("cand"), names_to="candID", values_to="simID") |>
              drop_na() |>
              mutate(simID=str_remove(simID, "IP_")),
            by=join_by(ens_id, candID)) |>
  mutate(modType=str_split_fixed(name, "_", 2)[,1]) |>
  mutate(name=if_else(grepl("cand", name), simID, name))

write_csv(resampleEval_df, "out/ensBlend_resample_performance_NEW_TEMP.csv")


# summarize and visualize -------------------------------------------------

resampleEval_df <- read_csv("out/ensBlend_resample_performance_NEW_TEMP.csv")
resampleEval_df |>
  group_by(modType, name, metric) |>
  summarise(mnRank=mean(rank, na.rm=T),
            prop1=mean(rank==1, na.rm=T)) |>
  arrange(mnRank) |>
  group_by(metric, modType) |>
  slice_head(n=1) |>
  print(n=50)

resampleEval_df |>
  group_by(type, nSim, D, modType, name, metric) |>
  summarise(mnRank=mean(rank, na.rm=T),
            prop1=mean(rank==1, na.rm=T)) |>
  mutate(mnTile=mnRank/(as.numeric(str_sub(nSim, 2, -1))+1)) |>
  ggplot(aes(mnTile, name, colour=type)) + 
  geom_point() + 
  facet_grid(metric~D*nSim)

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

ggpubr::ggarrange(p1, p2, nrow=1, common.legend=T)
ggsave("figs/pub_new/ens_nBetter_meanSkill_type.png", width=9, height=8)





resampleEval_df |>
  filter(nSim != "n20") |>
  filter(modType %in% c("ens", "cand")) |>
  # mean rank among weeks or farms
  group_by(metric, type, nSim, modType, name, sample) |> 
  summarise(mnSkill=mean(skill)) |> 
  arrange(nSim, type, sample, metric, mnSkill) |>
  group_by(nSim, type, sample, metric) |>
  mutate(percentile=percent_rank(mnSkill)) |>
  # filter(modType=="ens") |>
  filter(type=="site") |>
  mutate(nSim=factor(nSim, levels=paste0("n", c(3,5,10,15)), labels=c(3, 5, 10,15))) |>
  ggplot(aes(percentile, name)) + 
  ggdist::stat_histinterval(normalize="xy") + 
  facet_grid(metric~nSim) +
  theme_classic()

resampleEval_df |>
  filter(nSim != "n20") |>
  filter(type=="site") |>
  filter(modType %in% c("ens", "cand")) |>
  # mean rank among weeks or farms
  group_by(metric, type, nSim, modType, name, sample) |> 
  summarise(mnSkill=mean(skill)) |> 
  arrange(nSim, type, sample, metric, modType) |>
  group_by(nSim, type, sample, metric) |>
  summarise(nBetterThan=sum(last(mnSkill) > mnSkill, na.rm=T)) |>
  mutate(prEnsBetter=nBetterThan/as.numeric(str_sub(nSim, 2, -1)),
         nSim=factor(nSim, levels=paste0("n", c(3,5,10,15)), labels=c(3, 5, 10,15))) |>
  ggplot(aes(nBetterThan)) +
  geom_bar(colour="grey30", linewidth=0.1) +
  scale_x_continuous(breaks=0:15) +
  facet_grid(metric~nSim, scales="free") +
  theme_classic()
  
resampleEval_df |>
  filter(nSim != "n20") |>
  group_by(metric, type, nSim, modType, name, sample) |> 
  summarise(mnSkill=mean(skill)) |> 
  arrange(nSim, type, sample, metric, mnSkill) |>
  group_by(nSim, type, sample, metric) |>
  mutate(percentile=percent_rank(mnSkill),
         nSim=factor(nSim, levels=paste0("n", c(3,5,10,15)), labels=c(3, 5, 10,15))) |>
  ggplot(aes(percentile, as.character(sample), fill=modType)) + 
  geom_raster() + 
  scale_fill_manual(values=c("cand"="#a6cee3", "ens"="#33a02c", "avg"="#b2df8a")) +
  facet_grid(type*nSim~metric, scales="free_y") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())



resampleEval_df |>
  filter(nSim != "n20") |>
  group_by(metric, type, nSim, modType, name, sample) |> 
  summarise(mnSkill=mean(skill)) |> 
  arrange(nSim, type, sample, metric, mnSkill) |>
  group_by(nSim, type, sample, metric) |>
  mutate(percentile=percent_rank(mnSkill),
         nSim=factor(nSim, levels=paste0("n", c(3,5,10,15)), labels=c(3, 5, 10,15))) |>
  filter(modType != "cand") |>
  ggplot(aes(percentile, fill=modType)) + 
  geom_bar(position="dodge") + 
  # scale_fill_manual(values=c("cand"="#a6cee3", "ens"="#33a02c", "avg"="#b2df8a")) +
  facet_grid(type*nSim~metric, scales="free_y") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())




CV_df |> 
  # filter(sepaSiteNum < 5) |>
  arrange(sepaSiteNum, date) |>
  group_by(ens_id, sepaSiteNum) |>
  mutate(dayDiff=date - lag(date),
         gap=dayDiff > 28,
         series=dplyr::consecutive_id(gap)) |>
  ungroup() |>
  ggplot(aes(date)) + 
  geom_point(data=CV_df |> filter(ens_id==first(ens_id)),
             aes(y=licePerFish_rtrt), shape=1, size=0.25, colour="grey30") +
  geom_line(aes(y=ens_pred, group=paste(ens_id, series), colour=D), alpha=0.2) + 
  scale_colour_brewer(type="qual", palette=2) +
  facet_wrap(~sepaSiteNum)




bestOfEach_df <- resampleEval_df |> 
  filter(nSim != "n20") |>
  filter(modType != "avg") |> 
  # mean rank among weeks or farms
  group_by(metric, type, nSim, modType, name, sample) |> 
  summarise(mnRank=mean(rank), mdRank=median(rank),
            mnVal=mean(value), mnSkill=mean(skill)) |> 
  # select ensBlend, ensMean, and best constituent
  group_by(metric, type, nSim, modType, sample) |> 
  slice_max(mnSkill, with_ties=F) |>
  group_by(metric, type, nSim, sample) |>
  arrange(modType) |>
  mutate(nSim=factor(nSim, levels=paste0("n", c(3,5,10,15))),
         bestMod=if_else(first(mnSkill) > last(mnSkill), "cand", "ens"),
         bestMod=factor(bestMod, 
                        levels=c("cand", "ens"),
                        labels=c("Opt['Param']", "Ens['Blend']")))

p <- bestOfEach_df |> 
  arrange(metric, type, nSim, sample, name) |> 
  mutate(bestCand=str_sub(last(name), -2, -1),
         sim07_included=paste("3D.7,3", 
                              if_else(bestCand=="07" | bestCand=="03",
                                      "in resample",
                                      "NOT in resample"))) |> 
  group_by(metric, type, nSim, sample) |>
  slice_head(n=1) |>
  ungroup() |>
  mutate(metric=factor(metric, levels=c("RMSE", "r", "ROC_AUC"),
                       labels=c("RMSE", "Spearmans~~rho", "AUC['ROC']"))) |>
  ggplot(aes(nSim, fill=bestMod)) +
  geom_bar(position="fill") + 
  scale_fill_manual("Best model", values=c("#9FB6CC", "#ca0020"), 
                    labels=scales::label_parse()) +
  scale_y_continuous("Percentage of resamples", 
                     breaks=c(0, 0.5, 1),
                     labels=scales::label_percent()) +
  xlab("Number of parameterizations per resample") +
  facet_grid(metric~type, labeller=labeller(metric=label_parsed,
                                                      sim07_included=label_wrap_gen(12))) +
  theme(panel.grid.major.y=element_line(colour="grey80"),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())
p
ggsave("figs/pub/ensBlend_eval_PrEnsBest.png", p, width=6, height=8)


bestOfEach_df |>
  arrange(metric, type, nSim, sample, name) |> 
  mutate(bestCand=str_sub(last(name), -2, -1),
         sim04_included=paste("3D.7", 
                              if_else(bestCand=="07",
                                      "in resample",
                                      "NOT in resample"))) |> 
  arrange(metric, type, nSim, sample, modType) |> 
  group_by(metric, type, nSim, sample, sim04_included) |>
  summarise(bestMod=first(bestMod), 
            dSkill=(last(mnSkill)-first(mnSkill))/first(mnSkill)*100) |> 
  ggplot(aes(dSkill, nSim, fill=sim04_included)) + 
  ggdist::stat_halfeye(alpha=0.5) + 
  facet_grid(metric~.)



resampleEval_df |>
  mutate(candID=str_remove(name, "_pred")) |>
  left_join(mod_ids |> 
              select(ens_id, sim_cols) |>
              unnest(sim_cols) |>
              pivot_longer(starts_with("cand"), names_to="candID", values_to="simID") |>
              drop_na(),
            by=join_by(ens_id, candID)) |>
  filter(nSim=="n3") |>
  glimpse()


bestCandidates_df <- resampleEval_df |> 
  inner_join(bestOfEach_df |> filter(modType=="cand") |> select(metric, type, nSim, name, sample),
             by=join_by(metric, type, nSim, name, sample)) |>
  mutate(bestMod=case_when(metric %in% c("MAE", "RMSE") & ens_m_cand < 0 ~ "ens",
                           metric %in% c("MAE", "RMSE") & ens_m_cand > 0 ~ "cand",
                           !metric %in% c("MAE", "RMSE") & ens_m_cand > 0 ~ "ens",
                           !metric %in% c("MAE", "RMSE") & ens_m_cand < 0 ~ "cand"))
bestCandidates_df |>
  group_by(ens_id, metric, type, nSim, sample) |>
  summarise(mn_ens_m_cand=mean(ens_m_cand, na.rm=T)) |>
  ungroup() |>
  mutate(bestMod=case_when(metric %in% c("MAE", "RMSE") & mn_ens_m_cand < 0 ~ "ens",
                           metric %in% c("MAE", "RMSE") & mn_ens_m_cand > 0 ~ "cand",
                           !metric %in% c("MAE", "RMSE") & mn_ens_m_cand > 0 ~ "ens",
                           !metric %in% c("MAE", "RMSE") & mn_ens_m_cand < 0 ~ "cand")) |>
  ggplot(aes(mn_ens_m_cand, type, fill=nSim)) +
    geom_vline(xintercept=0) +
    geom_boxplot() + 
    facet_grid(metric~.)

bestCandidates_df |>
  group_by(ens_id, metric, type, nSim, sample) |>
  summarise(mn_ens_m_cand=mean(ens_m_cand, na.rm=T)) |>
  ungroup() |>
  mutate(bestMod=case_when(metric %in% c("MAE", "RMSE") & mn_ens_m_cand < 0 ~ "ens",
                           metric %in% c("MAE", "RMSE") & mn_ens_m_cand > 0 ~ "cand",
                           !metric %in% c("MAE", "RMSE") & mn_ens_m_cand > 0 ~ "ens",
                           !metric %in% c("MAE", "RMSE") & mn_ens_m_cand < 0 ~ "cand")) |>
  mutate(ens_improvement=if_else(metric %in% c("MAE", "RMSE"), -mn_ens_m_cand, mn_ens_m_cand)) |>
  ggplot(aes(ens_improvement, nSim, fill=nSim)) +
  geom_vline(xintercept=0) +
  ggdist::stat_halfeye() +
  scale_fill_viridis_d("Number of\nconstituent\nparameterizations", 
                         option="mako", end=0.8, direction=-1) +
  xlab("Ensemble improvement (+ = better)") +
  facet_grid(metric~type)

bestCandidates_df |>
  group_by(ens_id, metric, type, nSim, sample) |>
  summarise(mn_ens_m_cand=mean(ens_m_cand, na.rm=T)) |>
  ungroup() |>
  mutate(bestMod=case_when(metric %in% c("MAE", "RMSE") & mn_ens_m_cand < 0 ~ "ens",
                           metric %in% c("MAE", "RMSE") & mn_ens_m_cand > 0 ~ "cand",
                           !metric %in% c("MAE", "RMSE") & mn_ens_m_cand > 0 ~ "ens",
                           !metric %in% c("MAE", "RMSE") & mn_ens_m_cand < 0 ~ "cand")) |>
  ggplot(aes(nSim, fill=bestMod)) + 
  geom_bar(position="fill") +
  scale_fill_manual("Best model", values=c("#9FB6CC", "#ca0020"), 
                    labels=scales::label_parse()) +
  scale_y_continuous("Percentage of resamples", 
                     breaks=c(0, 0.5, 1),
                     labels=scales::label_percent()) +
  xlab("Number of parameterizations per resample") +
  facet_grid(.~metric) +
  theme(panel.grid.major.y=element_line(colour="grey80"),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())

p <- bestOfEach_df |>
  filter(metric=="ROC_AUC") |>
  ggplot(aes(modType, mnVal, group=paste(metric, type, nSim, sample), colour=bestMod)) + 
  geom_line(alpha=0.5) +
  scale_colour_manual("Best model", values=c("#9FB6CC", "#ca0020"), 
                      labels=scales::label_parse()) +
  facet_grid(type~nSim)
ggsave("figs/temp_p.png", p, width=10, height=10)  



cand_v_ens_df <- resampleEval_df |>
  mutate(modType=str_split_fixed(name, "_", 2)[,1]) |> 
  arrange(rank) |>
  group_by(ens_id, nSim, sample, metric, type, modType, group_num) |>
  slice_head(n=1) |>
  ungroup() 
  
cand_v_ens_df |>
  filter(modType=="cand") |>
  ggplot(aes(ens_m_cand, colour=nSim)) + 
  geom_vline(xintercept=0) +
  geom_density() +
  facet_grid(type~metric, scales="free")

prop_bar_df |>
  filter(metric=="RMSE") |>
  ggplot(aes(nSim, fill=modType)) + 
  geom_bar(position="fill") +
  scale_fill_manual("Best model", values=c("#9FB6CC", "#ca0020"), 
                    labels=scales::label_parse()) +
  scale_y_continuous("Percentage of resamples", 
                     breaks=c(0, 0.5, 1),
                     labels=scales::label_percent()) +
  xlab("Number of parameterizations per resample") +
  facet_grid(.~type) +
  theme(panel.grid.major.y=element_line(colour="grey80"),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())

lab_expressions <- c(expression(Ens['Blend']),
                     expression(Ens['Avg']),
                     expression(Constituent))
lab_expressions_07 <- c(expression(Ens['Blend']),
                        expression(Ens['Avg']),
                        expression('3D.7'),
                        expression(Constituent))
point_df <- resampleEval_df |> 
  filter(nSim != "n20") |>
  # mean rank among weeks or farms
  group_by(metric, type, nSim, modType, name, sample) |> 
  summarise(mnRank=mean(rank), mdRank=median(rank), value=mean(value), skill=mean(skill)) |>
  ungroup() |>
  mutate(modType=factor(modType, 
                        levels=c("ens", "avg", "cand"),
                        labels=c("Ens['Blend']", "Ens['Avg']", "Constituent")),
         type=factor(type, levels=c("site", "date"),
                     labels=paste("Mean among", c("farms", "weeks"))),
         sample=factor(sample, levels=1:100),
         nSim=factor(nSim, levels=paste0("n", c(3, 5, 10, 15, 20)),
                     labels=paste0("n: ", c(3, 5, 10, 15, 20)))) |>
  arrange(metric, type, nSim, sample, skill, desc(name)) |>
  group_by(metric, type, nSim, sample) |> 
  # mutate(bestMod=if_else(metric=="RMSE", first(modType), last(modType)),
  #        bestMod=factor(bestMod, levels=levels(modType))) |>
  mutate(bestMod=last(name),
         bestMod=case_when(bestMod=="ens_pred" ~ "Ens['Blend']",
                           bestMod=="avg_pred" ~ "Ens['Avg']",
                           bestMod=="sim_07" ~ "3D.7",
                           .default="Other"),
         bestMod=factor(bestMod, levels=c("Ens['Blend']", "Ens['Avg']", "3D.7", "Other"))) |>
  ungroup()
 


p <- point_df |> 
  filter(metric=="RMSE") |>
  ggplot(aes(value, sample)) +
  geom_line(data=point_df |> filter(metric=="RMSE") |> 
              group_by(metric, type, nSim, sample, modType) |>
              slice_max(skill) |>
              group_by(metric, type, nSim, sample) |> 
              slice_max(skill, n=2) |> ungroup(),
            aes(colour=bestMod)) +
  scale_colour_manual("Best model", values=c("#ca0020", "#3F6B99", "#9FB6CC"),
                      labels=lab_expressions_07) +
  ggnewscale::new_scale_colour() +
  geom_point(aes(colour=modType, shape=modType), size=2, alpha=0.75) +
  geom_point(data=point_df |> filter(metric=="RMSE") |> 
               group_by(metric, type, nSim, sample) |>
               slice_max(skill) |> ungroup(),
             aes(colour=modType, shape=modType), size=3) +
  scale_colour_manual("Model type", values=c("#ca0020", "#ca0020", "#9FB6CC"),
                      labels=lab_expressions) +
  scale_shape_manual("Model type", values=c("|", "X", "o"),
                     labels=lab_expressions) +
  labs(x="RMSE (mean)", y="Resample from 20 parameterizations") +
  facet_grid(nSim~type, scales="free_y", space="free_y") +
  theme(panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())
ggsave(glue("figs/pub_new/ensBlend_eval_RMSE.png"), p, width=6, height=10.5, dpi=300)


p <- point_df |> 
  filter(metric=="r") |>
  ggplot(aes(value, sample)) +
  geom_line(data=point_df |> filter(metric=="r") |> 
              group_by(metric, type, nSim, sample, modType) |>
              slice_max(value) |> ungroup(),
            aes(colour=bestMod)) +
  scale_colour_manual("Best model", values=c("#ca0020", "#ca0020", "#3F6B99", "#9FB6CC"),
                      labels=lab_expressions_07) +
  ggnewscale::new_scale_colour() +
  geom_point(aes(colour=modType, shape=modType), size=2, alpha=0.75) +
  geom_point(data=point_df |> filter(metric=="r") |> 
               group_by(metric, type, nSim, sample) |>
               slice_max(value) |> ungroup(),
             aes(colour=modType, shape=modType), size=3) +
  scale_colour_manual("Model type", values=c("#ca0020", "#ca0020", "#9FB6CC"),
                      labels=lab_expressions) +
  scale_shape_manual("Model type", values=c("|", "X", "o"),
                     labels=lab_expressions) +
  labs(x=expression(Spearmans~rho~~'(mean)'), y="Resample from 20 parameterizations") +
  facet_grid(nSim~type, scales="free_y", space="free_y") +
  theme(panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())
ggsave(glue("figs/pub_new/ensBlend_eval_rho.png"), p, width=6, height=10.5, dpi=300)


p <- point_df |> 
  filter(metric=="ROC_AUC") |>
  ggplot(aes(value, sample)) +
  geom_line(data=point_df |> filter(metric=="ROC_AUC") |> 
              group_by(metric, type, nSim, sample, modType) |>
              slice_max(value) |> ungroup() |>
              group_by(metric, type, nSim, sample) |> 
              slice_max(value, n=2),
            aes(colour=bestMod)) +
  scale_colour_manual("Best model", values=c("#ca0020", "#ca0020", "#3F6B99", "#9FB6CC"),
                      labels=lab_expressions_07) +
  ggnewscale::new_scale_colour() +
  geom_point(aes(colour=modType, shape=modType), size=2, alpha=0.75) +
  geom_point(data=point_df |> filter(metric=="ROC_AUC") |> 
               group_by(metric, type, nSim, sample) |>
               slice_max(value) |> ungroup(),
             aes(colour=modType, shape=modType), size=3) +
  scale_colour_manual("Model type", values=c("#ca0020", "#ca0020", "#9FB6CC"),
                      labels=lab_expressions) +
  scale_shape_manual("Model type", values=c("|", "X", "o"),
                     labels=lab_expressions) +
  scale_x_continuous(limits=c(0.69, 0.825), breaks=c(0.7, 0.75, 0.8)) +
  labs(x=expression(ROC['AUC']~'(median)'), y="Resample from 20 parameterizations") +
  facet_grid(nSim~type, scales="free_y", space="free_y") +
  theme(panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())
ggsave(glue("figs/pub_new/ensBlend_eval_ROC-AUC.png"), p, width=6, height=10.5, dpi=300)


test_df <- point_df |> 
  filter(nSim != "n: 20") |>
  group_by(metric, type, nSim, sample) |>
  slice_head(n=1) |>
  mutate(ensBest=bestMod=="Ens['Blend']") |>
  droplevels()
library(brms) 
out <- brm(ensBest ~ nSim * type * metric, 
           data=test_df, cores=4,
           family=bernoulli())
pred_df <- expand_grid(metric=unique(test_df$metric),
                       nSim=unique(test_df$nSim),
                       type=unique(test_df$type))
post_pred <- posterior_epred(out, pred_df)
pred_df <- pred_df |>
  mutate(post_mn=apply(post_pred, 2, mean),
         post_md=apply(post_pred, 2, median),
         post_lo=apply(post_pred, 2, quantile, probs=0.025),
         post_hi=apply(post_pred, 2, quantile, probs=0.975))
ggplot(pred_df, aes(nSim, post_mn, ymin=post_lo, ymax=post_hi, colour=nSim)) + 
  geom_pointrange(position=position_dodge(width=0.5)) +
  scale_colour_viridis_d("Number of\nconstituent\nparameterizations", 
                         option="mako", end=0.8, direction=-1) +
  scale_y_continuous("Pr(Ensemble outperforms all constituents)", limits=c(0,1),
                     breaks=c(0, 0.5, 1),
                     labels=scales::label_percent()) +
  xlab("Number of parameterizations per resample") +
  facet_grid(type~metric) +
  theme(panel.grid.major.y=element_line(colour="grey80"),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())




