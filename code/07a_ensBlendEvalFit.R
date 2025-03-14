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
ensFull_df <- read_csv("out/valid_df_2021-2024.csv") |>
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
ensFull_LatLon <- read_csv("out/valid_df_2021-2024.csv") |>
  select(rowNum, date, CV_k, sepaSite, sepaSiteNum, licePerFish_rtrt, starts_with("sim")) |>
  select(-contains("avg")) |>
  mutate(across(starts_with("sim_"), ~.x - mean(.x), .names="c_{.col}")) |>
  left_join(site_i) |>
  select(-sepaSite) |>
  arrange(rowNum)



# cross-validation --------------------------------------------------------

set.seed(66)
mods <- expand_grid(mod=c(paste0("sLonLatD", 3:4)),
                    nSims=c(paste0("n", c(2, 5, 10))),
                    sample=1:10)
mods <- mods |> 
  filter(mod=="sLonLatD3") |>
  rowwise() |>
  mutate(sim=list(sample(1:20, as.numeric(str_sub(nSims, 2, -1))))) |>
  ungroup() |>
  select(nSims, sample, sim) |>
  full_join(mods)
saveRDS(mods, "out/ensembles/ensBlend_EvalModSpecs.rds")

CV_ensBlend <- CV_ensAvg <- vector("list", length(folds))

for(k in (seq_along(folds))) {
  CV_k_ls <- CV_k_ensAvg_ls <- vector("list", nrow(mods))
  
  for(i in rev(1:nrow(mods))) {
    # EnsBlend setup
    if(grepl("sLonLat", mods$mod[i])) {
      recipe_i <- make_spline_recipe(ensFull_LatLon, 
                                     as.numeric(str_split_fixed(mods$mod[i], "D", 2)[2]), 
                                     sim_i$sim[mods$sim[[i]]])
      train_df <- bake(recipe_i, ensFull_LatLon |> filter(CV_k != folds[k]))
      test_df <- bake(recipe_i, ensFull_LatLon |> filter(CV_k == folds[k]))
      dat_rstan <- make_data_rstan_sLonLat(train_df)
      pars <- c("b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu",
                paste0("b_s_", c("easting", "northing", "easting_x_northing")))
    } else {
      train_df <- ensFull_df |> filter(CV_k != folds[k]) |>
        select(rowNum, sepaSite, sepaSiteNum, CV_k, date, licePerFish_rtrt,
               all_of(sim_i$sim[mods$sim[[i]]]), all_of(paste0("c_", sim_i$sim[mods$sim[[i]]])))
      test_df <- ensFull_df |> filter(CV_k == folds[k]) |>
        select(rowNum, sepaSite, sepaSiteNum, CV_k, date, licePerFish_rtrt,
               all_of(sim_i$sim[mods$sim[[i]]]), all_of(paste0("c_", sim_i$sim[mods$sim[[i]]])))
      dat_rstan <- make_data_rstan(train_df, R2D2_sd=F)
      pars <- c("b_p", "b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu",
                "b_p_uc", "r_grp", "sd_grp")
    }
    if(grepl("pRE_huRE", mods$mod[i])) {
      pars <- c(pars, "r_grp_hu", "sd_grp_hu") 
    }
    # fit EnsBlend
    fname <- glue("out/ensembles/ensBlend_{mods$mod[i]}_{mods$nSims[i]}-{mods$sample[i]}_CV-{folds[k]}")
    if(file.exists(glue("{fname}_stanfit.rds"))) {
      cat("File exists:", fname, "\n")
      out_ensBlend <- readRDS(glue("{fname}_stanfit.rds"))
    } else {
      stanMod <- str_split_fixed(mods$mod[i], "D", 2)[1]
      out_ensBlend <- stan(file=glue("code/stan/ensemble_mixture_model_{stanMod}.stan"),
                           model_name=mods$mod[i], data=dat_rstan,
                           chains=30, cores=30, iter=3100, warmup=3000,
                           # chains=3, cores=3, iter=3000, warmup=2000,
                           control=list(adapt_delta=0.95, max_treedepth=20),
                           pars=pars)
      saveRDS(out_ensBlend, glue("{fname}_stanfit.rds"))
      saveRDS(dat_rstan, glue("{fname}_standata.rds")) 
    }
    if(grepl("sLonLat", mods$mod[i])) {
      CV_k_ls[[i]] <- out_ensBlend |>
        make_predictions_ensBlend_sLonLat(test_df, mode="point_epred") |>
        colMeans() |>
        as_tibble() |>
        set_names(paste0("IP_", mods$mod[i], "_", mods$nSims[i], "_", mods$sample[i]))
    } else {
      CV_k_ls[[i]] <- out_ensBlend |>
        make_predictions_ensBlend(test_df, re=T, re_hu=grepl("huRE", mods$mod[i])) |>
        colMeans() |>
        as_tibble() |>
        set_names(paste0("IP_", mods$mod[i], "_", mods$nSims[i], "_", mods$sample[i]))
    }
    if(mods$mod[i]=="sLonLatD3") {
      # fit using avg AEIP
      dat_avg_df <- ensFull_df |>
        select(rowNum, sepaSite, sepaSiteNum, CV_k, date, licePerFish_rtrt,
               all_of(sim_i$sim[mods$sim[[i]]])) |>
        rowwise() |>
        mutate(sim_avg=mean(c_across(starts_with("sim")))) |>
        ungroup() |>
        mutate(c_sim_avg=c(scale(sim_avg))) |>
        select(-matches("sim_[0-9]"))
      train_df <- dat_avg_df |> filter(CV_k != folds[k])
      test_df <- dat_avg_df |> filter(CV_k == folds[k])
      fname_avg <- glue("out/ensembles/ensAvg_{mods$nSims[i]}-{mods$sample[i]}_CV-{folds[k]}")
      if(file.exists(fname_avg)) {
        out_sim <- readRDS(fname_avg)
      } else {
        dat_rstan <- train_df |> make_data_rstan()
        out_sim <- stan(file="code/stan/candidate_model.stan",
                        model_name=glue("avg-{folds[k]}"), data=dat_rstan,
                        chains=6, cores=6,iter=3000, warmup=2500,
                        pars=c("b_b0", "b_IP", "Intercept_hu", "b_hu", "sigma"))
        saveRDS(out_sim, fname_avg)
      }
      CV_k_ensAvg_ls[[i]] <- make_predictions_candidate(out_sim, test_df, "avg") |>
        colMeans() |>
        as_tibble() |>
        set_names(paste0("IP_avg", "_", mods$nSims[i], "_", mods$sample[i]))
    }
  }
  CV_ensBlend[[k]] <- bind_cols(test_df |> select(rowNum),
                                reduce(CV_k_ls, bind_cols))
  CV_ensAvg[[k]] <- bind_cols(test_df |> select(rowNum),
                                reduce(CV_k_ensAvg_ls, bind_cols))
}
reduce(CV_ensBlend, bind_rows) |>
  write_csv("out/ensembles/CV_ensBlend_EvalCV.csv")
reduce(CV_ensAvg, bind_rows) |>
  write_csv("out/ensembles/CV_ensAvg_EvalCV.csv")



# compile predictions -----------------------------------------------------

mod_ids <- readRDS("out/ensembles/ensBlend_EvalModSpecs.rds") |> 
  mutate(ens_id=paste("IP", str_sub(mod, -2, -1), nSims, sample, sep="_"),
         sim_cols=map(sim, 
                      ~tibble(s=paste0("IP_sim_", str_pad(.x, 2, "left", "0"))) |>
                        mutate(cand_id=paste0("cand", row_number())) |> 
                        pivot_wider(names_from=cand_id, values_from=s))) |>
  bind_rows(tibble(ens_id=paste0("IP_D", 3:10, "_n20_1"),
                   mod=paste0("sLonLatD", 3:10),
                   nSims="n20",
                   sample=1)) |>
  mutate(nSims=factor(nSims, levels=c("n2", "n5", "n10", "n20"))) |>
  arrange(nSims, sample, mod) |>
  mutate(ens_id_ord=factor(ens_id, levels=unique(ens_id)))
candidate_df <- read_csv("out/candidates/CV_candidate_predictions.csv") |>
  pivot_longer(starts_with("IP_"), names_to="s", values_to="cand_pred")
avg_df <- read_csv("out/ensembles/CV_ensAvg_EvalCV.csv") |>
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
  mutate(ens_id=paste0(str_remove(ens_id, "sLonLat"), "_1")) |>
  bind_cols(map_dfc(1:20, 
                    ~tibble(x=paste0("IP_sim_", str_pad(.x, 2, "left", "0"))) |>
                      set_names(paste0("cand", .x))))
CV_df <- read_csv("out/ensembles/CV_ensBlend_EvalCV.csv") |>
  rename_with(~str_remove(.x, "sLonLat")) |>
  inner_join(ensFull_LatLon |> 
              select(rowNum, licePerFish_rtrt, date, sepaSiteNum)) |>
  mutate(lice_g05=factor(licePerFish_rtrt^4 > 0.5)) |> 
  pivot_longer(starts_with("IP_D"), names_to="ens_id", values_to="ens_pred") |>
  left_join(mod_ids |>
              select(ens_id, sim_cols) |>
              unnest_wider(sim_cols)) |>
  bind_rows(ens20_df)

CV_df <- CV_df |>
  full_join(map(1:20, ~CV_df |> join_candidates(candidate_df, .x)) |> 
              reduce(full_join, by=join_by(rowNum, ens_id)), 
            by=join_by(rowNum, ens_id)) |>
  separate_wider_delim(ens_id, "_", names=c("x", "D", "nSim", "sample"), cols_remove=F) |>
  select(-x) |>
  full_join(avg_df, by=join_by(rowNum, nSim, sample))

library(yardstick)
RMSE_site <- CV_df |>
  group_by(sepaSiteNum, ens_id, D, nSim, sample) |>
  mutate(N=n()) |>
  filter(N >= 30) |>
  summarise(across(contains("pred"), ~rmse_vec(.x, truth=licePerFish_rtrt))) |>
  group_by(ens_id, D, nSim, sample) |>
  summarise(across(contains("pred"), median)) |>
  pivot_longer(contains("pred")) |>
  arrange(ens_id, name) |>
  group_by(ens_id) |>
  mutate(ens_m_cand=last(value)-value,
         ens_m_cand_pct=ens_m_cand/value * 100,
         rank=min_rank(value)) |>
  ungroup()

RMSE_date <- CV_df |>
  group_by(date, ens_id, D, nSim, sample) |>
  mutate(N=n()) |>
  filter(N >= 30) |>
  summarise(across(contains("pred"), ~rmse_vec(.x, truth=licePerFish_rtrt))) |>
  group_by(ens_id, D, nSim, sample) |>
  summarise(across(contains("pred"), median)) |>
  pivot_longer(contains("pred")) |>
  arrange(ens_id, name) |>
  group_by(ens_id) |>
  mutate(ens_m_cand=last(value)-value,
         ens_m_cand_pct=ens_m_cand/value * 100,
         rank=min_rank(value)) |>
  ungroup()

metric_df <- bind_rows(
  RMSE_site |> mutate(metric="RMSE", type="site"),
  RMSE_date |> mutate(metric="RMSE", type="date")
) |>
  mutate(nSim=factor(nSim, levels=levels(mod_ids$nSims)),
         ens_id=factor(ens_id, levels=levels(mod_ids$ens_id_ord)))



# summarize and visualize -------------------------------------------------

metric_df |> 
  filter(D=="D3") |>
  group_by(name) |>
  summarise(mnRank=mean(rank, na.rm=T),
            prop1=mean(rank==1, na.rm=T)) |> 
  arrange(mnRank)

metric_df |> 
  filter(D=="D3") |>
  group_by(nSim, name) |>
  summarise(mnRank=mean(rank, na.rm=T),
            prop1=mean(rank==1, na.rm=T)) |> 
  ungroup() |>
  arrange(nSim, mnRank) |>
  drop_na() |>
  print(n=50)

metric_df |> 
  filter(D=="D3") |> 
  mutate(model_type=str_sub(name, 1, 3)) |> 
  group_by(D, nSim, sample, metric, type, model_type) |> 
  slice_min(rank) |> 
  group_by(D, nSim, sample, metric, type) |> 
  arrange(type, name) |> 
  mutate(cand_rank=nth(rank, 2)) |> 
  ungroup() |> 
  filter(model_type != "can") |> 
  group_by(D, nSim, metric, model_type) |> 
  summarise(Pr_ens_better=mean(rank < cand_rank)) |>
  ungroup() |>
  arrange(model_type, nSim, metric)

temp_df <- metric_df |>
  filter(D=="D3") |>
  group_by(ens_id, type, metric) |>
  arrange(rank) |>
  mutate(bestMod=case_when(first(name)=="ens_pred" ~ "Ens['Blend']",
                           first(name)=="avg_pred" ~ "Ens['Mean']",
                           .default="Variant")) |>
  arrange(name) |>
  mutate(ensPct=100-last(percent_rank(rank))*100) |>
  ungroup() |>
  mutate(modType=case_when(name=="ens_pred" ~ "Ens['Blend']",
                           name=="avg_pred" ~ "Ens['Mean']",
                           .default="Variant"),
         sample=factor(sample, levels=10:1),
         type=factor(type, levels=c("site", "date"),
                     labels=paste("Median among", c("farms", "weeks"))),
         nSim=factor(nSim, levels=paste0("n", c(2, 5, 10, 20)),
                     labels=paste0("n: ", c(2, 5, 10, 20))))
lab_expressions <- c(expression(Ens['Blend']),
                     expression(Ens['Mean']),
                     expression(Variant))
p <- temp_df |>
  ggplot(aes(value, sample)) + 
  geom_line(data=temp_df |> filter(!grepl("ens|avg", name)) |> arrange(value),
            aes(colour=bestMod)) +
  scale_colour_manual("Best model", values=c("#ca0020", "#3F6B99", "#9FB6CC"),
                      labels=lab_expressions) +
  ggnewscale::new_scale_colour() +
  geom_point(aes(shape=modType, size=modType, colour=modType), alpha=0.75) +
  scale_colour_manual("Model type", values=c("#ca0020", "#3F6B99", "#3F6B99"),
                      labels=lab_expressions) +
  scale_shape_manual("Model type", values=c(1, 5, 1),
                     labels=lab_expressions) +
  scale_size_manual("Model type", values=c(2.5, 2, 1),
                    labels=lab_expressions) +
  facet_grid(nSim~type, scales="free_y", space="free_y") +
  labs(x="RMSE", y="Sample from pool of 20 variants") +
  theme(panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank())
ggsave("figs/pub/ensBlend_eval_RMSE.png", p, width=6, height=10.5)
