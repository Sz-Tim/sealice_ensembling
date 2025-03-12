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

CV_ensBlend <- vector("list", length(folds))

for(k in (seq_along(folds))) {
  CV_k_ls <- vector("list", nrow(mods))
  
  for(i in rev(1:nrow(mods))) {
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
    
    fname <- glue("out/ensembles/ensBlend_{mods$mod[i]}_{mods$nSims[i]}-{mods$sample[i]}_CV-{folds[k]}")
    if(file.exists(glue("{fname}_stanfit.rds"))) {
      cat("File exists:", fname, "\n")
      out_ensBlend <- readRDS(glue("{fname}_stanfit.rds"))
    } else {
      stanMod <- str_split_fixed(mods$mod[i], "D", 2)[1]
      out_ensBlend <- stan(file=glue("code/stan/ensemble_mixture_model_{stanMod}.stan"),
                           model_name=mods$mod[i], data=dat_rstan,
                           # chains=30, cores=30, iter=2250, warmup=2000,
                           chains=3, cores=3, iter=3000, warmup=2000,
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
  }
  CV_ensBlend[[k]] <- bind_cols(test_df |> select(rowNum),
                                reduce(CV_k_ls, bind_cols))
}
reduce(CV_ensBlend, bind_rows) |>
  write_csv("out/ensembles/CV_ensBlend_EvalCV.csv")



# compile predictions -----------------------------------------------------

mod_ids <- readRDS("out/ensembles/ensBlend_EvalModSpecs.rds") |> 
  mutate(ens_id=paste("IP", str_sub(mod, -2, -1), nSims, sample, sep="_"),
         nSims=factor(nSims, levels=c("n2", "n5", "n10", "n20")),
         sim_cols=map(sim, 
                      ~tibble(s=paste0("IP_sim_", str_pad(.x, 2, "left", "0"))) |>
                        mutate(cand_id=paste0("cand", row_number())) |> 
                        pivot_wider(names_from=cand_id, values_from=s))) |>
  arrange(nSims, sample, mod) |>
  mutate(ens_id_ord=factor(ens_id, levels=unique(ens_id)))
candidate_df <- read_csv("out/candidates/CV_candidate_predictions.csv") |>
  pivot_longer(starts_with("IP_"), names_to="s", values_to="cand_pred")
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
            by=join_by(rowNum, ens_id))

library(yardstick)
RMSE_site <- CV_df |>
  group_by(sepaSiteNum, ens_id) |>
  mutate(N=n()) |>
  filter(N >= 30) |>
  summarise(across(contains("pred"), ~rmse_vec(.x, truth=licePerFish_rtrt))) |>
  group_by(ens_id) |>
  summarise(across(contains("pred"), mean)) |>
  pivot_longer(contains("pred")) |>
  arrange(ens_id, name) |>
  group_by(ens_id) |>
  mutate(ens_m_cand=last(value)-value,
         ens_m_cand_pct=ens_m_cand/value * 100,
         rank=min_rank(value)) |>
  separate_wider_delim(ens_id, "_", names=c("x", "D", "nSim", "sample"), cols_remove=F) |>
  select(-x)

RMSE_date <- CV_df |>
  group_by(date, ens_id) |>
  mutate(N=n()) |>
  filter(N >= 30) |>
  summarise(across(contains("pred"), ~rmse_vec(.x, truth=licePerFish_rtrt))) |>
  group_by(ens_id) |>
  summarise(across(contains("pred"), mean)) |>
  pivot_longer(contains("pred")) |>
  arrange(ens_id, name) |>
  group_by(ens_id) |>
  mutate(ens_m_cand=last(value)-value,
         ens_m_cand_pct=ens_m_cand/value * 100,
         rank=min_rank(value)) |>
  separate_wider_delim(ens_id, "_", names=c("x", "D", "nSim", "sample"), cols_remove=F) |>
  select(-x)

metric_df <- bind_rows(
  RMSE_site |> mutate(metric="RMSE", type="site"),
  RMSE_date |> mutate(metric="RMSE", type="date")
) |>
  mutate(nSim=factor(nSim, levels=levels(mod_ids$nSims)),
         ens_id=factor(ens_id, levels=levels(mod_ids$ens_id_ord)))



# summarize and visualize -------------------------------------------------

metric_df |> 
  group_by(name) |>
  summarise(mnRank=mean(rank, na.rm=T),
            prop1=mean(rank==1, na.rm=T)) |> 
  arrange(mnRank)

metric_df |> 
  filter(nSim != "n20") |>
  group_by(nSim, name) |>
  summarise(mnRank=mean(rank, na.rm=T),
            prop1=mean(rank==1, na.rm=T)) |> 
  arrange(nSim, mnRank) |>
  print(n=50)

metric_df |>
  filter(nSim != "n20") |>
  ggplot(aes(ens_m_cand, D, fill=nSim)) + 
  geom_vline(xintercept=0) +
  geom_boxplot() + 
  facet_grid(type~metric, scales="free") 

metric_df |>
  ggplot(aes(ens_m_cand_pct, D, fill=nSim)) + 
  geom_vline(xintercept=0) +
  geom_boxplot() + 
  facet_grid(type~metric, scales="free") 

metric_df |>
  ggplot(aes(value, ens_id, colour=name=="ens_pred")) + 
  geom_point(shape=1) +
  facet_grid(nSim~metric*type, scales="free")

metric_df |>
  group_by(ens_id, type, metric) |>
  arrange(rank) |>
  mutate(ensBest=first(name)=="ens_pred") |>
  ungroup() |>
  filter(name != "ens_pred") |>
  ggplot(aes(ens_m_cand, ens_id, colour=ensBest)) + 
  geom_vline(xintercept=0) +
  geom_point(shape=1, alpha=0.75) +
  geom_line() +
  scale_colour_manual(values=c("grey50", "green4")) +
  facet_grid(nSim~type, scales="free")


metric_df |>
  group_by(ens_id, type, metric) |>
  arrange(rank) |>
  mutate(ensBest=first(name)=="ens_pred") |>
  arrange(name) |>
  mutate(ensPct=100-last(percent_rank(rank))*100) |>
  ungroup() |>
  ggplot(aes(value, ens_id, colour=ensPct)) + 
  geom_point(aes(shape=name=="ens_pred", size=name=="ens_pred"), alpha=0.75) +
  scale_colour_distiller("Ensemble\npercentile", 
                         limits=c(0,100), direction=1) +
  scale_shape_manual(values=c(1, 4)) +
  scale_size_manual(values=c(1, 2)) +
  facet_grid(nSim~type, scales="free")

temp_df <- metric_df |>
  filter(D=="D3") |>
  group_by(ens_id, type, metric) |>
  arrange(rank) |>
  mutate(ensBest=if_else(first(name)=="ens_pred", "Ens['Blend']", "Candidate")) |>
  arrange(name) |>
  mutate(ensPct=100-last(percent_rank(rank))*100) |>
  ungroup() |>
  mutate(modType=if_else(name=="ens_pred", "Ens['Blend']", "Candidate"), 
         sample=factor(sample, levels=10:1),
         type=factor(type, levels=c("site", "date"),
                     labels=paste("Mean among", c("farms", "weeks"))),
         nSim=factor(nSim, levels=paste0("n", c(2, 5, 10, 20)),
                     labels=paste0("n: ", c(2, 5, 10, 20))))
lab_expressions <- c(expression(Candidate),
                     expression(Ens['Blend']))
p <- temp_df |>
  ggplot(aes(value, sample, colour=ensBest)) + 
  geom_point(aes(shape=modType, size=modType), alpha=0.75) +
  geom_line(data=temp_df |> filter(name != "ens_pred")) +
  scale_colour_manual("Best performance", values=c("#7fcdbb", "#0c2c84"),
                      labels=lab_expressions) +
  scale_shape_manual("Model type", values=c(1, 4),
                     labels=lab_expressions) +
  scale_size_manual("Model type", values=c(1, 2),
                    labels=lab_expressions) +
  facet_grid(nSim~type, scales="free_y", space="free_y") +
  labs(x="RMSE", y="Sample from pool of 20 candidates")
ggsave("figs/pub/ensBlend_eval_RMSE.png", p, width=6, height=10.5)
