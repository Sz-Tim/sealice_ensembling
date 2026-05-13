# Project: Sealice IP Ensemble
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Blending ensemble

library(tidyverse)
library(glue)
library(rstan)
library(recipes)
library(doFuture)
rstan_options(auto_write=T)

source("code/00_fn.R")


# ensembling --------------------------------------------------------------

ensFull_df <- read_csv("out/valid_df_2021-2024_inclPostTreat.csv") |>
  mutate(across(starts_with("sim_"), ~.x - mean(.x), .names="c_{.col}"))
folds <- unique(ensFull_df$CV_k)

set.seed(66)
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
                                   "2D", "3D", "Null"))) |>
  mutate(n20=row_number() %in% 1:20)

site_i <- read_csv("data/farm_sites.csv")
ensFull_LatLon <- ensFull_df |>
  select(rowNum, date, CV_k, sepaSite, sepaSiteNum, licePerFish_rtrt, starts_with("sim")) |>
  select(-contains("avg")) |>
  mutate(across(starts_with("sim_"), ~.x - mean(.x), .names="c_{.col}")) |>
  left_join(site_i) |>
  select(-sepaSite) |>
  arrange(rowNum)

sim_ls <- list("n20"=filter(sim_i, n20)$sim)




# cross-validation --------------------------------------------------------

fit_only <- F
if(fit_only) plan(multicore, workers=30)

CV_k <- vector("list", length(folds))
# foreach(k=seq_along(folds)) %dofuture% {
for(k in seq_along(folds)) {
  test_rows <- which(ensFull_LatLon$CV_k == folds[k])
  recipe_i <- make_spline_recipe(ensFull_LatLon, 4, sim_ls[[1]])
  full_df <- bake(recipe_i, ensFull_LatLon)
  test_rows <- which(ensFull_LatLon$CV_k == folds[k])
  dat_rstan <- make_data_rstan_sLonLat_GQ(full_df, test_rows)
  pars <- c("GQ_Ypred")
  
  # fit EnsBlend
  fname <- glue("out/ensembles/TRT_sens/ensBlend_inclPostTreat_CV-{folds[k]}")
  if(!fit_only & file.exists(glue("{fname}_Ypred.rds"))) {
    cat("File exists:", fname, "\n")
    CV_k[[k]] <- readRDS(glue("{fname}_Ypred.rds"))
  } else {
    out_ensBlend <- stan(file="code/stan/ensBlend_model_GQ.stan",
                         model_name="postTreat", data=dat_rstan,
                         chains=3, cores=3, iter=3000, warmup=2000,
                         control=list(adapt_delta=0.95, max_treedepth=20),
                         pars=pars)
    saveRDS(dat_rstan, glue("{fname}_standata.rds")) 
    if(fit_only) {
      colMeans(rstan::extract(out_ensBlend, pars="GQ_Ypred")[[1]]) |>
        as_tibble() |>
        set_names("IP_D4_n20") |>
        mutate(rowNum=test_rows) |>
        saveRDS(glue("{fname}_Ypred.rds"))
    } else {
      CV_k[[k]] <- colMeans(rstan::extract(out_ensBlend, pars="GQ_Ypred")[[1]]) |>
        as_tibble() |>
        set_names("IP_D4_n20") |>
        mutate(rowNum=test_rows)
      CV_k[[k]] |>
        saveRDS(glue("{fname}_Ypred.rds"))
    }
  }
}

inner_join(ensFull_LatLon |> select(rowNum),
           bind_rows(CV_k)) |>
  rename_with(~str_remove(.x, "RE_GQ_"), starts_with("IP")) |>
  write_csv("out/ensembles/TRT_sens/CV_ensBlend_predictions.csv")



# full dataset ------------------------------------------------------------

recipe_i <- make_spline_recipe(ensFull_LatLon, 4, sim_ls[[1]])
full_df <- bake(recipe_i, ensFull_LatLon)
dat_rstan <- make_data_rstan_sLonLat(full_df)
pars <- c("b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu", "b_p",
          paste0("b_s_", c("easting", "northing", "easting_x_northing")))
if(!file.exists(glue("out/ensembles/TRT_sens/ensBlend_inclPostTreat_stanfit.rds"))) {
  out_ensBlend <- stan(file="code/stan/ensBlend_model.stan",
                       model_name="TRT_sens", data=dat_rstan,
                       chains=3, cores=3, iter=3000, warmup=2000, refresh=10,
                       control=list(adapt_delta=0.95, max_treedepth=20),
                       pars=pars)
  saveRDS(out_ensBlend, glue("out/ensembles/TRT_sens/ensBlend_inclPostTreat_stanfit.rds"))
  saveRDS(dat_rstan, glue("out/ensembles/TRT_sens/ensBlend_inclPostTreat_standata.rds"))
}
