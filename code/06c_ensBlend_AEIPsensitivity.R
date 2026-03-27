# Project: Sealice IP Ensemble
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Blending ensemble

library(tidyverse)
library(glue)
library(rstan)
library(recipes)
rstan_options(auto_write=T)

source("code/00_fn.R")


# ensembling --------------------------------------------------------------

# Full dataset: includes +/- 10% for survival and development times in AEIP calculation
ensFull_df <- read_csv("out/valid_sens_df_2021-2024.csv") |>
  group_by(survAdjust, devAdjust) |>
  mutate(across(starts_with("sim_"), ~.x - mean(.x), .names="c_{.col}")) |>
  ungroup()
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
  select(ends_with("Adjust"), rowNum, date, CV_k, sepaSite, sepaSiteNum, licePerFish_rtrt, starts_with("sim")) |>
  select(-contains("avg")) |>
  mutate(across(starts_with("sim_"), ~.x - mean(.x), .names="c_{.col}")) |>
  left_join(site_i) |>
  select(-sepaSite) |>
  arrange(rowNum)



# full dataset ------------------------------------------------------------

sens_df <- ensFull_LatLon |>
  slice_head(n=1, by=ends_with("Adjust")) |>
  select(ends_with("Adjust"))
  
sim_ls <- list("n20"=filter(sim_i, n20)$sim)

for(i in 1:nrow(sens_df)) {
  ensFull_LatLon_i <- ensFull_LatLon |>
    inner_join(sens_df[i,], by=join_by(survAdjust, devAdjust)) |>
    select(-ends_with("Adjust"))
  recipe_i <- make_spline_recipe(ensFull_LatLon_i, 
                                 4, 
                                 sim_ls[[1]])
  full_df <- bake(recipe_i, ensFull_LatLon)
  dat_rstan <- make_data_rstan_sLonLat(full_df)
  pars <- c("b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu", "b_p",
            paste0("b_s_", c("easting", "northing", "easting_x_northing")))
  if(!file.exists(glue("out/ensembles/AEIP_sens/ensBlend_s{sens_df[i,1]}_d{sens_df[i,2]}_stanfit.rds"))) {
    out_ensBlend <- stan(file="code/stan/ensBlend_model.stan",
                         model_name="AEIP_sens", data=dat_rstan,
                         chains=3, cores=3, iter=3000, warmup=2000, refresh=100,
                         # chains=3, cores=3, iter=3000, warmup=2000,
                         control=list(adapt_delta=0.95, max_treedepth=20),
                         pars=pars)
    saveRDS(out_ensBlend, glue("out/ensembles/AEIP_sens/ensBlend_s{sens_df[i,1]}_d{sens_df[i,2]}_stanfit.rds"))
    saveRDS(dat_rstan, glue("out/ensembles/AEIP_sens/ensBlend_s{sens_df[i,1]}_d{sens_df[i,2]}_standata.rds"))
  } else {
    Sys.sleep(3)
  }
}
