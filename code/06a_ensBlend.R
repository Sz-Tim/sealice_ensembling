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

# Full dataset
ensFull_df <- read_csv("out/valid_df_2021-2024.csv") |>
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
  mutate(n5=row_number() %in% sample(1:20, 5),
         n10=row_number() %in% sample(1:20, 10),
         n20=row_number() %in% 1:20,
         best5=row_number() %in% c(4, 10, 9, 8, 20))

site_i <- read_csv("data/farm_sites.csv")
ensFull_LatLon <- read_csv("out/valid_df_2021-2024.csv") |>
  select(rowNum, date, CV_k, sepaSite, sepaSiteNum, licePerFish_rtrt, starts_with("sim")) |>
  select(-contains("avg")) |>
  mutate(across(starts_with("sim_"), ~.x - mean(.x), .names="c_{.col}")) |>
  left_join(site_i) |>
  select(-sepaSite) |>
  arrange(rowNum)



# cross-validation --------------------------------------------------------

mods <- expand_grid(mod=c("ranef", paste0("sLonLatD", 3:10)),
                    nSims=c(paste0("n", c(5, 10, 20)), "best5"))

sim_ls <- list("n5"=filter(sim_i, n5)$sim,
               "n10"=filter(sim_i, n10)$sim,
               "n20"=filter(sim_i, n20)$sim,
               "best5"=filter(sim_i, best5)$sim)

CV_ensBlend <- vector("list", length(folds))

for(k in (seq_along(folds))) {
  CV_k_ls <- vector("list", nrow(mods))
  
  for(i in 1:nrow(mods)) {
    mod <- mods$mod[i]
    nSims <- mods$nSims[i]
    
    if(grepl("sLonLat", mod)) {
      recipe_i <- make_spline_recipe(ensFull_LatLon, 
                                     as.numeric(str_split_fixed(mod, "D", 2)[2]), 
                                     sim_ls[[nSims]])
      train_df <- bake(recipe_i, ensFull_LatLon |> filter(CV_k != folds[k]))
      test_df <- bake(recipe_i, ensFull_LatLon |> filter(CV_k == folds[k]))
      dat_rstan <- make_data_rstan_sLonLat(train_df)
      pars <- c("b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu",
                paste0("b_s_", c("easting", "northing", "easting_x_northing")))
    } else {
      train_df <- ensFull_df |> filter(CV_k != folds[k]) |>
        select(rowNum, sepaSite, sepaSiteNum, CV_k, date, licePerFish_rtrt,
               all_of(sim_ls[[nSims]]), all_of(paste0("c_", sim_ls[[nSims]])))
      test_df <- ensFull_df |> filter(CV_k == folds[k]) |>
        select(rowNum, sepaSite, sepaSiteNum, CV_k, date, licePerFish_rtrt,
               all_of(sim_ls[[nSims]]), all_of(paste0("c_", sim_ls[[nSims]])))
      dat_rstan <- make_data_rstan(train_df, R2D2_sd=F)
      pars <- c("b_p", "b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu",
                "b_p_uc", "r_grp", "sd_grp")
    }
    
    if(grepl("pRE_huRE", mod)) {
      pars <- c(pars, "r_grp_hu", "sd_grp_hu") 
    }
    
    fname <- glue("out/ensembles/ensBlend_{nSims}_{mod}_CV-{folds[k]}")
    if(file.exists(glue("{fname}_stanfit.rds"))) {
      cat("File exists:", fname, "\n")
      out_ensBlend <- readRDS(glue("{fname}_stanfit.rds"))
    } else {
      stanMod <- str_split_fixed(mod, "D", 2)[1]
      out_ensBlend <- stan(file=glue("code/stan/ensemble_mixture_model_{stanMod}.stan"),
                           model_name=mod, data=dat_rstan,
                           # chains=30, cores=30, iter=2250, warmup=2000,
                           chains=3, cores=3, iter=3000, warmup=2000,
                           control=list(adapt_delta=0.95, max_treedepth=20),
                           pars=pars)
      saveRDS(out_ensBlend, glue("{fname}_stanfit.rds"))
      saveRDS(dat_rstan, glue("{fname}_standata.rds")) 
    }
    
    if(grepl("sLonLat", mod)) {
      CV_k_ls[[i]] <- out_ensBlend |>
        make_predictions_ensBlend_sLonLat(test_df, mode="point_epred") |>
        colMeans() |>
        as_tibble() |>
        set_names(paste0("IP_", mod, "_", nSims))
    } else {
      CV_k_ls[[i]] <- out_ensBlend |>
        make_predictions_ensBlend(test_df, re=T, re_hu=grepl("huRE", mod)) |>
        colMeans() |>
        as_tibble() |>
        set_names(paste0("IP_", mod, "_", nSims)) 
    }
  }
  CV_ensBlend[[k]] <- bind_cols(test_df |> select(rowNum), 
                                reduce(CV_k_ls, bind_cols))
}
reduce(CV_ensBlend, bind_rows) |>
  write_csv("out/ensembles/CV_ensBlend_predictions.csv")




# full dataset ------------------------------------------------------------

mods <- expand_grid(mod=c("ranef", paste0("sLonLatD", 3:10)),
                    nSims=c(paste0("n", c(5, 10, 20)), "best5"))

sim_ls <- list("n5"=filter(sim_i, n5)$sim,
               "n10"=filter(sim_i, n10)$sim,
               "n20"=filter(sim_i, n20)$sim,
               "best5"=filter(sim_i, best5)$sim)

for(i in rev(1:nrow(mods))) {
  mod <- mods$mod[i]
  nSims <- mods$nSims[i]
  
  if(grepl("sLonLat", mod)) {
    recipe_i <- make_spline_recipe(ensFull_LatLon, 
                                   as.numeric(str_split_fixed(mod, "D", 2)[2]), 
                                   sim_ls[[nSims]])
    train_df <- bake(recipe_i, ensFull_LatLon)
    dat_rstan <- make_data_rstan_sLonLat(train_df)
    pars <- c("b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu",
              paste0("b_s_", c("easting", "northing", "easting_x_northing")))
  } else {
    train_df <- ensFull_df |> 
      select(rowNum, sepaSite, sepaSiteNum, CV_k, date, licePerFish_rtrt,
             all_of(sim_ls[[nSims]]), all_of(paste0("c_", sim_ls[[nSims]])))
    dat_rstan <- make_data_rstan(train_df, R2D2_sd=F)
    pars <- c("b_p", "b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu",
              "b_p_uc", "r_grp", "sd_grp")
  }
  if(grepl("pRE_huRE", mod)) {
    pars <- c(pars, "r_grp_hu", "sd_grp_hu") 
  } 
  if(!file.exists(glue("out/ensembles/ensBlend_{nSims}_{mod}_FULL_stanfit.rds"))) {
    stanMod <- str_split_fixed(mod, "D", 2)[1]
    out_ensBlend <- stan(file=glue("code/stan/ensemble_mixture_model_{stanMod}.stan"),
                         model_name=mod, data=dat_rstan,
                         # chains=30, cores=15, iter=2250, warmup=2000,
                         chains=3, cores=3, iter=3000, warmup=2000,
                         control=list(adapt_delta=0.95, max_treedepth=20),
                         pars=pars)
    saveRDS(out_ensBlend, glue("out/ensembles/ensBlend_{nSims}_{mod}_FULL_stanfit.rds"))
    saveRDS(dat_rstan, glue("out/ensembles/ensBlend_{nSims}_{mod}_FULL_standata.rds"))
  } else {
    cat("File exists:", nSims, mod, "\n")
    Sys.sleep(3)
  }
}
