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

sim_i <- read_csv("out/sim_2021-2024/sim_i.csv") |>
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
sims_3D <- filter(sim_i, lab_short=="3D")$sim
sims_2D <- filter(sim_i, lab_short=="2D")$sim

site_i <- read_csv("data/farm_sites.csv")
ensFull_LatLon <- read_csv("out/valid_df_2021-2024.csv") |>
  select(rowNum, date, CV_k, sepaSite, sepaSiteNum, licePerFish_rtrt, starts_with("sim")) |>
  select(-contains("avg")) |>
  mutate(across(starts_with("sim_"), ~.x - mean(.x), .names="c_{.col}")) |>
  left_join(site_i) |>
  select(-sepaSite) |>
  arrange(rowNum)


spline_rec_3D_ls <- map(
  1:10,
  ~recipe(licePerFish_rtrt ~ ., 
          data=ensFull_LatLon |>
            select(-any_of(sims_2D)) |>
            select(-any_of(paste0("c_", sims_2D)))) |>
    step_interact(~easting:northing) |>
    step_bs(easting, northing, easting_x_northing, deg_free=.x) |>
    update_role(c(rowNum, CV_k, date, sepaSiteNum, starts_with("c_")), new_role="id") |>
    update_role_requirements("id", bake=F) |>
    prep()
  )
saveRDS(spline_rec_3D_ls, "out/ensembles/recipe_sLonLat_3D.rds")

spline_rec_all_ls <- map(
  1:10,
  ~recipe(licePerFish_rtrt ~ ., data=ensFull_LatLon) |>
    step_interact(~easting:northing) |>
    step_bs(easting, northing, easting_x_northing, deg_free=.x) |>
    update_role(c(rowNum, CV_k, date, sepaSiteNum, starts_with("c_")), new_role="id") |>
    update_role_requirements("id", bake=F) |>
    prep()
)
saveRDS(spline_rec_all_ls, "out/ensembles/recipe_sLonLat_all.rds")



# cross-validation --------------------------------------------------------

mods <- expand_grid(mod=c("ranef", paste0("sLonLatD", 3:10)),
                    d=c("3D", "all"))
mods <- mods |> filter(d=="3D")

CV_ensMix <- vector("list", length(folds))

for(k in seq_along(folds)) {
  CV_k_ls <- vector("list", nrow(mods))
  
  for(i in 1:nrow(mods)) {
    mod <- mods$mod[i]
    d <- mods$d[i]
    
    if(grepl("sLonLat", mod)) {
      spline_deg_free <- as.numeric(str_split_fixed(mod, "D", 2)[2])
      
      if(d=="3D") {
        train_df <- spline_rec_3D_ls[[spline_deg_free]] |>
          bake(ensFull_LatLon |> 
                 filter(CV_k != folds[k]) |>
                 select(rowNum, easting, northing, sepaSiteNum, CV_k, date, licePerFish_rtrt,
                        any_of(sims_3D), 
                        any_of(paste0("c_", sims_3D))))
        test_df <- spline_rec_3D_ls[[spline_deg_free]] |>
          bake(ensFull_LatLon |> 
                 filter(CV_k == folds[k]) |>
                 select(rowNum, easting, northing, sepaSiteNum, CV_k, date, licePerFish_rtrt,
                        any_of(sims_3D), 
                        any_of(paste0("c_", sims_3D))))
      } else {
        train_df <- spline_rec_all_ls[[spline_deg_free]] |>
          bake(ensFull_LatLon |> 
                 filter(CV_k != folds[k]) |>
                 select(rowNum, easting, northing, sepaSiteNum, CV_k, date, licePerFish_rtrt,
                        any_of(c(sims_3D, sims_2D)), 
                        any_of(paste0("c_", c(sims_3D, sims_2D)))))
        test_df <- spline_rec_all_ls[[spline_deg_free]] |>
          bake(ensFull_LatLon |> 
                 filter(CV_k == folds[k]) |>
                 select(rowNum, easting, northing, sepaSiteNum, CV_k, date, licePerFish_rtrt,
                        any_of(c(sims_3D, sims_2D)), 
                        any_of(paste0("c_", c(sims_3D, sims_2D)))))
      }
    
      dat_rstan <- make_data_rstan_sLonLat(train_df)
      pars <- c("b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu",
                paste0("b_s_", c("easting", "northing", "easting_x_northing")))
      
    } else {
      train_df <- ensFull_df |>
        filter(CV_k != folds[k]) |>
        select(rowNum, sepaSite, sepaSiteNum, CV_k, date, licePerFish_rtrt,
               any_of(filter(sim_i, lab_short %in% c("3D", ifelse(d=="all", "2D", "")))$sim), 
               any_of(paste0("c_", filter(sim_i, lab_short %in% c("3D", ifelse(d=="all", "2D", "")))$sim)))
      test_df <- ensFull_df |>
        filter(CV_k == folds[k]) |>
        select(rowNum, sepaSite, sepaSiteNum, CV_k, date, licePerFish_rtrt, 
               any_of(filter(sim_i, lab_short %in% c("3D", ifelse(d=="all", "2D", "")))$sim), 
               any_of(paste0("c_", filter(sim_i, lab_short %in% c("3D", ifelse(d=="all", "2D", "")))$sim)))
      dat_rstan <- make_data_rstan(train_df, R2D2_sd=F)
      pars <- c("b_p", "b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu",
                "b_p_uc", "r_grp", "sd_grp")
    }
  
    if(grepl("pRE_huRE", mod)) {
      pars <- c(pars, "r_grp_hu", "sd_grp_hu") 
    }
  
    if(file.exists(glue("out/ensembles/ensMix_{d}_{mod}_CV-{folds[k]}_stanfit.rds"))) {
      cat("File exists:", d, mod, folds[k], "\n")
      out_ensMix <- readRDS(glue("out/ensembles/ensMix_{d}_{mod}_CV-{folds[k]}_stanfit.rds"))
    } else {
      stanMod <- str_split_fixed(mod, "D", 2)[1]
      out_ensMix <- stan(file=glue("code/stan/ensemble_mixture_model_{stanMod}.stan"),
                         model_name=mod, data=dat_rstan,
                         # chains=30, cores=30, iter=2250, warmup=2000,
                         chains=3, cores=3, iter=3000, warmup=2000,
                         control=list(adapt_delta=0.95, max_treedepth=20),
                         pars=pars)
      saveRDS(out_ensMix, glue("out/ensembles/ensMix_{d}_{mod}_CV-{folds[k]}_stanfit.rds"))
      saveRDS(dat_rstan, glue("out/ensembles/ensMix_{d}_{mod}_CV-{folds[k]}_standata.rds")) 
    }
  
    if(grepl("sLonLat", mod)) {
      CV_k_ls[[i]] <- make_predictions_ensMix_sLonLat(out_ensMix, test_df, mode="epred", iter=10) |>
        colMeans() |>
        as_tibble() |>
        set_names(paste0("IP_", mod, "_", d))
    } else {
      CV_k_ls[[i]] <- make_predictions_ensMix(out_ensMix, test_df, re=T, re_hu=grepl("huRE", mod)) |>
        colMeans() |>
        as_tibble() |>
        set_names(paste0("IP_", mod, "_", d)) 
    }
  }
  CV_ensMix[[k]] <- bind_cols(test_df |> select(rowNum), reduce(CV_k_ls, bind_cols))
}
reduce(CV_ensMix, bind_rows) |>
  write_csv("out/ensembles/CV_ensMix_predictions.csv")




# full dataset ------------------------------------------------------------

mods <- expand_grid(mod=c("ranef", paste0("sLonLatD", 3:10)),
                    d=c("3D", "all"))
mods <- mods |> filter(d=="3D")

for(i in 1:nrow(mods)) {
  mod <- mods$mod[i]
  d <- mods$d[i]
  
  if(grepl("sLonLat", mod)) {
    spline_deg_free <- as.numeric(str_split_fixed(mod, "D", 2)[2])
  
    if(d=="3D") {
      train_df <- spline_rec_3D_ls[[spline_deg_free]] |>
        bake(ensFull_LatLon |> 
               select(rowNum, easting, northing, sepaSiteNum, date, licePerFish_rtrt,
                      any_of(sims_3D), 
                      any_of(paste0("c_", sims_3D))))
    } else {
      train_df <- spline_rec_all_ls[[spline_deg_free]] |>
        bake(ensFull_LatLon |> 
               select(rowNum, easting, northing, sepaSiteNum, date, licePerFish_rtrt,
                      any_of(c(sims_3D, sims_2D)), 
                      any_of(paste0("c_", c(sims_3D, sims_2D)))))
    }
    dat_rstan <- make_data_rstan_sLonLat(train_df)
    pars <- c("b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu",
              paste0("b_s_", c("easting", "northing", "easting_x_northing")))
  } else {
    train_df <- ensFull_df |>
      select(rowNum, sepaSite, sepaSiteNum, date, licePerFish_rtrt,
             any_of(filter(sim_i, lab_short %in% c("3D", ifelse(d=="all", "2D", "")))$sim),
             any_of(paste0("c_", filter(sim_i, lab_short %in% c("3D", ifelse(d=="all", "2D", "")))$sim)))
    dat_rstan <- make_data_rstan(train_df, R2D2_sd=F)
    pars <- c("b_p", "b_b0", "b_IP", "b_hu", "sigma", "Intercept_hu",
              "b_p_uc", "r_grp", "sd_grp")
    if(grepl("pRE_huRE", mod)) {
      pars <- c(pars, "r_grp_hu", "sd_grp_hu") 
    } 
  }
  if(!file.exists(glue("out/ensembles/ensMix_{d}_{mod}_FULL_stanfit.rds"))) {
    stanMod <- str_split_fixed(mod, "D", 2)[1]
    out_ensMix <- stan(file=glue("code/stan/ensemble_mixture_model_{stanMod}.stan"),
                       model_name=mod, data=dat_rstan,
                       # chains=30, cores=15, iter=2250, warmup=2000,
                       chains=3, cores=3, iter=3000, warmup=2000,
                       control=list(adapt_delta=0.95, max_treedepth=20),
                       pars=pars)
    saveRDS(out_ensMix, glue("out/ensembles/ensMix_{d}_{mod}_FULL_stanfit.rds"))
    saveRDS(dat_rstan, glue("out/ensembles/ensMix_{d}_{mod}_FULL_standata.rds"))
  } else {
    cat("File exists:", d, mod, "\n")
    Sys.sleep(3)
  }
}
