# Project: Sealice IP Ensemble
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Individual variants and mean IP

library(tidyverse)
library(glue)
library(rstan)
rstan_options(auto_write=T)

source("code/00_fn.R")

# Full dataset
ensFull_df <- read_csv("out/valid_sens_df_2021-2024.csv") |>
  group_by(survAdjust, devAdjust) |>
  mutate(across(starts_with("sim_"), ~.x - mean(.x), .names="c_{.col}")) |>
  ungroup()
folds <- unique(ensFull_df$CV_k)

sim_i <- read_csv("out/sim_2021-2024/sim_i.csv") |>
  mutate(sim=paste0("sim_", i),
         lab_short=if_else(fixDepth, "2D", "3D")) |>
  group_by(lab_short) |>
  mutate(lab=paste0(lab_short, ".", row_number())) |>
  ungroup() 

sens_df <- ensFull_df |>
  slice_head(n=1, by=ends_with("Adjust")) |>
  select(ends_with("Adjust"))



# individual simulation models --------------------------------------------

# Cross-validation
CV_k <- vector("list", length(folds))
for(k in seq_along(folds)) {
  CV_j <- vector("list", nrow(sens_df))
  for(j in 1:nrow(sens_df)) {
    ensFull_df_j <- ensFull_df |>
      inner_join(sens_df[j,], by=join_by(survAdjust, devAdjust)) |>
      select(-ends_with("Adjust"))
    test_rows <- which(ensFull_df_j$CV_k == folds[k])
    CV_i <- vector("list", nrow(sim_i))
    for(i in 1:nrow(sim_i)) {
      sim <- sim_i$sim[i]
      dat_rstan <- ensFull_df_j |>
        select(rowNum, sepaSite, sepaSiteNum, date,
               licePerFish_rtrt, any_of(sim_i$sim), any_of(paste0("c_", sim_i$sim))) |>
        select(-any_of(sim_i$sim[-i]), -any_of(paste0("c_", sim_i$sim[-i]))) |>
        make_data_rstan_GQ(test_rows)
      fit_ijk <- glue("out/candidates/AEIP_sens/{sim}_s{sens_df[j,1]}_d{sens_df[j,2]}_CV-{folds[k]}_stanfit.rds")
      if(file.exists(fit_ijk)) {
        out_sim <- readRDS(fit_ijk)
      } else {
        out_sim <- stan(file="code/stan/constituent_model_GQ.stan",
                        model_name=glue("{sim}-{folds[k]}"), data=dat_rstan,
                        chains=6, cores=6,iter=3000, warmup=2500,
                        pars=c("b_b0", "b_IP", "Intercept_hu", "b_hu", "sigma", "GQ_Ypred"))
        saveRDS(out_sim, fit_ijk)
      }
      CV_i[[i]] <- colMeans(rstan::extract(out_sim, pars="GQ_Ypred")[[1]]) |>
        as_tibble() |>
        set_names(paste0("IP_", sim)) |>
        mutate(rowNum=test_rows)
    }  
    CV_j[[j]] <- reduce(CV_i, full_join, by=join_by(rowNum)) |>
      mutate(survAdjust=sens_df$survAdjust[j],
             devAdjust=sens_df$devAdjust[j])
  }
  CV_k[[k]] <- inner_join(ensFull_df |> select(rowNum, survAdjust, devAdjust),
                          bind_rows(CV_j),
                          by=join_by(rowNum, survAdjust, devAdjust))
}
reduce(CV_k, bind_rows) |>
  write_csv("out/candidates/AEIP_sens/CV_candidate_predictions.csv")




# unweighted mean models --------------------------------------------------

# Cross-validation
sim_avgs <- c("sim_avgAll", "sim_avg2D", "sim_avg3D")
CV_k <- vector("list", length(folds))
for(k in seq_along(folds)) {
  CV_j <- vector("list", nrow(sens_df))
  for(j in 1:nrow(sens_df)) {
    ensFull_df_j <- ensFull_df |>
      inner_join(sens_df[j,], by=join_by(survAdjust, devAdjust)) |>
      select(-ends_with("Adjust"))
    test_rows <- which(ensFull_df_j$CV_k == folds[k])
    CV_i <- vector("list", length(sim_avgs))
    for(i in 1:length(sim_avgs)) {
      sim <- sim_avgs[i]
      dat_rstan <- ensFull_df_j |>
        select(rowNum, sepaSite, sepaSiteNum, date, 
               licePerFish_rtrt, starts_with("sim_avg"), starts_with("c_sim_avg")) |>
        select(-any_of(sim_avgs[-i]), -any_of(paste0("c_", sim_avgs[-i]))) |>
        make_data_rstan_GQ(test_rows)
      fit_ijk <- glue("out/ensembles/AEIP_sens/{sim}_s{sens_df[j,1]}_d{sens_df[j,2]}_CV-{folds[k]}_stanfit.rds")
      if(file.exists(fit_ijk)) {
        out_sim <- readRDS(fit_ijk)
      } else {
        out_sim <- stan(file="code/stan/constituent_model_GQ.stan",
                        model_name=glue("{sim}-{folds[k]}"), data=dat_rstan,
                        chains=6, cores=6,iter=3000, warmup=2500,
                        pars=c("b_b0", "b_IP", "Intercept_hu", "b_hu", "sigma", "GQ_Ypred"))
        saveRDS(out_sim, fit_ijk)
      }
      CV_i[[i]] <- colMeans(rstan::extract(out_sim, pars="GQ_Ypred")[[1]]) |>
        as_tibble() |>
        set_names(paste0("IP_", sim)) |>
        mutate(rowNum=test_rows)
    }
    CV_j[[j]] <- reduce(CV_i, full_join, by=join_by(rowNum)) |>
      mutate(survAdjust=sens_df$survAdjust[j],
             devAdjust=sens_df$devAdjust[j])
  }
  CV_k[[k]] <- inner_join(ensFull_df |> select(rowNum, survAdjust, devAdjust), 
                          bind_rows(CV_j),
                          by=join_by(rowNum, survAdjust, devAdjust))
}
reduce(CV_k, bind_rows) |>
  write_csv("out/ensembles/AEIP_sens/CV_avg_predictions.csv")


