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
ensFull_df <- read_csv("out/valid_df_2021-2024.csv") |>
  mutate(across(starts_with("sim_"), ~.x - mean(.x), .names="c_{.col}"))
folds <- unique(ensFull_df$CV_k)

sim_i <- read_csv("out/sim_2021-2024/sim_i.csv") |>
  mutate(sim=paste0("sim_", i),
         lab_short=if_else(fixDepth, "2D", "3D")) |>
  group_by(lab_short) |>
  mutate(lab=paste0(lab_short, ".", row_number())) |>
  ungroup() 


# individual simulation models --------------------------------------------

CV_sims <- vector("list", length(folds))
for(k in seq_along(folds)) {
  test_rows <- which(ensFull_df$CV_k == folds[k])
  CV_k <- vector("list", nrow(sim_i))
  for(i in 1:nrow(sim_i)) {
    sim <- sim_i$sim[i]
    dat_rstan <- ensFull_df |>
      select(rowNum, sepaSite, sepaSiteNum, date, 
             licePerFish_rtrt, any_of(sim_i$sim), any_of(paste0("c_", sim_i$sim))) |>
        select(-any_of(sim_i$sim[-i]), -any_of(paste0("c_", sim_i$sim[-i]))) |>
      make_data_rstan_GQ(test_rows)
    fit_ik <- glue("out/candidates/{sim}_CV-{folds[k]}_stanfit.rds")
    if(file.exists(fit_ik)) {
      out_sim <- readRDS(fit_ik)
    } else {
      out_sim <- stan(file="code/stan/candidate_model_GQ.stan",
                      model_name=glue("{sim}-{folds[k]}"), data=dat_rstan,
                      chains=6, cores=2,iter=3000, warmup=2500,
                      pars=c("b_b0", "b_IP", "Intercept_hu", "b_hu", "sigma", "GQ_Ypred"))
      # saveRDS(out_sim, fit_ik)
    }
    CV_k[[i]] <- colMeans(rstan::extract(out_sim, pars="GQ_Ypred")[[1]]) |>
      as_tibble() |>
      set_names(paste0("IP_", sim)) |>
      mutate(rowNum=test_rows)
  }
  CV_sims[[k]] <- inner_join(ensFull_df |> select(rowNum), 
                             reduce(CV_k, full_join, by=join_by(rowNum)),
                             by=join_by(rowNum))
}
reduce(CV_sims, bind_rows) |>
  write_csv("out/candidates/CV_candidate_predictions.csv")



# unweighted mean models --------------------------------------------------

sim_avgs <- c("sim_avgAll", "sim_avg3D", "sim_avg2D")
CV_avgs <- vector("list", length(folds))
for(k in seq_along(folds)) {
  test_rows <- which(ensFull_df$CV_k == folds[k])
  CV_k <- vector("list", length(sim_avgs))
  for(i in 1:length(sim_avgs)) {
    sim <- sim_avgs[i]
    dat_rstan <- ensFull_df |>
      select(rowNum, sepaSite, sepaSiteNum, date, 
             licePerFish_rtrt, starts_with("sim_avg"), starts_with("c_sim_avg")) |>
      select(-any_of(sim_avgs[-i]), -any_of(paste0("c_", sim_avgs[-i]))) |>
      make_data_rstan_GQ(test_rows)
    fit_ik <- glue("out/ensembles/{sim}_CV-{folds[k]}_stanfit.rds")
    if(file.exists(fit_ik)) {
      out_sim <- readRDS(fit_ik)
    } else {
    out_sim <- stan(file="code/stan/candidate_model_GQ.stan",
                    model_name=glue("{sim}-{folds[k]}"), data=dat_rstan,
                    chains=6, cores=6,iter=3000, warmup=2500,
                    pars=c("b_b0", "b_IP", "Intercept_hu", "b_hu", "sigma", "GQ_Ypred"))
    saveRDS(out_sim, fit_ik)
    }
    CV_k[[i]] <- colMeans(rstan::extract(out_sim, pars="GQ_Ypred")[[1]]) |>
      as_tibble() |>
      set_names(paste0("IP_", sim)) |>
      mutate(rowNum=test_rows)
  }
  CV_avgs[[k]] <- inner_join(ensFull_df |> select(rowNum), 
                             reduce(CV_k, full_join, by=join_by(rowNum)),
                             by=join_by(rowNum))
}
reduce(CV_avgs, bind_rows) |>
  write_csv("out/ensembles/CV_avg_predictions.csv")


