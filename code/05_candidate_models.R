

library(tidyverse)
library(glue)
library(rstan)
rstan_options(auto_write=T)

source("code/00_fn.R")

# Full dataset
ensFull_df <- read_csv("out/valid_df.csv") |>
  mutate(across(starts_with("sim_"), ~.x - mean(.x), .names="c_{.col}"))
years <- unique(ensFull_df$year)

sim_i <- read_csv("out/sim_2019-2023/sim_i.csv") |>
  mutate(sim=paste0("sim_", i),
         lab_short=if_else(fixDepth, "2D", "3D")) |>
  group_by(lab_short) |>
  mutate(lab=paste0(lab_short, ".", row_number())) 


# individual simulation models --------------------------------------------

CV_sims <- vector("list", length(years))
for(yr in seq_along(years)) {
  train_df <- ensFull_df |>
    filter(year != years[yr]) |>
    select(rowNum, sepaSite, sepaSiteNum, date, 
           licePerFish_rtrt, any_of(sim_i$sim), any_of(paste0("c_", sim_i$sim)))
  test_df <- ensFull_df |>
    filter(year == years[yr]) |>
    select(rowNum, sepaSite, sepaSiteNum, date,
           licePerFish_rtrt, any_of(sim_i$sim), any_of(paste0("c_", sim_i$sim)))
  CV_yr <- vector("list", nrow(sim_i))
  for(i in 1:nrow(sim_i)) {
    sim <- sim_i$sim[i]
    dat_rstan <- train_df |>
      select(-any_of(sim_i$sim[-i]), -any_of(paste0("c_", sim_i$sim[-i]))) |>
      make_data_rstan()
    if(file.exists(glue("out/candidates/{sim}_CV-{years[yr]}_stanfit.rds"))) {
      out_sim <- readRDS(glue("out/candidates/{sim}_CV-{years[yr]}_stanfit.rds"))
    } else {
      out_sim <- stan(file="code/stan/candidate_model.stan",
                      model_name=glue("{sim}-{years[yr]}"), data=dat_rstan,
                      chains=6, cores=6,iter=3000, warmup=2500,
                      pars=c("b_b0", "b_IP", "Intercept_hu", "b_hu", "sigma"))
      saveRDS(out_sim, glue("out/candidates/{sim}_CV-{years[yr]}_stanfit.rds"))
    }
    CV_yr[[i]] <- make_predictions_candidate(out_sim, test_df, sim) |>
      colMeans() |>
      as_tibble() |>
      set_names(paste0("IP_", sim))
  }
  CV_sims[[yr]] <- bind_cols(test_df |> select(rowNum), reduce(CV_yr, bind_cols))
}
reduce(CV_sims, bind_rows) |>
  write_csv("out/candidates/CV_candidate_predictions.csv")




# unweighted mean models --------------------------------------------------


sim_avgs <- c("sim_avg3D", "sim_avg2D")
CV_avgs <- vector("list", length(years))
for(yr in seq_along(years)) {
  train_df <- ensFull_df |>
    filter(year != years[yr]) |>
    select(rowNum, sepaSite, sepaSiteNum, date, 
           licePerFish_rtrt, starts_with("sim_avg"), starts_with("c_sim_avg"))
  test_df <- ensFull_df |>
    filter(year == years[yr]) |>
    select(rowNum, sepaSite, sepaSiteNum, date,
           licePerFish_rtrt, starts_with("sim_avg"), starts_with("c_sim_avg"))
  CV_yr <- vector("list", length(sim_avgs))
  for(i in 1:length(sim_avgs)) {
    sim <- sim_avgs[i]
    dat_rstan <- train_df |>
      select(-any_of(sim_avgs[-i]), -any_of(paste0("c_", sim_avgs[-i]))) |>
      make_data_rstan()
    if(file.exists(glue("out/ensembles/{sim}_CV-{years[yr]}_stanfit.rds"))) {
      out_sim <- readRDS(glue("out/ensembles/{sim}_CV-{years[yr]}_stanfit.rds"))
    } else {
    out_sim <- stan(file="code/stan/candidate_model.stan",
                    model_name=glue("{sim}-{years[yr]}"), data=dat_rstan,
                    chains=6, cores=6,iter=3000, warmup=2500,
                    pars=c("b_b0", "b_IP", "Intercept_hu", "b_hu", "sigma"))
    saveRDS(out_sim, glue("out/ensembles/{sim}_CV-{years[yr]}_stanfit.rds"))
    }
    CV_yr[[i]] <- make_predictions_candidate(out_sim, test_df, sim) |>
      colMeans() |>
      as_tibble() |>
      set_names(paste0("IP_", sim))
  }
  CV_avgs[[yr]] <- bind_cols(test_df |> select(rowNum), reduce(CV_yr, bind_cols))
}
reduce(CV_avgs, bind_rows) |>
  write_csv("out/ensembles/CV_avg_predictions.csv")


