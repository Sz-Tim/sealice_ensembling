# Project: Sealice IP Ensemble
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Forecasting ensembles

# This script tunes and validates forecasting ensembles, selecting the best 
# performing model. Optimization is performed separately for predicting the
# mean number of lice per fish and the probability that the mean lice per fish
# exceeds the Code of Good Practice treatment threshold. Optimization is also
# performed separately for predictions 1 and5 weeks in the future.


# setup -------------------------------------------------------------------
library(tidyverse); library(glue); #library(scico)
library(tidymodels); #library(DALEXtra); library(butcher)
library(baguette); # bag_mars
# library(discrim); # naive_Bayes
# library(bonsai); # lightgbm
library(finetune)
library(future)
library(sevcheck)
theme_set(theme_bw() + theme(panel.grid=element_blank()))
options(tidymodels.dark = TRUE)

gridSize <- 100
cores <- 50







# compile dataset ---------------------------------------------------------

ensFull_df <- read_csv("out/valid_df_2021-2024_FULL.csv") |>
  mutate(lice_g05=factor(licePerFish_rtrt^4 > 0.5)) |>
  left_join(read_csv("data/farm_sites.csv"))

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

fit_df <- ensFull_df |>
  select(rowNum, sepaSite, sepaSiteNum, productionCycleNumber, CV_k, date,
         easting, northing, licePerFish_rtrt, lice_g05, 
         any_of(filter(sim_i, lab_short %in% c("3D", "2D"))$sim)) |>
  arrange(sepaSite, date) |>
  drop_na()
rm(ensFull_df); gc()



# Ensemble formulas -------------------------------------------------------

base_recipe <- recipe(licePerFish_rtrt ~ ., data=fit_df) |>
  step_normalize() |>
  step_harmonic(date, frequency=1, cycle_size=365.24) |>
  step_interact(terms= ~date_sin_1:date_cos_1) |>
  step_interact(terms= ~easting:northing) |>
  step_bs(easting, deg_free=tune("bs_E_df")) |>
  step_bs(northing, deg_free=tune("bs_N_df")) |>
  step_bs(easting_x_northing, deg_free=tune("bs_EN_df")) |>
  update_role(c(rowNum, starts_with("sepaSite"), productionCycleNumber, CV_k,
                lice_g05, licePerFish_rtrt), new_role="not used") |>
  update_role_requirements("not used", bake=F)

recipes <- list(
  n20=base_recipe,
  PCA_n20=base_recipe |>
    step_pca(starts_with("sim") & has_role("predictor"), num_comp=tune())
)



# Candidate ensemble models -----------------------------------------------
# Ridge, Elastic Net, Lasso, Random Forest, XGBoost, Neural Network, CART, KNN
licePerFish_models <- list(
  enet=linear_reg(penalty=tune(), mixture=tune()) |>
    set_engine("glmnet") |> set_mode("regression"),
  mlp=bag_mlp(hidden_units=tune(), penalty=tune(), epochs=tune()) |>
    set_engine("nnet") |> set_mode("regression"),
  rf=rand_forest(trees=tune(), min_n=tune()) |>
    set_engine("randomForest") |> set_mode("regression")
)

liceBinary_models <- list(
  enet=logistic_reg(penalty=tune(), mixture=tune()) |>
    set_engine("glmnet") |> set_mode("classification"),
  mlp=bag_mlp(hidden_units=tune(), penalty=tune(), epochs=tune()) |>
    set_engine("nnet") |> set_mode("classification"),
  rf=rand_forest(trees=tune(), min_n=tune()) |>
    set_engine("randomForest") |> set_mode("classification")
)




# Tune ensembles: Dates ---------------------------------------------------
advance <- c(1) # DEPRECATED: number of weeks to forecast in advance
for(i in seq_along(advance)) {
  indices_ls <- map(unique(fit_df$CV_k), ~list(analysis=which(fit_df$CV_k != .x), 
                                 assessment=which(fit_df$CV_k == .x)))
  splits <- lapply(indices_ls, make_splits, data=fit_df)
  folds <- manual_rset(splits, paste("Split", seq_along(splits)))
  
  # DEPRECATED: 
  # expanding window cross-validation: simulate weekly forecasts within each year
  #   Fit model each week using all non-focus years + that year up to the week
  #   and predict the on-farm lice i weeks into the future
  
  # folds <- unique(fit_df$CV_k)
  # folds_ls <- vector("list", length(folds))
  # for(k in seq_along(folds)) {
  #   folds_ls[[k]] <- make_splits(fit_df, list(analysis=which(fit_df$CV_k!=k), assessment=which(fit_df$CV_k==k)))
  # }
  # for(k in seq_along(folds)) {
  #   fit_df_k <- fit_df |>
  #     mutate(date_mod=if_else(CV_k==folds[k], date+dyears(12), date)) |>
  #     arrange(date_mod)
  #   folds_ls[[k]] <- sliding_period(fit_df_k, index=date_mod, period="week", lookback=Inf, complete=F,
  #                                   skip=max(which(year(unique(fit_df_k$date_mod)) < 2025))-advance[i],
  #                                   assess_start=advance[i], assess_stop=advance[i])
  #   rm(fit_df_k)
  # }
  # folds_merged <- reduce(folds_ls, bind_rows) |>
  #   filter(map_lgl(splits, ~length(.x$out_id)>0))
  # folds <- manual_rset(folds_merged$splits, paste0("Slice", str_pad(1:nrow(folds_merged), 3, "left", "0")))
  fold_rowNums <- folds |>
    mutate(rowNum=map(splits, ~.x$data$rowNum[.x$out_id]),
           .row=map(splits, ~.x$out_id)) |>
    select(id, .row, rowNum) |>
    unnest(c(".row", "rowNum"))
  
  j <- 1
  if(j==1) {
    # licePerFish
    if(get_os()=="windows") {
      plan(multisession, workers=cores)
    } else {
      plan(multicore, workers=cores)
    }
    licePerFish_wfs <- workflow_set(
      preproc=map(recipes, ~.x |> update_role(licePerFish_rtrt, new_role="outcome")),
      models=licePerFish_models
    ) |>
      filter(!grepl("^PCA_.*enet", wflow_id)) |>
      workflow_map("tune_grid",
                   resamples=folds, 
                   grid=gridSize,
                   metrics=metric_set(rmse, rsq),
                   control=control_grid(save_pred=T,
                                        save_workflow=T,
                                        parallel_over="everything"),
                   verbose=T)
    plan(sequential)
    cat(format(Sys.time(), "%F %T"), "  Finished licePerFish tuning, advance:", advance[i], "\n")
    autoplot(licePerFish_wfs) + scale_colour_brewer(type="qual", palette="Paired")
    ggsave(glue("figs/licePerFish_ranks_{advance[i]}wk.png"), width=15, height=5)
    
    for(m in c("rmse", "rsq")) {
      map(m, 
          ~rank_results(licePerFish_wfs, rank_metric=.x, select_best=TRUE) |>
            filter(.metric==.x) |>
            select(rank, .metric, mean, model, wflow_id, .config))
      map(m,
          ~rank_results(licePerFish_wfs, rank_metric=.x, select_best=TRUE) |>
            filter(.metric==.x) |>
            select(rank, .metric, mean, model, wflow_id, .config)) |>
        saveRDS(glue("out/ensembles/licePerFish_ranks_{advance[i]}wk_{m}.rds"))
      gc()
      
      ## Best fits
      licePerFish_best_mod <- rank_results(licePerFish_wfs, rank_metric=m, select_best=TRUE)
      licePerFish_best_wf <- licePerFish_wfs |> 
        extract_workflow(licePerFish_best_mod$wflow_id[1])
      licePerFish_best_results <- licePerFish_wfs |> 
        extract_workflow_set_result(id=licePerFish_best_mod$wflow_id[1]) |>
        select_best(metric=m)
      licePerFish_final_fit <- licePerFish_best_wf |>
        finalize_workflow(licePerFish_best_results) |>
        fit(data=fit_df)
      licePerFish_best_preds <- licePerFish_wfs |> 
        extract_workflow_set_result(id=licePerFish_best_mod$wflow_id[1]) |>
        collect_predictions() |>
        filter(.config==licePerFish_best_mod$.config[1]) |>
        left_join(fold_rowNums) 
      
      saveRDS(licePerFish_best_wf, glue("out/ensembles/licePerFish_best_wf_{advance[i]}wk_{m}.rds"))
      saveRDS(licePerFish_best_results, glue("out/ensembles/licePerFish_best_results_{advance[i]}wk_{m}.rds"))
      saveRDS(licePerFish_final_fit, glue("out/ensembles/licePerFish_best_fitted_{advance[i]}wk_{m}.rds"))
      write_csv(licePerFish_best_preds, glue("out/ensembles/CV_ensFc-{advance[i]}_{m}.csv")) 
    }
    rm(licePerFish_wfs); rm(licePerFish_best_wf); rm(licePerFish_best_results)
    rm(licePerFish_best_mod); rm(licePerFish_final_fit); rm(licePerFish_best_preds)
    gc()
  }
  
  
  
  j <- 1
  if(j==2) {
    # liceBinary
    if(get_os()=="windows") {
      plan(multisession, workers=cores)
    } else {
      plan(multicore, workers=cores)
    }
    liceBinary_wfs <- workflow_set(
      preproc=map(recipes, ~.x |> update_role(lice_g05, new_role="outcome")),
      models=liceBinary_models
    ) |>
      filter(!grepl("^PCA_.*enet", wflow_id)) |>
      workflow_map("tune_grid",
                   resamples=folds,
                   grid=gridSize,
                   metrics=metric_set(roc_auc, average_precision),
                   control=control_grid(save_pred=T,
                                        save_workflow=T,
                                        parallel_over="everything",
                                        event_level="second"),
                   verbose=T)
    plan(sequential)
    cat(format(Sys.time(), "%F %T"), "  Finished liceBinary tuning, advance:", advance[i], "\n")
    autoplot(liceBinary_wfs) + scale_colour_brewer(type="qual", palette="Paired")
    ggsave(glue("figs/liceBinary_ranks_{advance[i]}wk.png"), width=15, height=5)
    
    for(m in c("roc_auc", "average_precision")) {
      map(m,
          ~rank_results(liceBinary_wfs, rank_metric=.x, select_best=TRUE) |>
            filter(.metric==.x) |>
            select(rank, .metric, mean, model, wflow_id, .config)) |>
        saveRDS(glue("out/ensembles/liceBinary_ranks_{advance[i]}wk_{m}.rds"))
      # Best fits
      liceBinary_best_mod <- rank_results(liceBinary_wfs, rank_metric=m, select_best=TRUE)
      liceBinary_best_wf <- liceBinary_wfs |>
        extract_workflow(liceBinary_best_mod$wflow_id[1])
      liceBinary_best_results <- liceBinary_wfs |>
        extract_workflow_set_result(id=liceBinary_best_mod$wflow_id[1]) |>
        select_best(metric=m)
      liceBinary_final_fit <- liceBinary_best_wf |>
        finalize_workflow(liceBinary_best_results) |>
        fit(data=fit_df)
      liceBinary_best_preds <- liceBinary_wfs |>
        extract_workflow_set_result(id=liceBinary_best_mod$wflow_id[1]) |>
        collect_predictions() |>
        filter(.config==liceBinary_best_mod$.config[1]) |>
        left_join(fold_rowNums)
      
      saveRDS(liceBinary_best_wf, glue("out/ensembles/liceBinary_best_wf_{advance[i]}wk_{m}.rds"))
      saveRDS(liceBinary_best_results, glue("out/ensembles/liceBinary_best_results_{advance[i]}wk_{m}.rds"))
      saveRDS(liceBinary_final_fit, glue("out/ensembles/liceBinary_best_fitted_{advance[i]}wk_{m}.rds"))
      write_csv(liceBinary_best_preds, glue("out/ensembles/CV_ensFc-{advance[i]}_{m}.csv"))
    }
    rm(liceBinary_wfs); rm(liceBinary_best_wf); rm(liceBinary_best_results)
    rm(liceBinary_best_mod); rm(liceBinary_final_fit); rm(liceBinary_best_preds)
    gc()
  }
  
}







# Tune ensembles: Farms ---------------------------------------------------
indices_ls <- map(unique(fit_df$sepaSiteNum),
                  ~list(analysis=which(fit_df$sepaSiteNum != .x), 
                        assessment=which(fit_df$sepaSiteNum == .x)))
splits <- lapply(indices_ls, make_splits, data=fit_df)
folds <- manual_rset(splits, paste("Split", seq_along(splits)))
fold_rowNums <- folds |>
  mutate(rowNum=map(splits, ~.x$data$rowNum[.x$out_id]),
         .row=map(splits, ~.x$out_id)) |>
  select(id, .row, rowNum) |>
  unnest(c(".row", "rowNum"))

j <- 1
if(j==1) {
  # licePerFish
  if(get_os()=="windows") {
    plan(multisession, workers=cores)
  } else {
    plan(multicore, workers=cores)
  }
  licePerFish_wfs <- workflow_set(
    preproc=map(recipes, ~.x |> update_role(licePerFish_rtrt, new_role="outcome")),
    models=licePerFish_models
  ) |>
    filter(!grepl("^PCA_.*enet", wflow_id)) |>
    workflow_map("tune_grid",
                 resamples=folds, 
                 grid=gridSize,
                 metrics=metric_set(rmse, rsq),
                 control=control_grid(save_pred=T,
                                      save_workflow=T,
                                      parallel_over="everything"),
                 verbose=T)
  plan(sequential)
  cat(format(Sys.time(), "%F %T"), "  Finished licePerFish tuning", "\n")
  autoplot(licePerFish_wfs) + scale_colour_brewer(type="qual", palette="Paired")
  ggsave(glue("figs/licePerFish_ranks_CV-farms.png"), width=15, height=5)
  
  for(m in c("rmse", "rsq")) {
    map(m, 
        ~rank_results(licePerFish_wfs, rank_metric=.x, select_best=TRUE) |>
          filter(.metric==.x) |>
          select(rank, .metric, mean, model, wflow_id, .config))
    map(m,
        ~rank_results(licePerFish_wfs, rank_metric=.x, select_best=TRUE) |>
          filter(.metric==.x) |>
          select(rank, .metric, mean, model, wflow_id, .config)) |>
      saveRDS(glue("out/ensembles/licePerFish_ranks_CV-farms_{m}.rds"))
    gc()
    
    ## Best fits
    licePerFish_best_mod <- rank_results(licePerFish_wfs, rank_metric=m, select_best=TRUE)
    licePerFish_best_wf <- licePerFish_wfs |> 
      extract_workflow(licePerFish_best_mod$wflow_id[1])
    licePerFish_best_results <- licePerFish_wfs |> 
      extract_workflow_set_result(id=licePerFish_best_mod$wflow_id[1]) |>
      select_best(metric=m)
    licePerFish_final_fit <- licePerFish_best_wf |>
      finalize_workflow(licePerFish_best_results) |>
      fit(data=fit_df)
    licePerFish_best_preds <- licePerFish_wfs |> 
      extract_workflow_set_result(id=licePerFish_best_mod$wflow_id[1]) |>
      collect_predictions() |>
      filter(.config==licePerFish_best_mod$.config[1]) |>
      left_join(fold_rowNums) 
    
    saveRDS(licePerFish_best_wf, glue("out/ensembles/licePerFish_best_wf_CV-farms_{m}.rds"))
    saveRDS(licePerFish_best_results, glue("out/ensembles/licePerFish_best_results_CV-farms_{m}.rds"))
    saveRDS(licePerFish_final_fit, glue("out/ensembles/licePerFish_best_fitted_CV-farms_{m}.rds"))
    write_csv(licePerFish_best_preds, glue("out/ensembles/CV_ensFc-CV-farms_{m}.csv")) 
  }
  rm(licePerFish_wfs); rm(licePerFish_best_wf); rm(licePerFish_best_results)
  rm(licePerFish_best_mod); rm(licePerFish_final_fit); rm(licePerFish_best_preds)
  gc()
}



j <- 2
if(j==2) {
  # liceBinary
  if(get_os()=="windows") {
    plan(multisession, workers=cores)
  } else {
    plan(multicore, workers=cores)
  }
  liceBinary_wfs <- workflow_set(
    preproc=map(recipes, ~.x |> update_role(lice_g05, new_role="outcome")),
    models=liceBinary_models
  ) |>
    filter(!grepl("^PCA_.*enet", wflow_id)) |>
    workflow_map("tune_grid",
                 resamples=folds,
                 grid=gridSize,
                 metrics=metric_set(roc_auc, average_precision),
                 control=control_grid(save_pred=T,
                                      save_workflow=T,
                                      parallel_over="everything",
                                      event_level="second"),
                 verbose=T)
  plan(sequential)
  cat(format(Sys.time(), "%F %T"), "  Finished liceBinary tuning", "\n")
  autoplot(liceBinary_wfs) + scale_colour_brewer(type="qual", palette="Paired")
  ggsave(glue("figs/liceBinary_ranks_CV-farms.png"), width=15, height=5)
  
  for(m in c("roc_auc", "average_precision")) {
    map(m,
        ~rank_results(liceBinary_wfs, rank_metric=.x, select_best=TRUE) |>
          filter(.metric==.x) |>
          select(rank, .metric, mean, model, wflow_id, .config)) |>
      saveRDS(glue("out/ensembles/liceBinary_ranks_CV-farms_{m}.rds"))
    # Best fits
    liceBinary_best_mod <- rank_results(liceBinary_wfs, rank_metric=m, select_best=TRUE)
    liceBinary_best_wf <- liceBinary_wfs |>
      extract_workflow(liceBinary_best_mod$wflow_id[1])
    liceBinary_best_results <- liceBinary_wfs |>
      extract_workflow_set_result(id=liceBinary_best_mod$wflow_id[1]) |>
      select_best(metric=m)
    liceBinary_final_fit <- liceBinary_best_wf |>
      finalize_workflow(liceBinary_best_results) |>
      fit(data=fit_df)
    liceBinary_best_preds <- liceBinary_wfs |>
      extract_workflow_set_result(id=liceBinary_best_mod$wflow_id[1]) |>
      collect_predictions() |>
      filter(.config==liceBinary_best_mod$.config[1]) |>
      left_join(fold_rowNums)
    
    saveRDS(liceBinary_best_wf, glue("out/ensembles/liceBinary_best_wf_CV-farms_{m}.rds"))
    saveRDS(liceBinary_best_results, glue("out/ensembles/liceBinary_best_results_CV-farms_{m}.rds"))
    saveRDS(liceBinary_final_fit, glue("out/ensembles/liceBinary_best_fitted_CV-farms_{m}.rds"))
    write_csv(liceBinary_best_preds, glue("out/ensembles/CV_ensFc-CV-farms_{m}.csv"))
  }
  rm(liceBinary_wfs); rm(liceBinary_best_wf); rm(liceBinary_best_results)
  rm(liceBinary_best_mod); rm(liceBinary_final_fit); rm(liceBinary_best_preds)
  gc()
}




