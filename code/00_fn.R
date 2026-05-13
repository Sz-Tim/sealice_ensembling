# Project: Sealice IP Ensemble
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Helper functions





#' Clean lice data from MSS
#'
#' @param f_xlsx path to .xlsx with a sheet for each year (2017-2020)
#'
#' @return compiled dataframe
#' @export
clean_mss_lice_xlsx <- function(f_xlsx) {
  library(readxl)
  list(
    # 2017
    read_xlsx(f_xlsx, sheet=1, range="A3:FD84", col_types="text",
              col_names=c("siteNo", "siteName", "businessName", "businessNo",
                          c(outer(c("count_", "mitigation_", "action_"), 
                                  1:52, 
                                  paste0)))) |>
      pivot_longer(contains("_"), names_to=c(".value", "week"), names_sep="_") |>
      mutate(weeklyAverageAf=as.numeric(count),
             week=as.numeric(week),
             weekBeginning=ymd("2017-01-01") + (week-1)*7),
    # 2018
    read_xlsx(f_xlsx, sheet=2, range="A3:FD57", col_types="text",
              col_names=c("siteNo", "siteName", "businessName", "businessNo",
                          c(outer(c("count_", "mitigation_", "action_"), 
                                  1:52, 
                                  paste0)))) |>
      pivot_longer(contains("_"), names_to=c(".value", "week"), names_sep="_") |>
      mutate(weeklyAverageAf=as.numeric(count),
             week=as.numeric(week),
             weekBeginning=ymd("2018-01-01") + (week-1)*7),
    # 2019
    read_xlsx(f_xlsx, sheet=3, range="A3:FD83", col_types="text",
              col_names=c("siteNo", "siteName", "businessName", "businessNo",
                          c(outer(c("count_", "mitigation_", "action_"), 
                                  1:52, 
                                  paste0)))) |>
      pivot_longer(contains("_"), names_to=c(".value", "week"), names_sep="_") |>
      mutate(weeklyAverageAf=as.numeric(count),
             week=as.numeric(week),
             weekBeginning=ymd("2019-01-01") + (week-1)*7),
    # 2020
    read_xlsx(f_xlsx, sheet=4, range="A3:HH78", col_types="text",
              col_names=c("siteNo", "siteName", "businessName", "businessNo",
                          c(outer(c("count_", "mitigation_", "addInfo_", "action_"), 
                                  1:53, 
                                  paste0)))) |>
      select(-starts_with("addInfo_")) |>
      pivot_longer(contains("_"), names_to=c(".value", "week"), names_sep="_") |>
      mutate(weeklyAverageAf=as.numeric(count),
             week=as.numeric(week),
             weekBeginning=ymd("2020-01-01") + (week-1)*7)
  ) |>
    reduce(bind_rows) |>
    select(siteNo, siteName, weekBeginning, weeklyAverageAf)
}






make_spline_recipe <- function(ens_df, bs_deg_free=4, IP_sims_incl=paste0("sim_0", 1:5)) {
  recipe(licePerFish_rtrt ~ .,
         data=ens_df |>
           select(rowNum, CV_k, sepaSiteNum, date, easting, northing, licePerFish_rtrt,
                  all_of(IP_sims_incl), all_of(paste0("c_", IP_sims_incl)))) |>
    step_interact(~easting:northing) |>
    step_bs(easting, northing, easting_x_northing, deg_free=bs_deg_free) |>
    update_role(c(rowNum, CV_k, date, sepaSiteNum, starts_with("c_")), new_role="id") |>
    update_role_requirements("id", bake=F) |>
    prep()
} 







make_data_rstan <- function(df, R2D2_sd=FALSE) {
  library(tidyverse)
  zeros <- df$licePerFish_rtrt==0
  sim_names <- grep("^sim_", names(df), value=T)
  dat_rstan <- list(N=nrow(df),
                    Y=df$licePerFish_rtrt,
                    N_Ye0=sum(zeros),
                    indexes_Ye0=which(zeros),
                    N_Yg0=sum(!zeros),
                    indexes_Yg0=which(!zeros),
                    K_sims=length(sim_names),
                    X=as.matrix(df |> select(all_of(sim_names))),
                    Xc=as.matrix(df |> select(all_of(paste0("c_", sim_names)))),
                    N_groups=max(df$sepaSiteNum),
                    J_group=df$sepaSiteNum,
                    sim_names=sim_names,
                    R2D2_mean_R2=0.5,
                    R2D2_prec_R2=2,
                    R2D2_cons_D2=rep(0.1, length(sim_names)*ifelse(R2D2_sd, 2, 1)),
                    R2D2_mean_R2_grp=0.5,
                    R2D2_prec_R2_grp=2,
                    R2D2_cons_D2_grp=rep(0.1, length(sim_names)),
                    R2D2_mean_R2_hu=0.5,
                    R2D2_prec_R2_hu=2,
                    R2D2_cons_D2_hu=rep(0.5, length(sim_names)),
                    prior_only=0)
  return(dat_rstan)
}


make_data_rstan_GQ <- function(df, test_rows, R2D2_sd=FALSE) {
  library(tidyverse)
  zeros <- df$licePerFish_rtrt==0
  sim_names <- grep("^sim_", names(df), value=T)
  dat_rstan <- list(orig_df=df |> mutate(test_row=row_number() %in% test_rows),
                    N=nrow(df[-test_rows,]),
                    Y=df$licePerFish_rtrt[-test_rows],
                    N_Ye0=sum(zeros[-test_rows]),
                    indexes_Ye0=which(zeros[-test_rows]),
                    N_Yg0=sum(!zeros[-test_rows]),
                    indexes_Yg0=which(!zeros[-test_rows]),
                    K_sims=length(sim_names),
                    X=as.matrix(df |> select(any_of(sim_names)))[-test_rows,, drop=F],
                    Xc=as.matrix(df |> select(any_of(paste0("c_", sim_names))))[-test_rows,, drop=F],
                    N_groups=n_distinct(df$sepaSiteNum),
                    J_group=as.numeric(as.factor(df$sepaSiteNum))[-test_rows],
                    sim_names=sim_names,
                    R2D2_mean_R2=0.5,
                    R2D2_prec_R2=2,
                    R2D2_cons_D2=rep(0.1, length(sim_names)*ifelse(R2D2_sd, 2, 1)),
                    R2D2_mean_R2_grp=0.5,
                    R2D2_prec_R2_grp=2,
                    R2D2_cons_D2_grp=rep(0.1, length(sim_names)),
                    R2D2_mean_R2_hu=0.5,
                    R2D2_prec_R2_hu=2,
                    R2D2_cons_D2_hu=rep(0.5, length(sim_names)),
                    prior_only=0,
                    GQ_rows=test_rows,
                    GQ_N=nrow(df[test_rows,]),
                    GQ_Y=df$licePerFish_rtrt[test_rows],
                    GQ_X=as.matrix(df |> select(any_of(sim_names)))[test_rows,, drop=F],
                    GQ_Xc=as.matrix(df |> select(any_of(paste0("c_", sim_names))))[test_rows,, drop=F],
                    GQ_J_group=as.numeric(as.factor(df$sepaSiteNum))[test_rows])
  return(dat_rstan)
}



make_data_rstan_AR <- function(df, test_rows, R2D2_sd=FALSE, AR=1) {
  library(tidyverse)
  zeros <- df$licePerFish_rtrt==0
  sim_names <- grep("^sim_", names(df), value=T)
  dat_rstan <- list(orig_df=df |> mutate(test_row=row_number() %in% test_rows),
                    N=nrow(df[-test_rows,]),
                    Y=df$licePerFish_rtrt[-test_rows],
                    N_Ye0=sum(zeros[-test_rows]),
                    indexes_Ye0=which(zeros[-test_rows]),
                    N_Yg0=sum(!zeros[-test_rows]),
                    indexes_Yg0=which(!zeros[-test_rows]),
                    K_sims=length(sim_names),
                    X=as.matrix(df |> select(any_of(sim_names)))[-test_rows,, drop=F],
                    Xc=as.matrix(df |> select(any_of(paste0("c_", sim_names))))[-test_rows,, drop=F],
                    N_groups=n_distinct(df$sepaSiteNum),
                    J_group=as.numeric(as.factor(df$sepaSiteNum))[-test_rows],
                    sim_names=sim_names,
                    R2D2_mean_R2=0.5,
                    R2D2_prec_R2=2,
                    R2D2_cons_D2=rep(0.1, length(sim_names)*ifelse(R2D2_sd, 2, 1)),
                    R2D2_mean_R2_grp=0.5,
                    R2D2_prec_R2_grp=2,
                    R2D2_cons_D2_grp=rep(0.1, length(sim_names)),
                    R2D2_mean_R2_hu=0.5,
                    R2D2_prec_R2_hu=2,
                    R2D2_cons_D2_hu=rep(0.5, length(sim_names)),
                    Kar=AR,
                    J_lag=df$AR_lag[-test_rows],
                    prior_only=0,
                    GQ_rows=test_rows,
                    GQ_N=nrow(df[test_rows,]),
                    GQ_Y=df$licePerFish_rtrt[test_rows],
                    GQ_X=as.matrix(df |> select(any_of(sim_names)))[test_rows,, drop=F],
                    GQ_Xc=as.matrix(df |> select(any_of(paste0("c_", sim_names))))[test_rows,, drop=F],
                    GQ_J_group=as.numeric(as.factor(df$sepaSiteNum))[test_rows])
  return(dat_rstan)
}



#' Create list of data for rstan 
#'
#' @param df Dataframe that includes columns for easting/northing/interaction 
#' splines, as created by recipes::step-bs()
#'
#' @return List of data
#' @export
make_data_rstan_sLonLat <- function(df) {
  library(tidyverse)
  zeros <- df$licePerFish_rtrt==0
  sim_names <- grep("^sim_", names(df), value=T)
  dat_rstan <- list(N=nrow(df),
                    Y=df$licePerFish_rtrt,
                    N_Ye0=sum(zeros),
                    indexes_Ye0=which(zeros),
                    N_Yg0=sum(!zeros),
                    indexes_Yg0=which(!zeros),
                    K_sims=length(sim_names),
                    X=as.matrix(df |> select(any_of(sim_names))),
                    Xc=as.matrix(df |> select(any_of(paste0("c_", sim_names)))),
                    N_groups=n_distinct(df$sepaSiteNum),
                    J_group=as.numeric(as.factor(df$sepaSiteNum)),
                    K_knots=sum(grepl("easting_bs_", names(df))),
                    Xs_easting=df |> 
                      slice_head(n=1, by=sepaSiteNum) |>
                      select(starts_with("easting_bs_")) |>
                      as.matrix(), 
                    Xs_northing=df |> 
                      slice_head(n=1, by=sepaSiteNum) |>
                      select(starts_with("northing_bs_")) |>
                      as.matrix(), 
                    Xs_easting_x_northing=df |> 
                      slice_head(n=1, by=sepaSiteNum) |>
                      select(starts_with("easting_x_northing_bs_")) |>
                      as.matrix(), 
                    sim_names=sim_names,
                    R2D2_mean_R2_hu=0.5,
                    R2D2_prec_R2_hu=2,
                    R2D2_cons_D2_hu=rep(0.5, length(sim_names)),
                    prior_only=0)
  return(dat_rstan)
}


make_data_rstan_sLonLat_GQ <- function(df, test_rows) {
  library(tidyverse)
  zeros <- df$licePerFish_rtrt==0
  sim_names <- grep("^sim_", names(df), value=T)
  dat_rstan <- list(orig_df=df |> mutate(test_row=row_number() %in% test_rows),
                    N=nrow(df[-test_rows,]),
                    Y=df$licePerFish_rtrt[-test_rows],
                    N_Ye0=sum(zeros[-test_rows]),
                    indexes_Ye0=which(zeros[-test_rows]),
                    N_Yg0=sum(!zeros[-test_rows]),
                    indexes_Yg0=which(!zeros[-test_rows]),
                    K_sims=length(sim_names),
                    X=as.matrix(df |> select(any_of(sim_names)))[-test_rows,],
                    Xc=as.matrix(df |> select(any_of(paste0("c_", sim_names))))[-test_rows,],
                    N_groups=n_distinct(df$sepaSiteNum),
                    J_group=as.numeric(as.factor(df$sepaSiteNum))[-test_rows],
                    K_knots=sum(grepl("easting_bs_", names(df))),
                    Xs_easting=df |> 
                      slice_head(n=1, by=sepaSiteNum) |>
                      select(starts_with("easting_bs_")) |>
                      as.matrix(), 
                    Xs_northing=df |> 
                      slice_head(n=1, by=sepaSiteNum) |>
                      select(starts_with("northing_bs_")) |>
                      as.matrix(), 
                    Xs_easting_x_northing=df |> 
                      slice_head(n=1, by=sepaSiteNum) |>
                      select(starts_with("easting_x_northing_bs_")) |>
                      as.matrix(), 
                    sim_names=sim_names,
                    R2D2_mean_R2_hu=0.5,
                    R2D2_prec_R2_hu=2,
                    R2D2_cons_D2_hu=rep(0.5, length(sim_names)),
                    prior_only=0,
                    GQ_rows=test_rows,
                    GQ_N=nrow(df[test_rows,]),
                    GQ_Y=df$licePerFish_rtrt[test_rows],
                    GQ_X=as.matrix(df |> select(any_of(sim_names)))[test_rows,],
                    GQ_Xc=as.matrix(df |> select(any_of(paste0("c_", sim_names))))[test_rows,],
                    GQ_J_group=as.numeric(as.factor(df$sepaSiteNum))[test_rows]
  )
  return(dat_rstan)
}



make_predictions_ensBlend <- function(out, newdata, iter=2000, seed=NULL, mode="epred", re=TRUE, re_hu=FALSE) {
  library(tidyverse); library(rstan)
  # hurdle component is fitted with centered IP
  
  set.seed(ifelse(is.null(seed), runif(1), seed))
  b_b0 <- rstan::extract(out, pars="b_b0")[[1]]
  b_IP <- rstan::extract(out, pars="b_IP")[[1]]
  Intercept_hu <- rstan::extract(out, pars="Intercept_hu")[[1]]
  b_hu <- rstan::extract(out, pars="b_hu")[[1]]
  sigma <- rstan::extract(out, pars="sigma")[[1]]
  
  iters <- sample.int(nrow(b_b0), iter, replace=T)
  preds <- matrix(0, nrow=iter, ncol=nrow(newdata))
  dat_ls <- make_data_rstan(newdata)
  ensIP <- matrix(0, nrow=iter, ncol=nrow(newdata))
  hu_RE <- matrix(0, nrow=iter, ncol=nrow(newdata))
  
  if(re) {
    r_grp <- rstan::extract(out, pars="r_grp")[[1]]
    for(i in 1:ncol(preds)) {
      ensIP[,i] <- dat_ls$X[i,,drop=F] %*% t(r_grp[iters, dat_ls$J_group[i],,drop=T])
    }
  } else {
    b_p <- rstan::extract(out, pars="b_p")[[1]]
    for(i in 1:ncol(preds)) {
      ensIP[,i] <- dat_ls$X[i,,drop=F] %*% t(b_p[iters,,drop=F])
    }
  }
  if(re_hu) {
    r_grp_hu <- rstan::extract(out, pars="r_grp_hu")[[1]]
    for(i in 1:ncol(preds)) {
      hu_RE[,i] <- cbind(1, dat_ls$Xc[i,,drop=F]) %*% t(r_grp_hu[iters, dat_ls$J_group[i],])
    }
  }
  
  for(i in 1:ncol(preds)) {
    mu <- b_b0[iters] + b_IP[iters] * c(ensIP[,i])
    hu <- Intercept_hu[iters] + c(dat_ls$Xc[i,,drop=F] %*% t(b_hu[iters,,drop=F])) + hu_RE[,i]
    # hu = pr(0)
    if(mode=="epred") {
      preds[,i] <- (1-brms::inv_logit_scaled(hu)) * exp(mu)
    } else if(mode=="predict") {
      preds[,i] <- (1 - rbinom(iter, 1, brms::inv_logit_scaled(hu))) * 
        rlnorm(iter, mu, sigma[iters]) 
    }
  }
  return(preds)
}



make_predictions_ensBlend_sLonLat <- function(out, newdata, iter=2000, seed=NULL, mode="point_epred") {
  library(tidyverse); library(rstan)
  # hurdle component is fitted with centered IP
  
  set.seed(ifelse(is.null(seed), runif(1), seed))
  b_b0 <- rstan::extract(out, pars="b_b0")[[1]]
  b_IP <- rstan::extract(out, pars="b_IP")[[1]]
  b_s_easting <- rstan::extract(out, pars="b_s_easting")[[1]]
  b_s_northing <- rstan::extract(out, pars="b_s_northing")[[1]]
  b_s_easting_x_northing <- rstan::extract(out, pars="b_s_easting_x_northing")[[1]]
  Intercept_hu <- rstan::extract(out, pars="Intercept_hu")[[1]]
  b_hu <- rstan::extract(out, pars="b_hu")[[1]]
  sigma <- rstan::extract(out, pars="sigma")[[1]]
  
  iters <- sample.int(nrow(b_b0), iter, replace=T)
  preds <- matrix(0, nrow=iter, ncol=nrow(newdata))
  dat_ls <- make_data_rstan_sLonLat(newdata)
  ensIP <- matrix(0, nrow=iter, ncol=nrow(newdata))
  hu_RE <- matrix(0, nrow=iter, ncol=nrow(newdata))
  
  # ARE GROUPS ORDERED CORRECTLY???
  # b_p_uc <- b_p <- array(dim=list(nrow(newdata), iter, dim(b_hu)[2]))
  b_p_uc <- b_p <- array(dim=list(dat_ls$N_groups, iter, dim(b_hu)[2]))
  for(k in 1:dim(b_s_easting)[3]) {
    b_p_uc[,,k] <- dat_ls$Xs_easting %*% t(b_s_easting[iters,,k]) + 
      dat_ls$Xs_northing %*% t(b_s_northing[iters,,k]) +
      dat_ls$Xs_easting_x_northing %*% t(b_s_easting_x_northing[iters,,k])
  }
  for(i in 1:dim(b_p_uc)[1]) {
    for(j in 1:dim(b_p_uc)[2]) {
      b_p[i,j,] <- exp(b_p_uc[i,j,])/sum(exp(b_p_uc[i,j,]))
    }
  }
  if(mode=="b_p") {
    return(b_p)
  }
  
  for(i in 1:ncol(preds)) {
    ensIP[,i] <- dat_ls$X[i,,drop=F] %*% t(b_p[dat_ls$J_group[i],,])
  }
  if(mode=="ensIP") {
    return(ensIP)
  }
  
  for(i in 1:ncol(preds)) {
    mu <- b_b0[iters] + b_IP[iters] * c(ensIP[,i])
    hu <- Intercept_hu[iters] + c(dat_ls$Xc[i,,drop=F] %*% t(b_hu[iters,,drop=F])) + hu_RE[,i]
    # hu = pr(0)
    if(mode=="point_epred") {
      preds[,i] <- (1-brms::inv_logit_scaled(hu)) * exp(mu)
    } else if(mode=="point_predict") {
      preds[,i] <- (1 - rbinom(iter, 1, brms::inv_logit_scaled(hu))) * 
        rlnorm(iter, mu, sigma[iters]) 
    } 
  }
  return(preds)
}



make_predictions_ensBlend_sLonLat_RE <- function(out, newdata, iter=2000, seed=NULL, mode="point_epred", type="site") {
  library(tidyverse); library(rstan)
  # hurdle component is fitted with centered IP
  
  set.seed(ifelse(is.null(seed), runif(1), seed))
  b_b0 <- rstan::extract(out, pars="b_b0")[[1]]
  b_IP <- rstan::extract(out, pars="b_IP")[[1]]
  Intercept_hu <- rstan::extract(out, pars="Intercept_hu")[[1]]
  b_hu <- rstan::extract(out, pars="b_hu")[[1]]
  sigma <- rstan::extract(out, pars="sigma")[[1]]

  if(type=="site") {
    b_p <- rstan::extract(out, pars="b_p")[[1]]
  } else if(type=="surface") {
    b_s_easting <- rstan::extract(out, pars="b_s_easting")[[1]]
    b_s_northing <- rstan::extract(out, pars="b_s_northing")[[1]]
    b_s_easting_x_northing <- rstan::extract(out, pars="b_s_easting_x_northing")[[1]]
  }
    
  iters <- sample.int(nrow(b_b0), iter, replace=T)
  preds <- matrix(0, nrow=iter, ncol=nrow(newdata))
  dat_ls <- make_data_rstan_sLonLat(newdata)
  ensIP <- matrix(0, nrow=iter, ncol=nrow(newdata))
  hu_RE <- matrix(0, nrow=iter, ncol=nrow(newdata))
  
  if(type=="surface") {
    b_p_uc <- b_p <- array(dim=list(dat_ls$N_groups, iter, dim(b_hu)[2]))
    for(k in 1:dim(b_s_easting)[3]) {
      b_p_uc[,,k] <- dat_ls$Xs_easting %*% t(b_s_easting[iters,,k]) + 
        dat_ls$Xs_northing %*% t(b_s_northing[iters,,k]) +
        dat_ls$Xs_easting_x_northing %*% t(b_s_easting_x_northing[iters,,k])
    }
    for(i in 1:dim(b_p_uc)[1]) {
      for(j in 1:dim(b_p_uc)[2]) {
        b_p[i,j,] <- exp(b_p_uc[i,j,])/sum(exp(b_p_uc[i,j,]))
      }
    }
  }
  if(mode=="b_p") {
    return(b_p)
  }
  
  if(type=="site") {
    for(i in 1:ncol(preds)) {
      ensIP[,i] <- dat_ls$X[i,,drop=F] %*% t(b_p[iters,dat_ls$J_group[i],])
    }
  } else if(type=="surface") {
    for(i in 1:ncol(preds)) {
      ensIP[,i] <- dat_ls$X[i,,drop=F] %*% t(b_p[dat_ls$J_group[i],,])
    }
  }
  if(mode=="ensIP") {
    return(ensIP)
  }
  
  for(i in 1:ncol(preds)) {
    mu <- b_b0[iters] + b_IP[iters] * c(ensIP[,i])
    hu <- Intercept_hu[iters] + c(dat_ls$Xc[i,,drop=F] %*% t(b_hu[iters,,drop=F])) + hu_RE[,i]
    # hu = pr(0)
    if(mode=="point_epred") {
      preds[,i] <- (1-brms::inv_logit_scaled(hu)) * exp(mu)
    } else if(mode=="point_predict") {
      preds[,i] <- (1 - rbinom(iter, 1, brms::inv_logit_scaled(hu))) * 
        rlnorm(iter, mu, sigma[iters]) 
    } 
  }
  return(preds)
}





make_predictions_candidate <- function(out, newdata, sim, iter=2000, seed=NULL, re=FALSE) {
  library(tidyverse); library(rstan)
  
  newdata <- newdata |> 
    select(sepaSite, sepaSiteNum, date, licePerFish_rtrt, matches(sim))
  
  set.seed(ifelse(is.null(seed), runif(1), seed))
  b_b0 <- rstan::extract(out, pars="b_b0")[[1]]
  b_IP <- rstan::extract(out, pars="b_IP")[[1]]
  Intercept_hu <- rstan::extract(out, pars="Intercept_hu")[[1]]
  b_hu <- rstan::extract(out, pars="b_hu")[[1]]
  sigma <- rstan::extract(out, pars="sigma")[[1]]
  
  iters <- sample.int(nrow(b_b0), iter, replace=T)
  preds <- matrix(0, nrow=iter, ncol=nrow(newdata))
  dat_ls <- make_data_rstan(newdata)
  hu_RE <- matrix(0, nrow=iter, ncol=nrow(newdata))
  
  if(re) {
    r_grp_hu <- rstan::extract(out, pars="r_grp_hu")[[1]]
    for(i in 1:ncol(preds)) {
      hu_RE[,i] <- cbind(1, dat_ls$Xc[i,,drop=F]) %*% t(r_grp_hu[iters, dat_ls$J_group[i],])
    }
  } 
  
  for(i in 1:ncol(preds)) {
    mu <- b_b0[iters] + b_IP[iters] * dat_ls$X[i,]
    hu <- Intercept_hu[iters] + b_hu[iters,] * dat_ls$Xc[i,] + hu_RE[,i]
    # hu = pr(0)
    preds[,i] <- (1-brms::inv_logit_scaled(hu)) * exp(mu)
  }
  return(preds)
}


ens_parallel <- function(z_df_i, z_mx, p_dir) {
  prog <- progressor(along=1:max(z_df_i$ii))
  foreach(k=1:max(z_df_i$ii), .combine=rbind, .inorder=TRUE, 
          .options.future=list(globals=structure(TRUE, add=c("z_mx", "p_dir", "z_df_i")))) %dofuture% {
            prog()
            rows <- which(z_df_i$ii==k)
            p <- readRDS(glue("{p_dir}i_{as.integer(z_df_i$i[rows[1]])}.rds"))
            ens_k <- z_mx[rows,,drop=F] %*% p
            t(apply(ens_k, 1, function(x) c(mean(x), sd(x))))
          }
}




join_candidates <- function(data, candidate_df, cand_id) {
  cand_id_col <- paste0("cand", cand_id)
  cand_pred_col <- paste0("cand_pred", cand_id)
  data |>
    select(all_of(c("rowNum", "ens_id", cand_id_col))) |>
    set_names("rowNum", "ens_id", "s") |>
    left_join(candidate_df, by=join_by(s, rowNum)) |>
    rename_with(~paste0(cand_pred_col), all_of("cand_pred")) |>
    select(all_of(c("rowNum", "ens_id", cand_pred_col)))
}

join_candidates_siteNum <- function(data, candidate_df, cand_id) {
  cand_id_col <- paste0("cand", cand_id)
  cand_pred_col <- paste0("cand_pred", cand_id)
  data |>
    select(all_of(c("sepaSiteNum", "ens_id", cand_id_col))) |>
    set_names("sepaSiteNum", "ens_id", "s") |>
    left_join(candidate_df, by=join_by(s, sepaSiteNum)) |>
    rename_with(~paste0(cand_pred_col), all_of("cand_pred")) |>
    select(all_of(c("sepaSiteNum", "ens_id", cand_pred_col)))
}

join_candidates_date <- function(data, candidate_df, cand_id) {
  cand_id_col <- paste0("cand", cand_id)
  cand_pred_col <- paste0("cand_pred", cand_id)
  data |>
    select(all_of(c("date", "ens_id", cand_id_col))) |>
    set_names("date", "ens_id", "s") |>
    left_join(candidate_df, by=join_by(s, date)) |>
    rename_with(~paste0(cand_pred_col), all_of("cand_pred")) |>
    select(all_of(c("date", "ens_id", cand_pred_col)))
}








metric_plot_base <- function(data, theme="ms", colours) {
  
  if(theme=="ms") {
    plot_theme <- theme(panel.grid.major.y=element_line(colour="grey85", linewidth=0.4),
                        panel.grid.minor.y=element_line(colour="grey90", linewidth=0.2),
                        axis.title.x=element_blank(),
                        axis.title.y=element_text(size=9),
                        axis.text.x=element_text(vjust=1, size=8),
                        axis.text.y=element_text(size=8),
                        legend.position="none")
  } else if(theme=="talk") {
    plot_theme <- theme(panel.grid.major.y=element_line(colour="grey85", linewidth=0.4),
                        panel.grid.minor.y=element_line(colour="grey90", linewidth=0.2),
                        axis.title.x=element_blank(),
                        axis.title.y=element_text(size=9),
                        axis.text.x=element_text(vjust=1, size=8),
                        axis.text.y=element_text(size=8),
                        legend.position="none")
  } else {
    plot_theme <- theme_classic()
  }
  
  ggplot(data=data) + 
    geom_point(aes(type, value, colour=lab_short, shape=lab_short, size=lab_short, alpha=lab_short), stroke=0.7) +
    geom_text(data=data |> filter(sim=="sim_07"), 
              aes(type, value, colour=lab_short, label=paste0(lab, "%->% ''")), parse=T,
              hjust=1, size=1.5, nudge_x=-0.1) +
    scale_colour_manual(values=colours) +
    scale_shape_manual(values=c(1, 1, 1, 4, 3)) +
    scale_size_manual(values=c(rep(2.5, 3), rep(1, 2))) +
    scale_alpha_manual(values=c(1, 1, 1, 1, 1)) +
    plot_theme
}


metric_plot_AEIPsens_base <- function(data, theme="ms", colours) {
  
  if(theme=="ms") {
    plot_theme <- theme(panel.grid.major.y=element_line(colour="grey85", linewidth=0.4),
                        panel.grid.minor.y=element_line(colour="grey90", linewidth=0.2),
                        axis.title.x=element_blank(),
                        axis.title.y=element_text(size=9),
                        axis.text.x=element_text(vjust=1, size=8),
                        axis.text.y=element_text(size=8),
                        legend.position="none")
  } else if(theme=="talk") {
    plot_theme <- theme(panel.grid.major.y=element_line(colour="grey85", linewidth=0.4),
                        panel.grid.minor.y=element_line(colour="grey90", linewidth=0.2),
                        axis.title.x=element_blank(),
                        axis.title.y=element_text(size=9),
                        axis.text.x=element_text(vjust=1, size=8),
                        axis.text.y=element_text(size=8),
                        legend.position="none")
  } else {
    plot_theme <- theme_classic()
  }
  
  # adjusted for no EnsFcst, 'run' on x-axis; multi-panel with 'type'~'metric``
  ggplot(data=data) + 
    geom_point(aes(run, value, colour=lab_short, shape=lab_short, size=lab_short), stroke=0.7) +
    geom_line(aes(run, value, colour=lab_short, linewidth=lab_short, group=sim)) +
    scale_colour_manual(values=colours) +
    scale_shape_manual(values=c(1, 1, 1, 4, 3) |> set_names(names(colours)[c(1:3,8:7)])) +
    scale_size_manual(values=c(rep(2.5, 3), rep(1, 2)) |> set_names(names(colours)[c(1:3,7:8)])) +
    scale_linewidth_manual(values=c(rep(0.5, 3), rep(0.1, 2)) |> set_names(names(colours)[c(1:3,7:8)])) +
    plot_theme
}



metric_plot_TRTsens_base <- function(data, theme="ms", colours) {
  
  if(theme=="ms") {
    plot_theme <- theme(panel.grid.major.y=element_line(colour="grey85", linewidth=0.3),
                        panel.grid.minor.y=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_text(size=9),
                        axis.text.x=element_text(vjust=1, size=8),
                        axis.text.y=element_text(size=8),
                        legend.position="none")
  } else if(theme=="talk") {
    plot_theme <- theme(panel.grid.major.y=element_line(colour="grey85", linewidth=0.3),
                        panel.grid.minor.y=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_text(size=9),
                        axis.text.x=element_text(vjust=1, size=8),
                        axis.text.y=element_text(size=8),
                        legend.position="none")
  } else {
    plot_theme <- theme_classic()
  }
  
  # adjusted for no EnsFcst, 'run' on x-axis; multi-panel with 'type'~'metric``
  ggplot(data=data) + 
    geom_point(aes(dataset, value, colour=lab_short, shape=lab_short, size=lab_short), stroke=0.7) +
    geom_line(aes(dataset, value, colour=lab_short, linewidth=lab_short, group=sim)) +
    scale_colour_manual(values=colours) +
    scale_shape_manual(values=c(1, 1, 1, 4, 3) |> set_names(names(colours)[c(1:3,8:7)])) +
    scale_size_manual(values=c(rep(2.5, 3), rep(1, 2)) |> set_names(names(colours)[c(1:3,7:8)])) +
    scale_linewidth_manual(values=c(rep(0.5, 3), rep(0.1, 2)) |> set_names(names(colours)[c(1:3,7:8)])) +
    plot_theme
}



summarise_param_posterior <- function(post_df, sim_key, param) {
  post_df |>
    inner_join(sim_key |> 
                 select(Simulation, any_of(param)) |>
                 rename_with(~"this_param", .cols=2), 
               by=join_by(Simulation)) |>
    drop_na() |>
    group_by(iter, rowNum) |>
    mutate(p=p/sum(p)) |>
    summarise(param_post=sum(p*this_param)) |>
    ungroup() |>
    set_names(c("iter", "rowNum", param))
}




get_param_scale2 <- function(param, sim_key, type="c") {
  library(scales)
  zero_to_one_params <- c("eggTemp_fn", "mortSal_fn", "fixDepth")
  swim_params <- c("swimUpSpeedMean", "swimDownSpeedMean")
  salinity_params <- c("salinityThreshMin", "salinityThreshMax")
  light_params <- c("lightThreshNauplius", "lightThreshCopepodid")
  diffusion_params <- c("D_h", "D_hVert")
  n_breaks <- 5
  lims <- range(sim_key[[param]], na.rm=T)
  if(param %in% zero_to_one_params) {
    if(grepl("fn", param)) {
      pal <- scale_fill_viridis_c(expression('Function'), option="turbo", limits=c(0, 1),
                                  breaks=seq(0, 1, by=0.1), labels=c("constant", rep("", 9), "logistic"))
    } else {
      pal <- scale_fill_viridis_c(expression('Dimensions'), option="turbo", limits=c(0, 1),
                                  breaks=seq(0, 1, by=0.1), labels=c("3D", rep("", 9), "2D"))
    }
  }
  if(param %in% swim_params) {
    pal <- scale_fill_viridis_c(expression(cm %.% s^-1), option="mako", 
                                limits=lims,
                                breaks=breaks_extended(n_breaks))#,
                                # labels = label_scientific(digits = 3)) 
  }
  if(param %in% salinity_params) {
    pal <- scale_fill_distiller("psu", palette="Blues", 
                                limits=lims,
                                breaks=breaks_extended(n_breaks)) 
  }
  if(param %in% light_params) {
    pal <- scale_fill_viridis_c(expression(mu*'mol' %.%~'m'^-2 %.% s^-1), 
                                limits=lims,
                                breaks=breaks_extended(n_breaks), 
                                labels = label_scientific(digits = 3)) 
  }
  if(param %in% diffusion_params) {
    pal <- scale_fill_viridis_c(expression(m^2 %.% s^-2), option="cividis", 
                                breaks=breaks_extended(n_breaks),
                                limits=log10(range(sim_key[[param]][sim_key[[param]]>0], na.rm=T)),
                                labels=label_math(expr = 10^.x, format = force))
  } 
  if(param == "viableDegreeDays") {
    pal <- scale_fill_viridis_c(expression(degree*C %.% d), option="rocket", 
                                limits=lims,
                                breaks=breaks_extended(n_breaks))
  }
  return(pal)
}




make_param_map_plot <- function(param_df, site_i, mesh_land, sim_key) {
  this_param <- param_df$var[1]
  if(grepl("D_h", this_param)) {
    param_df$post_mn <- log10(param_df$post_mn)
  }
  this_scale <- get_param_scale(this_param, sim_key)
  param_df |>
    ggplot() + 
    geom_raster(aes(easting, northing, fill=post_mn)) +
    stat_contour(aes(easting, northing, z=post_mn), colour="white", linewidth=0.1, breaks=this_scale$breaks) +
    geom_sf(data=mesh_land, fill="grey90", colour="grey40", linewidth=0.1) +
    geom_point(data=site_i, aes(easting, northing), shape=1, colour="black", size=0.4) +
    this_scale + 
    scale_x_continuous(limits=range(site_i$easting)*c(0.96, 1), oob=scales::oob_keep) +
    scale_y_continuous(limits=range(site_i$northing), oob=scales::oob_keep) +
    facet_wrap(~var_pretty, labeller=label_wrap_gen(19)) +
    theme(axis.text=element_blank(),
          axis.title=element_blank(),
          axis.ticks=element_blank(),
          legend.position="bottom",
          legend.title.position="top",
          legend.title=element_text(size=9, hjust=0.5),
          legend.box.margin=margin(0,0,0,0),
          legend.margin=margin(0,0,0,0),
          legend.key.height=unit(0.2, "cm"),
          legend.key.width=unit(0.8, "cm"),
          legend.text=element_text(size=6),
          panel.spacing=unit(0, 'cm'))
}





plot_metric_ordered <- function(df, m, highlight="sim_04") {
  df_highlight <- df |>
    filter(sim %in% highlight)
  df |>
    arrange({{m}}) |>
    mutate(sim=factor(sim, levels=unique(sim))) |>
    ggplot(aes({{m}}, sim)) + 
    geom_vline(data=df_highlight, aes(xintercept={{m}}, colour=sim), linetype=2) +
    geom_hline(data=df_highlight, aes(yintercept=sim, colour=sim), linetype=2) +
    geom_point() + 
    theme(legend.position="none")
}





# TODO: Hacky -- remove... solving an issue with biotrackR versions.
set_biotracker_properties2 <- function (properties_file_path = NULL, coordOS = "true", mesh1 = "/home/sa04ts/hydro/meshes/WeStCOMS2_mesh.nc", 
                                        mesh1Type = "FVCOM", mesh1Domain = "westcoms2", datadir = "/home/sa04ts/hydro/WeStCOMS2/Archive/", 
                                        datadirPrefix = "netcdf_", datadirSuffix = "", mesh2 = "", 
                                        mesh2Type = "FVCOM", mesh2Domain = "westcoms2", datadir2 = "/home/sa04ts/hydro/WeStCOMS2/Archive/", 
                                        datadir2Prefix = "netcdf_", datadir2Suffix = "", sitefile = "../../data/farm_sites.csv", 
                                        sitefileEnd = "../../data/farm_sites.csv", siteDensityPath = "", 
                                        daylightPath = "", verboseSetUp = "false", start_ymd = 20190401, 
                                        numberOfDays = 7, checkOpenBoundaries = "false", openBoundaryThresh = 500, 
                                        duplicateLastDay = "true", recordsPerFile1 = 25, dt = 3600, 
                                        maxDepth = 10000, parallelThreads = 4, parallelThreadsHD = 1, 
                                        releaseScenario = 1, releaseInterval = 1, nparts = 1, setStartDepth = "true", 
                                        startDepth = 1, fixDepth = "false", stepsPerStep = 30, variableDh = "false", 
                                        variableDhV = "false", D_h = 0.1, D_hVert = 0.001, salinityThreshMin = 23, 
                                        salinityThreshMax = 31, swimLightLevel = "true", lightThreshCopepodid = 2.06e-05, 
                                        lightThreshNauplius = 0.392, swimUpSpeedMean = NULL, swimUpSpeedStd = NULL, 
                                        swimUpSpeedCopepodidMean = -5e-04, swimUpSpeedCopepodidStd = 1e-04, 
                                        swimUpSpeedNaupliusMean = -0.00025, swimUpSpeedNaupliusStd = 5e-05, 
                                        swimDownSpeedMean = NULL, swimDownSpeedStd = NULL, swimDownSpeedCopepodidMean = 0.001, 
                                        swimDownSpeedCopepodidStd = 2e-04, swimDownSpeedNaupliusMean = 0.001, 
                                        swimDownSpeedNaupliusStd = 2e-04, passiveSinkingIntercept = 0.001527, 
                                        passiveSinkingSlope = -1.68e-05, eggTemp_fn = "constant", 
                                        eggTemp_b = "", mortSal_fn = "constant", mortSal_b = "", 
                                        viabletime = -1, maxParticleAge = -1, viableDegreeDays = 40, 
                                        maxDegreeDays = 150, recordImmature = "false", recordPsteps = "false", 
                                        splitPsteps = "false", pstepsInterval = 168, pstepsMaxDepth = 10000, 
                                        recordVertDistr = "false", vertDistrInterval = 1, vertDistrMax = 20, 
                                        recordMovement = "false", recordElemActivity = "false", recordConnectivity = "true", 
                                        connectImmature = "false", connectDepth1_min = 0, connectDepth1_max = 10000, 
                                        connectDepth2_min = 10000, connectDepth2_max = 10000, connectivityInterval = 24, 
                                        connectivityThresh = 100, recordLocations = "false", recordArrivals = "false") 
{
  params <- c(coordOS = coordOS, mesh1 = mesh1, mesh1Type = mesh1Type, 
              mesh1Domain = mesh1Domain, datadir = datadir, datadirPrefix = datadirPrefix, 
              datadirSuffix = datadirSuffix, mesh2 = mesh2, mesh2Type = mesh2Type, 
              mesh2Domain = mesh2Domain, datadir2 = datadir2, datadir2Prefix = datadir2Prefix, 
              datadir2Suffix = datadir2Suffix, sitefile = sitefile, 
              sitefileEnd = sitefileEnd, siteDensityPath = siteDensityPath, 
              daylightPath = daylightPath, verboseSetUp = verboseSetUp, 
              start_ymd = start_ymd, numberOfDays = numberOfDays, checkOpenBoundaries = checkOpenBoundaries, 
              openBoundaryThresh = openBoundaryThresh, duplicateLastDay = duplicateLastDay, 
              recordsPerFile1 = recordsPerFile1, dt = dt, maxDepth = maxDepth, 
              parallelThreads = parallelThreads, parallelThreadsHD = parallelThreadsHD, 
              releaseScenario = releaseScenario, releaseInterval = releaseInterval, 
              nparts = nparts, setStartDepth = setStartDepth, startDepth = startDepth, 
              fixDepth = fixDepth, stepsPerStep = stepsPerStep, variableDh = variableDh, 
              variableDhV = variableDhV, D_h = D_h, D_hVert = D_hVert, 
              salinityThreshMin = salinityThreshMin, salinityThreshMax = salinityThreshMax, 
              swimLightLevel = swimLightLevel, lightThreshCopepodid = lightThreshCopepodid, 
              lightThreshNauplius = lightThreshNauplius, swimUpSpeedMean = swimUpSpeedMean, 
              swimUpSpeedStd = swimUpSpeedStd, swimUpSpeedCopepodidMean = ifelse(is.null(swimUpSpeedMean), 
                                                                                 swimUpSpeedCopepodidMean, swimUpSpeedMean), swimUpSpeedCopepodidStd = ifelse(is.null(swimUpSpeedStd), 
                                                                                                                                                              swimUpSpeedCopepodidStd, swimUpSpeedStd), swimUpSpeedNaupliusMean = ifelse(is.null(swimUpSpeedMean), 
                                                                                                                                                                                                                                         swimUpSpeedNaupliusMean, swimUpSpeedMean), swimUpSpeedNaupliusStd = ifelse(is.null(swimUpSpeedStd), 
                                                                                                                                                                                                                                                                                                                    swimUpSpeedNaupliusStd, swimUpSpeedStd), swimDownSpeedMean = swimDownSpeedMean, 
              swimDownSpeedStd = swimDownSpeedStd, swimDownSpeedCopepodidMean = ifelse(is.null(swimDownSpeedMean), 
                                                                                       swimDownSpeedCopepodidMean, swimDownSpeedMean), swimDownSpeedCopepodidStd = ifelse(is.null(swimDownSpeedStd), 
                                                                                                                                                                          swimDownSpeedCopepodidStd, swimDownSpeedStd), swimDownSpeedNaupliusMean = ifelse(is.null(swimDownSpeedMean), 
                                                                                                                                                                                                                                                           swimDownSpeedNaupliusMean, swimDownSpeedMean), swimDownSpeedNaupliusStd = ifelse(is.null(swimDownSpeedStd), 
                                                                                                                                                                                                                                                                                                                                            swimDownSpeedNaupliusStd, swimDownSpeedStd), passiveSinkingIntercept = passiveSinkingIntercept, 
              passiveSinkingSlope = passiveSinkingSlope, eggTemp_fn = eggTemp_fn, 
              eggTemp_b = eggTemp_b, mortSal_fn = mortSal_fn, mortSal_b = mortSal_b, 
              viabletime = viabletime, maxParticleAge = maxParticleAge, 
              viableDegreeDays = viableDegreeDays, maxDegreeDays = maxDegreeDays, 
              recordImmature = recordImmature, recordPsteps = recordPsteps, 
              splitPsteps = splitPsteps, pstepsInterval = pstepsInterval, 
              pstepsMaxDepth = pstepsMaxDepth, recordVertDistr = recordVertDistr, 
              vertDistrInterval = vertDistrInterval, vertDistrMax = vertDistrMax, 
              recordMovement = recordMovement, recordElemActivity = recordElemActivity, 
              recordConnectivity = recordConnectivity, connectImmature = connectImmature, 
              connectDepth1_min = connectDepth1_min, connectDepth1_max = connectDepth1_max, 
              connectDepth2_min = connectDepth2_min, connectDepth2_max = connectDepth2_max, 
              connectivityInterval = connectivityInterval, connectivityThresh = connectivityThresh, 
              recordLocations = recordLocations, recordArrivals = recordArrivals)
  if (params["eggTemp_b"] == "") params <- params[-which(names(params) == "eggTemp_b")]
  if(params["mortSal_b"]=="") params <- params[-which(names(params)=="mortSal_b")]
  properties_out <- paste(names(params), params, sep = "=", 
                          collapse = "\n")
  if (!is.null(properties_file_path)) {
    cat(str_replace_all(str_replace_all(properties_out, "\\\\", 
                                        "\\\\\\\\"), "\\ ", "\\\\\\\\ "), "\n", file = properties_file_path)
  }
  return(properties_out)
}

