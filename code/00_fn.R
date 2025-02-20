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









metric_plot_base <- function(data, theme="ms") {
  
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
    geom_point(aes(type, value, colour=lab_short, shape=lab_short, size=lab_short, alpha=lab_short)) +
    scale_x_discrete(labels=label_wrap_gen(width=10)) +
    scale_colour_manual(values=c("black", "red",
                                 scico(2, begin=0.2, end=0.7, palette="broc", direction=1),
                                 scico(2, begin=0.2, end=0.7, palette="broc", direction=1),
                                 "grey50", "grey50")) +
    scale_shape_manual(values=c(1, 1, 5, 5, 1, 1, 3, 4)) +
    scale_size_manual(values=c(rep(2.5, 4), rep(1, 2), 1.5, 1.5)) +
    scale_alpha_manual(values=c(1, 1, 1, 1, 0.5, 0.5, 1, 1)) +
    facet_wrap(~metric, labeller=label_parsed) +
    plot_theme
}
