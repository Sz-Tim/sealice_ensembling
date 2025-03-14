# Project: Sealice IP Ensemble
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Compile output from long runs

# These simulations are intended for hourly visualizations.
# Output files are created for each hour, with a column for each simulation and 
# a row for each occupied WeStCOMS element.



# setup -------------------------------------------------------------------

library(tidyverse); library(glue)
library(doFuture)
library(furrr)
library(carrier)
library(sf)
library(rstan)
library(recipes)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
library(biotrackR) # devtools::install_github("Sz-Tim/biotrackR")
source("code/00_fn.R")
theme_set(theme_bw() + theme(panel.grid=element_blank()))

out_dir <- "out/sim_2023-MarMay"
p_dir <- "out/ensembles/p_meshCentroids/"
sim_i <- read_csv(glue("{out_dir}/sim_i.csv")) |>
  mutate(sim=paste0("sim_", i))




# extract output ----------------------------------------------------------

# extract biotracker output: sim_[0-9][0-9].tar.gz
if(FALSE) {
  f_tgz <- dirf(out_dir, "tar.gz")
  walk(f_tgz, ~untar(.x, exdir=str_remove(.x, ".tar.gz")))
}




# ensemble weights --------------------------------------------------------

mod <- "n20_sLonLatD3"
out_ensBlend <- readRDS(glue("out/ensembles/ensBlend_{mod}_FULL_stanfit.rds"))
dat_ensBlend <- readRDS(glue("out/ensembles/ensBlend_{mod}_FULL_standata.rds"))
meshCentroid_df <- st_read("data/WeStCOMS2_mesh.gpkg") |>
  st_centroid() |>
  sevcheck::add_lonlat(drop_geom=T) |>
  as_tibble() |>
  select(i, lon, lat) |>
  rename(sepaSiteNum=i, easting=lon, northing=lat) |>
  bind_cols(map_dfc(dat_ensBlend$sim_names, ~tibble(0) |> setNames(.x))) |>
  bind_cols(map_dfc(dat_ensBlend$sim_names, ~tibble(0) |> setNames(paste0("c_", .x)))) |>
  mutate(rowNum=1,
         date=ymd("2023-01-01"),
         CV_k=1,
         licePerFish_rtrt=1)
ensBlend_rec <- make_spline_recipe(meshCentroid_df, 3, sim_i$sim[1:20])
# plan(multisession, workers=18)
# foreach(i=1:nrow(meshCentroid_df),
#         .options.future=list(globals=structure(TRUE, add="p_dir"))) %dofuture% {
for(i in 1:nrow(meshCentroid_df)) {
  t(make_predictions_ensBlend_sLonLat(out_ensBlend,
                                      newdata=bake(ensBlend_rec, meshCentroid_df[i,]),
                                      iter=2000, seed=1003, mode="b_p")[1,,]) |>
    saveRDS(glue("{p_dir}/i_{meshCentroid_df$sepaSiteNum[i]}.rds"))
}
plan(sequential)






# particle densities ------------------------------------------------------

plan(multisession, workers=18)
date_seq <- seq(ymd("2023-04-01"), ymd("2023-05-31"), by=1) |> str_remove_all("-")
ps_lims <- tibble(ens_mn=c(0,0),
                  ens_sd=c(0,0),
                  ens_CL005=c(0,0),
                  ens_CL995=c(0,0),
                  ens_CI99width=c(0,0))
for(i in 1:length(date_seq)) {
  ps_i <- load_psteps_simSets(out_dir, 
                              st_read("data/WeStCOMS2_mesh.gpkg") |>
                                mutate(vol_top50m=area*pmin(depth, 50)) |>
                                st_drop_geometry() |> 
                                select(i, area, vol_top50m), 
                              sim_i, ncores=1, liceScale=1, 
                              stage=paste0("Mature_", date_seq[i]), per_m2=TRUE, trans="4th_rt")
  i_ts <- grep("^t_", names(ps_i), value=T)
  cat("Starting", date_seq[i])
  for(j in seq_along(i_ts)) {
    ps_j <- ps_i |> 
      filter(sim %in% dat_ensMixRE$sim_names) |>
      filter(!is.na(i)) |>
      select("sim", "i", all_of(i_ts[j])) |>
      drop_na() |>
      pivot_wider(names_from=sim, values_from=starts_with("t_")) 
    if(any(!is.na(ps_j$i))) {
      ps_j_mx <- ps_j |> select(all_of(dat_ensMixRE$sim_names)) |> as.matrix()
      ps_j_mx[is.na(ps_j_mx)] <- 0
      
      # calculate ensIP in parallel 
      ensIP <- foreach(k=1:nrow(ps_j), .combine=rbind, .inorder=TRUE, 
                       .options.future=list(globals=structure(TRUE, add=c("ps_j_mx", "p_dir", "ps_j")))) %dofuture% {
        ens_k <- ps_j_mx[k,,drop=F] %*% readRDS(glue("{p_dir}i_{as.integer(ps_j$i[k])}.rds"))
        c(mean(ens_k), sd(ens_k), quantile(ens_k, probs=c(0.005, 0.995)))
      }
      gc()
      
      ps_j <- ps_j |>
        mutate(ens_mn=ensIP[,1],
               ens_sd=ensIP[,2],
               ens_CL005=ensIP[,3],
               ens_CL995=ensIP[,4],
               ens_CI99width=ensIP[,4]-ensIP[,3]) |>
        select(i, starts_with("ens_"))
      gc()
      ps_j |>
        saveRDS(glue("{out_dir}/processed/hourly/Mature_{date_seq[i]}_{i_ts[j]}.rds"))
      ps_lims$ens_mn <- range(c(ps_lims$ens_mn, range(ps_j$ens_mn)))
      ps_lims$ens_sd <- range(c(ps_lims$ens_sd, range(ps_j$ens_sd)))
      ps_lims$ens_CL005 <- range(c(ps_lims$ens_CL005, range(ps_j$ens_CL005)))
      ps_lims$ens_CL995 <- range(c(ps_lims$ens_CL995, range(ps_j$ens_CL995)))
      ps_lims$ens_CI99width <- range(c(ps_lims$ens_CI99width, range(ps_j$ens_CI99width))) 
    }
    cat("", j)
  }
  cat("\n")
  plan(sequential)
  gc()
}
saveRDS(ps_lims, glue("{out_dir}/processed/hourly_Mature_pslims.rds"))





