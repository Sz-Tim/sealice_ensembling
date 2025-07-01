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
library(progressr)
source("code/00_fn.R")
theme_set(theme_bw() + theme(panel.grid=element_blank()))
handlers(global = TRUE)
options(future.globals.maxSize=5e9)

out_dir <- "D:/sealice_ensembling/out/sim_2024-MarMay"
p_dir <- "out/ensembles/p_meshCentroids/"
sim_i <- read_csv(glue("{out_dir}/sim_i.csv")) |>
  mutate(sim=paste0("sim_", i))
site_i <- read_csv("data/farm_sites_2024-MarMay.csv")
mesh_i <- st_read("data/WeStCOMS2_mesh.gpkg") |>
  st_drop_geometry() |> 
  select(i, area, depth)
date_range <- c(ymd("2024-04-01"), ymd("2024-05-31"))




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
                              sim_i, ncores=18, liceScale=1, 
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



# vertical distributions at farms -----------------------------------------

farm_mesh_i <- read_csv("data/farm_sites_100m_meshElems.csv")
unique_elems <- sort(unique(farm_mesh_i$i))
farm_areas <- left_join(farm_mesh_i, mesh_i) |>
  group_by(sepaSite) |>
  summarise(area=sum(area),
            depth=max(depth))

dates <- seq(date_range[1], date_range[2], by=1)
for(i in seq_along(dates)) {
  
  cat("Starting", i, "of", length(dates), "\n")
  
  z_df <- suppressMessages(
    load_vertDistr_simSets(out_dir, mesh_i |> select(i),
                           sim_i,
                           liceScale=1,
                           ncores=18,
                           stage="Mature",
                           date_rng=c(dates[i], dates[i]+3))
                           # date_rng=rep(dates[i], 2))
    ) |>
    filter(i %in% unique_elems) |>
    pivot_wider(names_from="sim", values_from="value") |>
    arrange(i, z)
  
  if(nrow(z_df) == 0) next
  
  vertDist_df <- left_join(farm_mesh_i, z_df, by=join_by(i), relationship="many-to-many") |> 
    filter(!is.na(hour)) |>
    group_by(sepaSite, z, hour) |> 
    summarise(across(starts_with("sim"), ~sum(.x, na.rm=T))) |>
    ungroup() |>
    complete(sepaSite, z, hour) |>
    mutate(across(contains("sim"), ~if_else(is.na(.x), 0, .x)),
           time=ymd("2024-03-17") + dhours(as.numeric(hour))) |>
    left_join(farm_areas) |>
    mutate(across(starts_with("sim"), ~.x/area)) |>
    filter(z <= depth)
  
  # TODO: Need to rework for site-level distributions
  # plan(multisession, workers=ncores)
  # z_mx <- z_df |> select(all_of(dat_ensBlend$sim_names)) |> as.matrix()
  # z_mx[is.na(z_mx)] <- 0
  # z_df_i <- z_df |> select(i) |> mutate(ii=as.numeric(factor(i, levels=unique(i))))
  # 
  # cat("  Ensembling", i, "of", length(dates), "-", nrow(z_df), "rows,", max(z_df_i$ii), "elements\n")
  # ensIP <- ens_parallel(z_df_i, z_mx, p_dir)
  # plan(sequential)
  # gc()
  # 
  # z_df <- z_df |>
  #   select(i, z, hour) |>
  #   mutate(ens_mn=ensIP[,1],
  #          ens_sd=ensIP[,2])
  # saveRDS(z_df, glue("{out_dir}/processed/verticalDistributions_{dates[i]}.rds"))
  # gc()
  
  most_IP <- vertDist_df |> 
    group_by(sepaSite) |> 
    pivot_longer(starts_with("sim_")) |> 
    summarise(value=mean(value)) |> 
    slice_max(value, n=4, with_ties=F)
  example_sites <- c("MCLN1", "FFMC32", "TAB1", "INV1")
  
  vertDist_df |>
    filter(sepaSite %in% most_IP$sepaSite) |>
    pivot_longer(starts_with("sim")) |>
    group_by(sepaSite, z, time) |>
    summarise(mn=mean(value)) |>
    ungroup() |> 
    ggplot(aes(time, z, fill=mn)) + 
    geom_raster() + 
    scale_y_reverse() +
    scale_fill_viridis_c("Copepodids/m3/h", 
                         option="turbo", end=0.95) +
    facet_wrap(~sepaSite, nrow=4)
  
  vertDist_df |>
    # filter(sepaSite %in% most_IP$sepaSite) |>
    filter(sepaSite %in% example_sites) |>
    pivot_longer(starts_with("sim")) |>
    ggplot(aes(time, z, fill=value)) + 
    geom_raster() + 
    scale_y_reverse() +
    scale_fill_viridis_c("Copepodids/m3/h", 
                         option="turbo", end=0.95) +
    facet_grid(name~sepaSite)
  
  vertDist_df |>
    filter(sepaSite %in% most_IP$sepaSite) |>
    pivot_longer(starts_with("sim")) |>
    group_by(sepaSite, z, time) |>
    summarise(mn=mean(value)) |>
    group_by(sepaSite, time) |>
    mutate(prop=mn/sum(mn)) |>
    ungroup() |> 
    ggplot(aes(time, z, fill=prop)) + 
    geom_raster() + 
    scale_y_reverse() +
    scale_fill_viridis_c("Proportion\nof copepodids", na.value=NA,
                         option="turbo", limits=c(0, 1), end=0.95) +
    facet_wrap(~sepaSite, nrow=4)
  
  vertDist_df |>
    # filter(sepaSite %in% most_IP$sepaSite) |>
    filter(sepaSite %in% example_sites) |>
    pivot_longer(starts_with("sim")) |>
    group_by(sepaSite, name, time) |>
    mutate(prop=value/sum(value)) |>
    ungroup() |> 
    ggplot(aes(time, z, fill=prop)) + 
    geom_raster() + 
    scale_y_reverse() +
    scale_fill_viridis_c("Proportion\nof copepodids", na.value=NA,
                         option="turbo", limits=c(0, 1), end=0.95) +
    facet_grid(name~sepaSite)
  
  vertDist_df |>
    # filter(sepaSite %in% most_IP$sepaSite) |>
    filter(sepaSite == "MCLN1") |>
    pivot_longer(starts_with("sim")) |>
    group_by(sepaSite, name, time) |>
    mutate(prop=value/sum(value)) |>
    ungroup() |> 
    ggplot(aes(time, z, fill=prop)) + 
    geom_raster() + 
    scale_y_reverse() +
    scale_fill_viridis_c("Proportion\nof copepodids", na.value=NA,
                         option="turbo", limits=c(0, 1), end=0.95) +
    facet_wrap(~name, nrow=5)
  
}

