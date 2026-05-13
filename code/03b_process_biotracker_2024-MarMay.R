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
# out_dir <- "out/sim_2021-2024"
# p_dir <- "out/ensembles/p_meshCentroids/"
p_dir <- "E:/sealice_ensembling/out/ensembles/p_meshCentroids/"
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

mod <- "n20_D4"
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
ensBlend_rec <- make_spline_recipe(meshCentroid_df, 
                                   as.numeric(str_split_fixed(mod, "D", 2)[,2]), 
                                   sim_i$sim[1:20])
# plan(multisession, workers=18)
plan(multicore, workers=18)
foreach(i=1:nrow(meshCentroid_df),
        .options.future=list(globals=structure(TRUE, add="p_dir"))) %dofuture% {
# for(i in 1:nrow(meshCentroid_df)) {
    t(make_predictions_ensBlend_sLonLat_RE(out_ensBlend,
                                           newdata=bake(ensBlend_rec, meshCentroid_df[i,]),
                                           iter=2000, seed=1003, 
                                           mode="b_p", type="surface")[1,,]) |>
    saveRDS(glue("{p_dir}/i_{meshCentroid_df$sepaSiteNum[i]}.rds"))
}
plan(sequential)





# vertical distributions at farms -----------------------------------------

mod <- "n20_D4"
dat_ensBlend <- readRDS(glue("out/ensembles/ensBlend_{mod}_FULL_standata.rds"))
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
                           # date_rng=c(dates[i], dates[i]+1))
                           date_rng=rep(dates[i], 2))
    ) |>
    filter(i %in% unique_elems) |>
    pivot_wider(names_from="sim", values_from="value") |>
    arrange(i, z)
  
  if(nrow(z_df) == 0) next
  
  plan(multisession, workers=18)
  z_mx <- z_df |> select(all_of(dat_ensBlend$sim_names)) |> as.matrix()
  z_mx[is.na(z_mx)] <- 0
  z_df_i <- z_df |> select(i) |> mutate(ii=as.numeric(factor(i, levels=unique(i))))

  cat("  Ensembling", i, "of", length(dates), "-", nrow(z_df), "rows,", max(z_df_i$ii), "elements\n")
  ensIP <- ens_parallel(z_df_i, z_mx, p_dir)
  plan(sequential)
  gc()

  z_df2 <- z_df |>
    mutate(sim_ens_mn=ensIP[,1],
           sim_ens_sd=ensIP[,2],
           sim_avg=rowMeans(z_mx)) |>
    left_join(farm_mesh_i, y=_, by=join_by(i), relationship="many-to-many") |> 
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
  saveRDS(z_df2, glue("{out_dir}/processed/verticalDistributions_{dates[i]}.rds"))
  gc()
  
}



# vertDist_df <- left_join(farm_mesh_i, z_df, by=join_by(i), relationship="many-to-many") |> 
#   filter(!is.na(hour)) |>
#   group_by(sepaSite, z, hour) |> 
#   summarise(across(starts_with("sim"), ~sum(.x, na.rm=T))) |>
#   ungroup() |>
#   complete(sepaSite, z, hour) |>
#   mutate(across(contains("sim"), ~if_else(is.na(.x), 0, .x)),
#          time=ymd("2024-03-17") + dhours(as.numeric(hour))) |>
#   left_join(farm_areas) |>
#   mutate(across(starts_with("sim"), ~.x/area)) |>
#   filter(z <= depth)

most_IP <- z_df2 |>
  group_by(sepaSite) |>
  pivot_longer(starts_with("sim_")) |>
  summarise(value=mean(value)) |>
  slice_max(value, n=4, with_ties=F)
example_sites <- c("MCLN1", "FFMC32", "TAB1", "INV1")
# 
# vertDist_df |>
#   filter(sepaSite %in% most_IP$sepaSite) |>
#   pivot_longer(starts_with("sim")) |>
#   group_by(sepaSite, z, time) |>
#   summarise(mn=mean(value)) |>
#   ungroup() |> 
#   ggplot(aes(time, z, fill=mn)) + 
#   geom_raster() + 
#   scale_y_reverse() +
#   scale_fill_viridis_c("Copepodids/m3/h", 
#                        option="turbo", end=0.95) +
#   facet_wrap(~sepaSite, nrow=4)
# 
z_df2 |>
  # filter(sepaSite %in% most_IP$sepaSite) |>
  filter(sepaSite %in% example_sites) |>
  ungroup() |>
  ggplot(aes(time, z, fill=sim_ens_mn)) +
  geom_raster() +
  scale_y_reverse() +
  scale_fill_viridis_c("Copepodids/m3/h",
                       option="turbo", end=0.95) +
  facet_wrap(~sepaSite, nrow=5)

z_df2 |>
  # filter(sepaSite %in% most_IP$sepaSite) |>
  filter(sepaSite %in% example_sites) |>
  ungroup() |>
  select(-sim_ens_sd) |>
  pivot_longer(starts_with("sim")) |>
  ggplot(aes(time, z, fill=value)) +
  geom_raster() +
  scale_y_reverse() +
  scale_fill_viridis_c("Copepodids/m3/h",
                       option="turbo", end=0.95) +
  facet_grid(sepaSite~name)
# 
# z_df2 |>
#   # filter(sepaSite %in% most_IP$sepaSite) |>
#   # filter(sepaSite %in% example_sites) |>
#   group_by(time, sepaSite) |>
#   mutate(sim_ens_p=sim_ens_mn/sum(sim_ens_mn, na.rm=T)) |>
#   ungroup() |>
#   ggplot(aes(time, z, fill=sim_ens_p)) +
#   geom_raster() +
#   scale_y_reverse() +
#   scale_fill_viridis_c("Copepodid distr.",
#                        option="turbo", end=0.95) +
#   facet_wrap(~sepaSite, nrow=5)
# 
# z_df2 |>
#   filter(z < 30) |>
#   filter(sepaSite %in% most_IP$sepaSite) |>
#   # filter(sepaSite %in% example_sites) |>
#   select(-matches("sim_[0-2]"), -"sim_ens_sd") |>
#   pivot_longer(starts_with("sim")) |>
#   ggplot(aes(time, z, fill=value)) +
#   geom_raster() +
#   scale_y_reverse() +
#   scale_fill_viridis_c("Copepodids/m3/h",
#                        option="turbo", end=0.95) +
#   facet_grid(name~sepaSite)
# 
# vertDist_df |>
#   filter(sepaSite %in% most_IP$sepaSite) |>
#   pivot_longer(starts_with("sim")) |>
#   group_by(sepaSite, z, time) |>
#   summarise(mn=mean(value)) |>
#   group_by(sepaSite, time) |>
#   mutate(prop=mn/sum(mn)) |>
#   ungroup() |> 
#   ggplot(aes(time, z, fill=prop)) + 
#   geom_raster() + 
#   scale_y_reverse() +
#   scale_fill_viridis_c("Proportion\nof copepodids", na.value=NA,
#                        option="turbo", limits=c(0, 1), end=0.95) +
#   facet_wrap(~sepaSite, nrow=4)
# 
# vertDist_df |>
#   # filter(sepaSite %in% most_IP$sepaSite) |>
#   filter(sepaSite %in% example_sites) |>
#   pivot_longer(starts_with("sim")) |>
#   group_by(sepaSite, name, time) |>
#   mutate(prop=value/sum(value)) |>
#   ungroup() |> 
#   ggplot(aes(time, z, fill=prop)) + 
#   geom_raster() + 
#   scale_y_reverse() +
#   scale_fill_viridis_c("Proportion\nof copepodids", na.value=NA,
#                        option="turbo", limits=c(0, 1), end=0.95) +
#   facet_grid(name~sepaSite)
# 
# vertDist_df |>
#   # filter(sepaSite %in% most_IP$sepaSite) |>
#   filter(sepaSite == "MCLN1") |>
#   pivot_longer(starts_with("sim")) |>
#   group_by(sepaSite, name, time) |>
#   mutate(prop=value/sum(value)) |>
#   ungroup() |> 
#   ggplot(aes(time, z, fill=prop)) + 
#   geom_raster() + 
#   scale_y_reverse() +
#   scale_fill_viridis_c("Proportion\nof copepodids", na.value=NA,
#                        option="turbo", limits=c(0, 1), end=0.95) +
#   facet_wrap(~name, nrow=5)
