# Project: Sealice IP Ensemble
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Compile output from long runs

# These simulations are intended for visualizations of 2023 vertical distributions.



# setup -------------------------------------------------------------------

library(tidyverse) 
library(glue)
library(doFuture)
library(furrr)
library(carrier)
library(sf)
library(rstan)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
library(biotrackR) # devtools::install_github("Sz-Tim/biotrackR")
library(progressr)
theme_set(theme_bw() + theme(panel.grid=element_blank()))
handlers(global = TRUE)
options(future.globals.maxSize=5e9)

mesh_i <- st_read("data/WeStCOMS2_mesh.gpkg") |>
  mutate(vol_top50m=area*pmin(depth, 50)) |> 
  st_drop_geometry() |> 
  select(i, area, vol_top50m)
out_dir <- "D:/sealice_ensembling/out/sim_2023/"
sim_i <- read_csv(glue("{out_dir}/sim_i.csv")) |>
  mutate(sim=paste0("sim_", i)) 
site_i <- read_csv("data/farm_sites_2023.csv")
dat_ensMixRE <- readRDS(glue("out/ensembles/ensMix_all_sLonLatD4_FULL_standata.rds"))
p_dir <- "E:/sealice_ensembling/out/ensembles/p_meshCentroids/"

ncores <- 20



# calculate ensemble distributions ----------------------------------------

dates <- seq(ymd("2023-01-01"), ymd("2023-12-31"), by=1)
for(i in seq_along(dates)) {
  
  cat("Starting", i, "of", length(dates), "\n")
  if(file.exists(glue("{out_dir}/processed/verticalDistributions_{dates[i]}.rds"))) next
  
  z_df <- suppressMessages(
    load_vertDistr_simSets(out_dir, mesh_i |> select(i),
                           sim_i,
                           liceScale=1,
                           ncores=ncores,
                           stage="Mature",
                           date_rng=rep(dates[i], 2))) |>
    pivot_wider(names_from="sim", values_from="value") |>
    arrange(i, z)
  
  if(nrow(z_df) == 0) next
  
  plan(multisession, workers=ncores)
  z_mx <- z_df |> select(all_of(dat_ensMixRE$sim_names)) |> as.matrix()
  z_mx[is.na(z_mx)] <- 0
  z_df_i <- z_df |> select(i) |> mutate(ii=as.numeric(factor(i, levels=unique(i))))
  
  cat("  Ensembling", i, "of", length(dates), "-", nrow(z_df), "rows,", max(z_df_i$ii), "elements\n")
  ensIP <- ens_parallel(z_df_i, z_mx, p_dir)
  plan(sequential)
  gc()
  
  z_df <- z_df |>
    select(i, z, hour) |>
    mutate(ens_mn=ensIP[,1],
           ens_sd=ensIP[,2])
  saveRDS(z_df, glue("{out_dir}/processed/verticalDistributions_{dates[i]}.rds"))
  gc()
}





# summarise by region -----------------------------------------------------

mesh_sf <- st_read("data/WeStCOMS2_mesh.gpkg") |> select(i, geom)
linnhe_mesh <- mesh_sf |> 
  st_crop(c(xmin=150000, xmax=220000, ymin=710000, ymax=785000))
skye_mesh <- mesh_sf |> 
  st_crop(c(xmin=100000, xmax=198000, ymin=780000, ymax=920000))

f <- dirf(glue("{out_dir}/processed"), "verticalDistributions_")
z_ls <- z_linnhe <- z_skye <- vector("list", length(f))
for(d in 1:length(f)) {
  z_df <- readRDS(f[d])
  gc()

  z_linnhe[[d]] <- z_df |>
    filter(i %in% linnhe_mesh$i) |>
    mutate(time=ymd_hms("2023-01-01 00:00:00") + dhours(as.numeric(hour) - 1),
           day=date(time)) |>
    group_by(day, z) |>
    summarise(N=sum(ens_mn)/24,
              mean_ens_sd=mean(ens_sd),
              median_ens_sd=median(ens_sd)) |>
    group_by(day) |>
    mutate(prop=N/sum(N)) |>
    ungroup()
  gc()

  z_skye[[d]] <- z_df |>
    filter(i %in% skye_mesh$i) |>
    mutate(time=ymd_hms("2023-01-01 00:00:00") + dhours(as.numeric(hour) - 1),
           day=date(time)) |>
    group_by(day, z) |>
    summarise(N=sum(ens_mn)/24,
              mean_ens_sd=mean(ens_sd),
              median_ens_sd=median(ens_sd)) |>
    group_by(day) |>
    mutate(prop=N/sum(N)) |>
    ungroup()
  gc()

  z_ls[[d]] <- z_df |>
    mutate(time=ymd_hms("2023-01-01 00:00:00") + dhours(as.numeric(hour) - 1),
           day=date(time)) |>
    group_by(day, z) |>
    summarise(N=sum(ens_mn)/24,
              mean_ens_sd=mean(ens_sd),
              median_ens_sd=median(ens_sd)) |>
    group_by(day) |>
    mutate(prop=N/sum(N)) |>
    ungroup()
  rm(z_df)
  gc()
}
z_linnhe |>
  reduce(bind_rows) |>
  saveRDS(glue("{out_dir}/processed/summary_daily_z_Linnhe.rds"))
z_skye |>
  reduce(bind_rows) |>
  saveRDS(glue("{out_dir}/processed/summary_daily_z_Skye.rds"))
z_ls |>
  reduce(bind_rows) |>
  saveRDS(glue("{out_dir}/processed/summary_daily_z_WeStCOMS.rds"))



z_df <- map_dfr(c("WeStCOMS", "Linnhe", "Skye"),
                ~readRDS(glue("{out_dir}/processed/summary_daily_z_{.x}.rds")) |>
                  mutate(region=.x))
