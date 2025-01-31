# Project: Sealice IP Ensemble
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Compile output from long runs




# setup -------------------------------------------------------------------

library(tidyverse); library(glue)
library(furrr)
library(carrier)
library(sf)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
library(biotrackR) # devtools::install_github("Sz-Tim/biotrackR")
theme_set(theme_bw() + theme(panel.grid=element_blank()))

mesh_sf <- st_read("data/WeStCOMS2_mesh.gpkg") |> mutate(vol_top50m=area*pmin(depth, 50))
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_i <- mesh_sf |> st_drop_geometry() |> select(i, area, vol_top50m)
out_dir <- "out/sim_2021-2024"
sim_i <- read_csv(glue("{out_dir}/sim_i.csv")) |>
  mutate(sim=paste0("sim_", i))
site_i <- read_csv("data/farm_sites.csv")
init_df <- full_join(
  site_i,
  read_csv("data/lice_daily_2021-01-01_2024-09-30.csv", skip=1,
             col_names=c("sepaSite", paste0("d_", 0:1736)))
) |>
  pivot_longer(starts_with("d_"), names_to="date", values_to="density") |>
  mutate(date=ymd("2021-01-01") + as.numeric(str_sub(date, 3, -1))) |>
  left_join(read_csv("data/lice_biomass_2017-01-01_2024-09-30.csv") |>
              rename(fishTonnes=actualBiomassOnSiteTonnes))




# extract output ----------------------------------------------------------

# extract biotracker output: sim_[0-9][0-9].tar.gz
if(FALSE) {
  f_tgz <- dirf(out_dir, "tar.gz")
  walk(f_tgz, ~untar(.x, exdir=str_remove(.x, ".tar.gz")))
}



# particle densities ------------------------------------------------------

ps_wide_rtrt <- load_psteps_simSets(out_dir, mesh_i, sim_i, ncores=30, liceScale=1/168, trans="4th_rt")
saveRDS(ps_wide_rtrt, glue("{out_dir}/processed/psteps_wide_rtrt.rds"))


# Averages: 4th root scale
ps_wide_rtrt <- readRDS(glue("{out_dir}/processed/psteps_wide_rtrt.rds"))

ps_long_rtrt <- ps_wide_rtrt |>
  pivot_longer(starts_with("t_"), names_to="date", values_to="rtrt_N_m2") |>
  mutate(date=first(init_df$date) + dhours(as.numeric(str_sub(date, 3, -1))))

tictoc::tic()
ps_avg_rtrt <- ps_long_rtrt |>
  calc_psteps_avg("rtrt_N_m2", ncores=30, mesh_sf=mesh_sf)
tictoc::toc()
saveRDS(ps_avg_rtrt, glue("{out_dir}/processed/psteps_avg_rtrt.rds"))

rm(list=ls() |> grep("ps_.*_rtrt", x=_, value=T))
gc()




ps_df <- readRDS(glue("{out_dir}/processed/psteps_wide_rtrt.rds"))
ps_ts <- grep("^t_", names(ps_df), value=T)
for(i in seq_along(ps_ts)) {
  ps_df |> 
    select(sim, i, all_of(ps_ts[i])) |>
    drop_na() |>
    arrange(sim, i) |>
    pivot_wider(names_from=sim, values_from=ps_ts[i]) |>
    mutate(across(starts_with("sim"), ~replace_na(.x, 0))) |>
    saveRDS(paste0("out/sim_2019-2023/processed/weekly/Mature_", ps_ts[i], ".rds"))
}





# connectivity ------------------------------------------------------------

# Mean hourly connectivity (day total / 24)
plan(multisession, workers=40)
c_long <- future_map_dfr(dirrf(out_dir, "connectivity.*csv"),
                         ~load_connectivity(.x, site_i$sepaSite, liceScale=1/24) |>
                           mutate(sim=str_sub(str_split_fixed(.x, "sim_", 3)[,3], 1, 2)))
plan(sequential)

# # filter so that only sites with fish are included
# active_df <- init_df |> filter(fishTonnes > 0) |> select(sepaSite, date)
# c_long <- c_long |>
#   inner_join(active_df, by=join_by(date, destination==sepaSite)) 

site_areas <- read_csv("data/farm_sites_100m_areas.csv") |> 
  select(sepaSite, area_m2) |> 
  rename(area=area_m2)

# mean hourly IP for each day
c_daily <- list(
  c_long |> 
    calc_influx(destination, value, sim, date) |> 
    rename(sepaSite=destination) |>
    left_join(site_areas) |>
    mutate(influx_m2=influx/area) |> 
    select(-area),
  c_long |> 
    calc_outflux(source, value, sim, date, 
                 dest_areas=site_areas |> rename(destination=sepaSite)) |> 
    rename(sepaSite=source),
  c_long |> 
    calc_self_infection(source, destination, value, sim, date) |> rename(sepaSite=source) |>
    left_join(site_areas) |>
    mutate(self_m2=self/area) |> 
    select(-area)
) |>
  reduce(full_join) |>
  left_join(init_df) |>
  # mutate(fishTonnes=if_else(is.na(fishTonnes), 0, fishTonnes)) |>
  # filter(fishTonnes > 0) |>
  complete(sepaSite, sim, date, 
           fill=list(influx=0, influx_m2=0, N_influx=0,
                     outflux=0, outflux_m2=0, N_outflux=0,
                     self=0, self_m2=0))

saveRDS(c_daily, glue("{out_dir}/processed/connectivity_day.rds"))




c_daily |> 
  ggplot(aes(date, influx_m2^0.25, colour=sim)) + 
  geom_line(linewidth=0.2, alpha=0.25) + 
  facet_wrap(~sepaSite) + 
  theme_classic()

c_daily |> 
  ggplot(aes(date, outflux_m2-self_m2, colour=sim)) + 
  geom_line() + 
  facet_wrap(~sepaSite) + 
  theme_classic()


c_daily |>
  filter(sepaSite=="CALL1")
