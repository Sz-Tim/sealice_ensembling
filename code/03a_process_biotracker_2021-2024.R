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
  read_csv("data/lice_daily_2021-01-01_2024-12-31.csv", skip=1,
             col_names=c("sepaSite", paste0("d_", 0:1736)))
) |>
  pivot_longer(starts_with("d_"), names_to="date", values_to="density") |>
  mutate(date=ymd("2021-01-01") + as.numeric(str_sub(date, 3, -1))) |>
  left_join(read_csv("data/lice_biomass_2017-01-01_2024-12-31.csv") |>
              rename(fishTonnes=actualBiomassOnSiteTonnes))



# farm conditions ---------------------------------------------------------

env_df <- dirf(glue("{out_dir}/sim_01/"), "siteConditions") |>
  map_dfr(~read_csv(.x, show_col_types=F) |>
            select(site, depth, u, v, w, uv, salinity, temperature) |>
            mutate(date=ymd(str_split_fixed(basename(.x), "_", 3)[,2])))
saveRDS(env_df, glue("{out_dir}/processed/env_df.rds"))

env_df <- readRDS(glue("{out_dir}/processed/env_df.rds"))
p <- env_df |> 
  select(-depth) |> 
  pivot_longer(2:7) |> 
  mutate(name=case_when(name=="salinity" ~ "Salinity (psu)",
                        name=="temperature" ~ "Temperature (C)",
                        name=="u" ~ "Eastward water velocity (m/s)",
                        name=="v" ~ "Northward water velocity (m/s)",
                        name=="w" ~ "Upward water velocity (m/s)",
                        name=="uv" ~ "Current magnitude (m/s)")) |>
  mutate(name=factor(name)) |>
  mutate(name=lvls_reorder(name, c(5, 4, 2, 3, 6, 1))) |>
  ggplot(aes(date, value, group=site)) + 
  geom_line(alpha=0.05, colour="#084594") + 
  facet_wrap(~name, scales="free_y", ncol=2) +
  theme(axis.title=element_blank(),
        panel.grid.major.x=element_line(linewidth=0.2, colour="grey90"),
        panel.grid.major.y=element_line(linewidth=0.2, colour="grey90"))
ggsave("figs/pub_new/env_vars.png", p, width=9, height=10)




# particle densities ------------------------------------------------------

ps_wide_rtrt <- load_psteps_simSets(out_dir, mesh_i, sim_i, ncores=5, 
                                    liceScale=1/168, trans="4th_rt")
saveRDS(ps_wide_rtrt, glue("{out_dir}/processed/psteps_wide_rtrt.rds"))


# Averages: 4th root scale
ps_wide_rtrt <- readRDS(glue("{out_dir}/processed/psteps_wide_rtrt.rds"))

ps_long_rtrt <- ps_wide_rtrt |>
  pivot_longer(starts_with("t_"), names_to="date", values_to="rtrt_N_m2") |>
  mutate(date=first(init_df$date) + dhours(as.numeric(str_sub(date, 3, -1))))

tictoc::tic()
ps_avg_rtrt <- ps_long_rtrt |>
  calc_psteps_avg("rtrt_N_m2", ncores=5, mesh_sf=mesh_sf)
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
    saveRDS(glue("{out_dir}/processed/weekly/Mature_{ps_ts[i]}.rds"))
}





# connectivity ------------------------------------------------------------

# Mean hourly connectivity (day total / 24)
plan(multisession, workers=20)
c_long <- future_map_dfr(dirrf(out_dir, "connectivity.*csv"),
                         ~load_connectivity(.x, site_i$sepaSite, liceScale=1/24) |>
                           mutate(sim=str_sub(str_split_fixed(.x, "sim_", 3)[,3], 1, 2)))
plan(sequential)

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
  left_join(init_df |> select(sepaSite, date, weeklyAverageAf, fishTonnes), 
            by=join_by(sepaSite, date)) |>
  complete(sepaSite, sim, date, 
           fill=list(influx=0, influx_m2=0, N_influx=0,
                     outflux=0, outflux_m2=0, N_outflux=0,
                     self=0, self_m2=0))

saveRDS(c_daily, glue("{out_dir}/processed/connectivity_day.rds"))

