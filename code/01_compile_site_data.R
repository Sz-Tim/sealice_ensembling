# Project: Sealice IP Ensemble
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Compile fish farm data


# setup
library(tidyverse); library(glue)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
library(biotrackR) # devtools::install_github("Sz-Tim/biotrackR")
library(zoo)
library(janitor)
library(sf)
library(raster)
library(spaths)
library(readxl)
conflicted::conflicts_prefer(dplyr::select, 
                             dplyr::filter, 
                             tidyr::extract)
theme_set(theme_bw() + theme(panel.grid=element_blank()))
source("code/00_fn.R")



# load data ---------------------------------------------------------------

mesh.fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")

# Biomass is managed by SEPA, while sea lice is managed by MSS
# Separate site names are used, and there is not always total overlap
# SEPA sites are directly connected to licences, so they are used as the 
# locations, with the nearest sea lice counts applied to each

# Download data
download_lice_counts("2021-03-29", today(), "data/aquaculture_scot/ms_sea_lice_2021-present.csv")
download_fish_biomass("2017-01-01", today(), "data/aquaculture_scot/biomass_monthly_reports.csv")
download_ms_site_details("data/aquaculture_scot/ms_site_details.csv")
download_sepa_licenses("data/aquaculture_scot/se_licence_conditions.csv")

# Process datasets
sepa.locs <- read_csv("data/aquaculture_scot/se_licence_conditions.csv") |>
  filter(salmon==TRUE & waterType=="Seawater") |>
  st_as_sf(coords=c("easting", "northing"), crs=27700, remove=F) %>%
  mutate(dist=as.numeric(st_distance(., mesh.fp)[,1])) |>
  filter(dist < 1000) |>
  rename(sepaSite=sepaSiteId) |>
  mutate(operator=if_else(operator=="Wester Ross Salmon Ltd", "Wester Ross Fisheries Ltd", operator)) |>
  st_drop_geometry() |>
  # relocate sites that are outside WeStCOMS2 mesh due to coastal simplification
  mutate(easting=if_else(sepaSite=="RONAR1", 161050, easting),
         northing=if_else(sepaSite=="RONAR1", 855080, northing)) |>
  mutate(easting=if_else(sepaSite=="PLOI1", 206760, easting),
         northing=if_else(sepaSite=="PLOI1", 916230, northing)) |>
  st_as_sf(coords=c("easting", "northing"), crs=27700, remove=F) |>
  select(sepaSite, operator, easting, northing, geometry)
mss.locs <- read_csv("data/aquaculture_scot/ms_site_details.csv") |>
  filter(aquacultureType=="Fish" & waterType=="Seawater") |>
  st_as_sf(coords=c("easting", "northing"), crs=27700, remove=F) %>%
  mutate(dist=as.numeric(st_distance(., mesh.fp)[,1])) |>
  filter(dist < 1000) |>
  rename(mssID=marineScotlandSiteId) |>
  mutate(operator=if_else(operator=="Bakkafrost Scotland", "Bakkafrost Scotland Ltd", operator),
         operator=if_else(operator=="Landcatch Natural Selection Ltd", "Landcatch Ltd", operator)) |>
  select(mssID, operator, easting, northing, geometry)
sites.i <- sepa.locs %>%
  mutate(mssID=mss.locs$mssID[st_nearest_feature(., mss.locs)]) |>
  left_join(mss.locs |>
              st_drop_geometry() |>
              select(mssID, operator) |> 
              rename(mss_operator=operator))
sites.i$dist <- apply(st_distance(sites.i, mss.locs), 1, min)
sites_validation <- sites.i |>
  filter(operator == mss_operator) |>
  slice_min(dist, by=mssID) |>
  select(sepaSite, mssID, easting, northing, geometry)

biomass_df <- read_csv("data/aquaculture_scot/biomass_monthly_reports.csv") |>
  mutate(year=year(date)) |>
  filter(waterType=="Seawater",
         between(date, ymd("2017-01-01"), ymd("2024-12-31"))) |>
  select(-siteName, -easting, -northing) |>
  inner_join(sites.i |> st_drop_geometry(), by="sepaSite") |>
  arrange(sepaSite, date)

# The 'CLEAN' xlsx has replaced <, >, etc:
# >X<Y = mean(x,y)
# <X = X/2
# <CoGP = 0.1
# >X = X+1
# 998 = 0
lice_compiled <- clean_mss_lice_xlsx("data/aquaculture_scot/Sea+Lice+data+2017+-+2020.xlsx") |>
  bind_rows(read_xlsx("data/aquaculture_scot/Sea+Lice+Publication+until+week+12%2C++Mar+2021_CLEAN.xlsx") |>
              janitor::clean_names(case="small_camel") |>
              rename(weeklyAverageAf=avgAfNoFish,
                     mitigation=mitigations) |>
              select(starts_with("site"), 
                     weekBeginning, weeklyAverageAf, mitigation) |>
              mutate(weekBeginning=date(weekBeginning))) |>
  bind_rows(read_csv("data/aquaculture_scot/ms_sea_lice_2021-present.csv") |>
              select(starts_with("site"), 
                     weekBeginning, weeklyAverageAf, mitigation)) |>
  mutate(siteNo=str_to_upper(siteNo)) |>
  group_by(siteNo, weekBeginning) |>
  summarise(weeklyAverageAf=mean(weeklyAverageAf, na.rm=T),
            mitigation=paste(unique(mitigation), collapse=", ")) |>
  ungroup() |>
  mutate(mitigation=str_remove(mitigation, ", NA") |>
           str_remove(",$") |>
           str_replace(",,", ","),
         mitigation=na_if(mitigation, "NA"))
# Align with SEPA sites for lice release
lice_df <- lice_compiled |>
  filter(siteNo %in% sites.i$mssID) |>
  right_join(sites.i |> st_drop_geometry() |> select(sepaSite, mssID), y=_, 
             by=c("mssID"="siteNo"), relationship="many-to-many")
# Keep as-is for validation
lice_validation <- lice_compiled |> 
  rename(mssID=siteNo) |>
  filter(!is.na(weeklyAverageAf)) |>
  inner_join(sites_validation |> st_drop_geometry()) |>
  filter(weekBeginning > "2019-04-01") |>
  arrange(sepaSite, weekBeginning) |>
  select(sepaSite, mssID, weekBeginning, weeklyAverageAf, mitigation)
write_csv(lice_validation, "data/lice_data_nonInterpolated.csv")



# merge datasets ----------------------------------------------------------

# Reporting of lice counts was only required if averageAf exceeded a threshold:
#  - 2017.Jul - 2019.Jun.09: 3
#  - 2019.Jun.10 - 2021.Mar.28: 2
# https://www.gov.scot/publications/fish-health-inspectorate-sea-lice-information/
# 
# If no value was given, the smaller of median(day) and the median of (values 
# below the threshold from 2021-Apr to present) was used
default_df <- bind_cols(
  lice_df |>
    filter(weekBeginning > "2021-03-31" &
             weeklyAverageAf < 2) |>
    mutate(week=week(weekBeginning)) |> 
    group_by(week) |> 
    summarise(default_2=median(weeklyAverageAf, na.rm=T)),
  lice_df |>
    filter(weekBeginning > "2021-03-31" &
             weeklyAverageAf < 3) |>
    mutate(week=week(weekBeginning)) |> 
    group_by(week) |> 
    summarise(default_3=median(weeklyAverageAf, na.rm=T)) |>
    select(default_3),
  lice_df |>
    filter(weekBeginning > "2021-03-31") |>
    mutate(week=week(weekBeginning)) |> 
    group_by(week) |> 
    summarise(default_all=median(weeklyAverageAf, na.rm=T)) |>
    select(default_all)
)
default_df <- bind_rows(
  default_df,
  default_df |> filter(week %in% c(1, 52)) |>
    summarise(week=53,
              across(starts_with("default"), mean))
)

biomass_interp <- left_join(
  expand_grid(sepaSite=unique(biomass_df$sepaSite),
              date=seq(ymd("2017-01-01"), ymd("2024-12-31"), by=1)) |>
    mutate(year=year(date), month=month(date), day=day(date)),
  biomass_df |>
    st_drop_geometry() |>
    mutate(year=year(date),
           month=month(date)) |>
    select(sepaSite,
           year, month,
           actualBiomassOnSiteTonnes) |>
    complete(sepaSite, nesting(year, month), fill=list(actualBiomassOnSiteTonnes=0)) |>
    mutate(day=1)) |>
  # assume 0 biomass continues through the end of the month
  group_by(sepaSite, year, month) |>
  mutate(actualBiomassOnSiteTonnes=if_else(is.na(actualBiomassOnSiteTonnes) & 
                                             max(actualBiomassOnSiteTonnes, na.rm=T)==0, 
                                           0, 
                                           actualBiomassOnSiteTonnes)) |>
  # interpolate biomass to day:
  #   0 until first non-0 value, linear between values, constant to end
  group_by(sepaSite) |>
  mutate(actualBiomassOnSiteTonnes=as.numeric(na.fill(zoo(actualBiomassOnSiteTonnes), 
                                                      c(0, "extend", "extend")))) |>
  filter(any(actualBiomassOnSiteTonnes > 0)) |>
  ungroup()

lice_interp <- lice_df |>
  mutate(week=week(weekBeginning),
         year=year(weekBeginning)) |>
  right_join(expand_grid(sepaSite=unique(biomass_df$sepaSite),
                        year=2017:2024,
                        week=1:53)) |>
  mutate(weekBeginning=if_else(is.na(weekBeginning),
                               ymd(paste0(year, "-01-01"))+7*(week-1),
                               weekBeginning)) |>
  left_join(default_df, by=join_by(week)) |>
  mutate(weeklyAverageAf=case_when(
    !is.na(weeklyAverageAf) ~ weeklyAverageAf,
    is.na(weeklyAverageAf) & weekBeginning < ymd("2019-06-10") ~ default_3,
    is.na(weeklyAverageAf) & between(weekBeginning, ymd("2019-06-10"), ymd("2021-03-28")) ~ default_2,
    is.na(weeklyAverageAf) & weekBeginning > ymd("2021-03-28") ~ default_all
  )) |>
  select(sepaSite, weekBeginning, weeklyAverageAf) 

combo.df <- full_join(
  biomass_interp,
  lice_interp |>
    mutate(year=year(weekBeginning),
           month=month(weekBeginning),
           day=day(weekBeginning)) |>
    select(sepaSite, 
           year, month, day,
           weeklyAverageAf)
) |>
  group_by(sepaSite) |>
  mutate(weeklyAverageAf=as.numeric(na.fill(zoo(weeklyAverageAf), c(NA, "extend", NA)))) |>
  ungroup()
  

out.df <- combo.df |>
  mutate(weeklyAverageAf=pmax(weeklyAverageAf, 0.01)) |>
  mutate(total_AF=actualBiomassOnSiteTonnes * 240 * weeklyAverageAf) |>
  arrange(sepaSite, date) |> 
  group_by(sepaSite) |>
  mutate(N_obs=sum(actualBiomassOnSiteTonnes > 0, na.rm=T)) |>
  filter(N_obs > 0) |>
  ungroup() |>
  filter(date <= rollforward(max(biomass_df$date)))
  

# save output -------------------------------------------------------------

out.df |>
  select(sepaSite, date, weeklyAverageAf, actualBiomassOnSiteTonnes, total_AF) |>
  write_csv(glue("data/lice_biomass_{min(out.df$date)}_{max(out.df$date)}.csv"))
out.df |>
  filter(date >= "2021-01-01") |>
  group_by(sepaSite) |>
  summarise() |>
  left_join(sites.i |> st_drop_geometry() |> select(sepaSite, easting, northing)) |>
  write_csv("data/farm_sites.csv")
out.df |>
  filter(date >= "2021-01-01") |>
  mutate(date.c=str_remove_all(date, "-")) |>
  select(sepaSite, date.c, total_AF) |>
  pivot_wider(names_from=date.c, values_from=total_AF) |>
  filter(sepaSite %in% read_csv("data/farm_sites.csv")$sepaSite) |>
  mutate(across(where(is.numeric), ~replace_na(.x, 0))) |>
  write_csv(glue("data/lice_daily_2021-01-01_{max(out.df$date)}.csv"))

out.df |>
  filter(year(date) == 2023) |>
  group_by(sepaSite) |>
  summarise(active=any(actualBiomassOnSiteTonnes>0)) |>
  ungroup() |>
  filter(active) |>
  select(-active) |>
  left_join(sites.i |> st_drop_geometry() |> select(sepaSite, easting, northing)) |>
  write_csv(glue("data/farm_sites_2023.csv"))
out.df |>
  filter(year(date) == 2023) |>
  mutate(date.c=str_remove_all(date, "-")) |>
  select(sepaSite, date.c, total_AF) |>
  pivot_wider(names_from=date.c, values_from=total_AF) |>
  filter(sepaSite %in% read_csv("data/farm_sites_2023.csv")$sepaSite) |>
  mutate(across(where(is.numeric), ~replace_na(.x, 0))) |>
  write_csv(glue("data/lice_daily_2023.csv"))

out.df |>
  filter(between(date, ymd("2023-03-17"), ymd("2023-05-31"))) |>
  group_by(sepaSite) |>
  summarise(active=any(actualBiomassOnSiteTonnes>0)) |>
  ungroup() |>
  filter(active) |>
  select(-active) |>
  left_join(sites.i |> st_drop_geometry() |> select(sepaSite, easting, northing)) |>
  write_csv(glue("data/farm_sites_2023-MarMay.csv"))
out.df |>
  filter(between(date, ymd("2023-03-17"), ymd("2023-05-31"))) |>
  mutate(date.c=str_remove_all(date, "-")) |>
  select(sepaSite, date.c, total_AF) |>
  pivot_wider(names_from=date.c, values_from=total_AF) |>
  filter(sepaSite %in% read_csv("data/farm_sites_2023-MarMay.csv")$sepaSite) |>
  mutate(across(where(is.numeric), ~replace_na(.x, 0))) |>
  write_csv(glue("data/lice_daily_2023-MarMay.csv"))

out.df |>
  filter(between(date, ymd("2024-03-17"), ymd("2024-05-31"))) |>
  group_by(sepaSite) |>
  summarise(active=any(actualBiomassOnSiteTonnes>0)) |>
  ungroup() |>
  filter(active) |>
  select(-active) |>
  left_join(sites.i |> st_drop_geometry() |> select(sepaSite, easting, northing)) |>
  write_csv(glue("data/farm_sites_2024-MarMay.csv"))
out.df |>
  filter(between(date, ymd("2024-03-17"), ymd("2024-05-31"))) |>
  mutate(date.c=str_remove_all(date, "-")) |>
  select(sepaSite, date.c, total_AF) |>
  pivot_wider(names_from=date.c, values_from=total_AF) |>
  filter(sepaSite %in% read_csv("data/farm_sites_2024-MarMay.csv")$sepaSite) |>
  mutate(across(where(is.numeric), ~replace_na(.x, 0))) |>
  write_csv(glue("data/lice_daily_2024-MarMay.csv"))



# farm catchment footprints and volumes
mesh_sf <- st_read("data/WeStCOMS2_mesh.gpkg")
connect_dist <- c(100, 250, 500) # radius in m
for(i in connect_dist) {
  mesh_sf |>
    st_intersection(read_csv("data/farm_sites.csv") |>
                      st_as_sf(coords=c("easting", "northing"), crs=27700) |>
                      st_buffer(dist=i)) %>%
    mutate(area=st_area(.)) |>
    st_drop_geometry() |>
    mutate(vol=area*depth,
           vol_30m=area*pmin(depth, 30)) |>
    group_by(sepaSite) |>
    summarise(area_m2=as.numeric(sum(area)),
              vol_m3=as.numeric(sum(vol)),
              vol30m_m3=as.numeric(sum(vol_30m))) |>
    left_join(read_csv("data/farm_sites.csv")) |>
    write_csv(glue("data/farm_sites_{i}m_areas.csv"))
}

