# Project: Sealice IP Ensemble
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Prepare ensemble data

# This script creates the dataset for fitting and validating the ensemble models.
# It reads the non-interpolated reported fish and lice data, identifies starts
# of production cycles, and subsets records from the start of each cycle to the 
# first application of a lice treatment. The daily infection pressure is then 
# read in for each simulation, followed by lagging, reduction due to daily
# mortality rates, and summed to estimate the resulting relative number of adult
# lice predicted by each simulation for each day. This is aligned with the data
# from the farms and stored for use in the next scripts.


# setup -------------------------------------------------------------------
library(tidyverse); library(glue)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck"); # dirrf() and get_lags()


site_i <- read_csv("data/farm_sites.csv") |>
  left_join(read_csv("data/farm_sites_100m_areas.csv") |> select(sepaSite, area_m2) |> rename(area=area_m2))
init_df <- full_join(
  site_i,
  read_csv("data/lice_daily_2021-01-01_2024-12-31.csv", skip=1,
           col_names=c("sepaSite", paste0("d_", 0:1736)))
) |>
  pivot_longer(starts_with("d_"), names_to="date", values_to="density") |>
  mutate(date=ymd("2021-01-01") + as.numeric(str_sub(date, 3, -1))) |>
  left_join(read_csv("data/lice_biomass_2017-01-01_2024-12-31.csv") |>
              rename(fishTonnes=actualBiomassOnSiteTonnes))
out_dir <- "out/sim_2021-2024"
sim_i <- read_csv(glue("{out_dir}/sim_i.csv")) |>
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
                                   "2D", "3D", "Null")))


# load influx -------------------------------------------------------------

# read daily influx for all simulations
c_daily <- readRDS(glue("{out_dir}/processed/connectivity_day.rds")) |>
  mutate(sim=paste0("sim_", sim)) |>
  arrange(sim, sepaSite, date)



# identify production cycle starts to first lice treatments ---------------

newFarms_df <- init_df |>
  # for each farm, identify 'new' cycles by fishTonnes == 0 | fishTonnes < 25% of max
  group_by(sepaSite) |>
  mutate(fish_propMax=fishTonnes/max(fishTonnes, na.rm=T)) |>
  arrange(sepaSite, date) |>
  mutate(fishIn=!(is.na(fishTonnes) | fishTonnes==0 | fish_propMax < 0.2)) |>
  # find minimum of each period with 0-few fish
  group_by(sepaSite) |>
  mutate(fishInID=consecutive_id(fishIn)) |>
  filter(!fishIn) |>
  group_by(sepaSite, fishInID) |>
  slice_min(fishTonnes) |>
  slice_tail(n=1) |> # for consecutive 0s, select last day
  group_by(sepaSite) |>
  mutate(productionCycleNumber=row_number()) |>
  ungroup() |>
  select(sepaSite, date, productionCycleNumber) |>
  # join full dataset to fill in dates
  full_join(init_df, y=_) |>
  group_by(sepaSite) |>
  fill(productionCycleNumber) |>
  filter(!is.na(productionCycleNumber)) |>
  filter(fishTonnes > 0) |>
  select(sepaSite, date, productionCycleNumber, fishTonnes) |>
  # join with non-interpolated lice data
  inner_join(read_csv("data/lice_data_nonInterpolated.csv") |>
               rename(date=weekBeginning,
                      licePerFish=weeklyAverageAf) |>
               select(sepaSite, date, licePerFish, mitigation)) |>
  # identify first application of lice treatment (mitigation)
  group_by(sepaSite, productionCycleNumber) |>
  fill(mitigation, .direction="down") |>
  mutate(nDays=as.numeric(date - first(date))) |>
  ungroup() |>
  # take production cycle start until first treatment
  filter(is.na(mitigation))



# calculate relative expected adults from IP ------------------------------

maxLag <- 150 # presumed effective on-fish life span
defaults <- list(survival=c(0.970, 0.998, 0.997, 0.942),
                 devDays=c(15, 20, 20))
sensitivity_df <- expand_grid(survAdjust=c(0.9, 1, 1.1),
            devAdjust=c(0.9, 1, 1.1)) |>
  rowwise() |>
  mutate(survRates=list(plogis(qlogis(defaults$survival)*survAdjust)),
         devDays=list(c(defaults$devDays*devAdjust, maxLag-sum(defaults$devDays*devAdjust)))) |>
  ungroup()
valid_ls <- vector("list", nrow(sensitivity_df))
for(i in 1:nrow(sensitivity_df)) {
  surv <- c(1, rep(sensitivity_df$survRates[[i]], times=sensitivity_df$devDays[[i]]))
  surv_ls <- map(1:length(surv), ~prod(surv[1:.x])) |>
    setNames(paste0("IP", 0:(length(surv)-1)))
  valid_ls[[i]] <- c_daily |>
    select(sim, sepaSite, date, influx_m3) |>
    # fill in all dates for all sites (c_daily is sparse with 0's omitted)
    full_join(expand_grid(date=seq(ymd("2019-04-01"), ymd("2024-12-31"), by=1), 
                          sepaSite=unique(newFarms_df$sepaSite),
                          sim=unique(c_daily$sim))) |>
    mutate(influx_m3=replace_na(influx_m3, 0)) |>
    rename(IP=influx_m3) |>
    arrange(sim, sepaSite, date) |>
    # for each date, get lagged influx for up to 150 days 
    group_by(sim, sepaSite) |>
    get_lags(IP, n=maxLag) |>
    ungroup() |>
    # join only dates in newFarms_df
    inner_join(newFarms_df |> select(sepaSite, date, nDays)) |>
    rename(IP0=IP) |>
    # apply mortality for the corresponding number of days for each lag
    mutate(across(starts_with("IP"), ~.x*prod(surv_ls[[cur_column()]]))) |>
    # remove non-adult stages
    select(-any_of(paste0("IP", 0:(sum(sensitivity_df$devDays[[i]][1:2])-1)))) |>
    pivot_longer(starts_with("IP")) |>
    # remove IP from before the start of the production cycle
    filter(as.numeric(str_remove(name, "IP")) <= nDays) |>
    # sum (mortality reduced) IP for all days = relative number of predicted adults
    group_by(sim, sepaSite, date) |>
    summarise(IP_adult=sum(value)) |>
    ungroup() |>
    inner_join(newFarms_df) |>
    filter(between(date, min(c_daily$date), max(c_daily$date))) |>
    filter(date >= "2021-04-01") |>
    mutate(survAdjust=sensitivity_df$survAdjust[i],
           devAdjust=sensitivity_df$devAdjust[i])
}

valid_sens_df <- bind_rows(valid_ls)


# final processing --------------------------------------------------------

CV_dateSplits <- seq(min(valid_sens_df$date)-1, max(valid_sens_df$date), length.out=11)
valid_sens_df <- valid_sens_df |>
  # 4th root transformations due to right skew + 0's
  mutate(IP_adult_rtrt=IP_adult^0.25,
         licePerFish_rtrt=licePerFish^0.25) |>
  select(survAdjust, devAdjust, sim, sepaSite, productionCycleNumber, date, licePerFish_rtrt, IP_adult_rtrt) |>
  pivot_wider(names_from=sim, values_from=IP_adult_rtrt)|>
  # binary treatment by Code of Good Practice
  mutate(lice_g05=licePerFish_rtrt^4 > 0.5,
         year=year(date)) |>
  arrange(survAdjust, devAdjust, date, sepaSite) |>
  drop_na() |>
  group_by(survAdjust, devAdjust) |>
  mutate(rowNum=row_number(),
         sepaSiteNum=as.numeric(factor(sepaSite))) |>
  rowwise() |>
  mutate(sim_avgAll=mean(c_across(any_of(filter(sim_i, lab_short %in% c("2D", "3D"))$sim))),
         sim_avg3D=mean(c_across(any_of(filter(sim_i, lab_short=="3D")$sim))),
         sim_avg2D=mean(c_across(any_of(filter(sim_i, lab_short=="2D")$sim)))) |>
  mutate(CV_k=sum(date>CV_dateSplits)) |>
  ungroup()

write_csv(valid_sens_df, "out/valid_sens_df_2021-2024.csv")

valid_df <- valid_sens_df |> filter(survAdjust == 1 & devAdjust ==1)
write_csv(valid_df, "out/valid_df_2021-2024.csv")

