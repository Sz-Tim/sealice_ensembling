# Project: Sealice IP Ensemble
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Short simulations


# setup
library(tidyverse); library(glue)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
library(biotrackR) # devtools::install_github("Sz-Tim/biotrackR")
library(doFuture)
theme_set(theme_bw() + theme(panel.grid=element_blank()))
source("code/00_fn.R")



# define parameters -------------------------------------------------------

cores_per_sim <- 25
parallel_sims <- 1
start_date <- "2024-03-17"
end_date <- "2024-05-31"
nDays <- length(seq(ymd(start_date), ymd(end_date), by=1))

set.seed(1001)

os <- get_os()
dirs <- switch(
  get_os(),
  linux=list(proj=getwd(),
             mesh="/home/sa04ts/hydro/meshes",
             hydro="/home/sa04ts/hydro/WeStCOMS2/Archive",
             jdk="/home/sa04ts/.jdks/jdk-23.0.1/bin/java",
             jar="/home/sa04ts/biotracker/biotracker_v1-0-0.jar",
             out=glue("{getwd()}/out/sim_2024-MarMay")),
  windows=list(proj=getwd(),
               mesh="E:/hydro",
               hydro="E:/hydro/WeStCOMS2/Archive",
               jdk="C:/Users/sa04ts/.jdks/openjdk-23.0.2/bin/java",
               jar="C:/Users/sa04ts/OneDrive - SAMS/Projects/03_packages/biotracker/out/biotracker_v1-0-0.jar",
               out=glue("D:/sealice_ensembling/out/sim_2024-MarMay"))
)

n_sim3D <- 16
n_sim2D <- 4
light_mx <- MASS::mvrnorm(n_sim3D, c(0,0), matrix(c(1, 0.8, 0.8, 1), nrow=2))
swim_mx <- MASS::mvrnorm(n_sim3D, c(0,0), matrix(c(1, 0.8, 0.8, 1), nrow=2))
sim.i <- bind_rows(
  tibble(fixDepth="false",
         D_hVert=exp(runif(n_sim3D, log(1e-5), log(1e-1))),
         mortSal_fn=sample(c("constant", "logistic"), n_sim3D, replace=T),
         eggTemp_fn=sample(c("constant", "logistic"), n_sim3D, replace=T),
         salinityThreshMin=runif(n_sim3D, 20, 28),
         salinityThreshMax=pmin(salinityThreshMin + runif(n_sim3D, 0.1, 6), 32),
         lightThreshCopepodid=qunif(pnorm(light_mx[,1]), (2e-6)^0.5, (2e-4)^0.5)^2,
         lightThreshNauplius=qunif(pnorm(light_mx[,2]), (0.05)^0.5, (0.5)^0.5)^2,
         swimUpSpeedMean=-(qunif(pnorm(swim_mx[,1]), (1e-4)^0.5, (2e-2)^0.5))^2,
         swimDownSpeedMean=(qunif(pnorm(swim_mx[,2]), (1e-4)^0.5, (2e-2)^0.5))^2),
  expand_grid(mortSal_fn=c("constant", "logistic"),
              eggTemp_fn=c("constant", "logistic")) |>
    mutate(fixDepth="true",
           D_h=exp(runif(n_sim2D, log(1e-3), log(1e1))))
) |>
  mutate(across(where(is.numeric), ~if_else(is.na(.x), 0, .x))) |>
  mutate(D_h=exp(runif(n_sim2D+n_sim3D, log(1e-3), log(1e1))),
         viableDegreeDays=runif(n_sim2D+n_sim3D, 30, 50)) |>
  mutate(across(where(is.numeric), ~signif(.x, 5))) |>
  mutate(i=str_pad(row_number(), 2, "left", "0"),
         outDir=glue("{dirs$out}/sim_{i}/"))
write_csv(sim.i, glue("{dirs$out}/sim_i.csv")) 
sim_seq <- 1:nrow(sim.i)


# set properties ----------------------------------------------------------

walk(sim_seq, ~dir.create(sim.i$outDir[.x], showWarnings=F))

walk(sim_seq,
     ~set_biotracker_properties2(
       # run settings
       properties_file_path=glue("{dirs$out}/sim_{sim.i$i[.x]}.properties"),
       parallelThreads=cores_per_sim,
       parallelThreadsHD=6,
       start_ymd=as.numeric(str_remove_all(start_date, "-")),
       numberOfDays=nDays,
       nparts=100,
       checkOpenBoundaries="true",
       # meshes and environment
       mesh1=glue("{dirs$mesh}/WeStCOMS2_mesh.nc"),
       mesh1Domain="westcoms2",
       datadir=glue("{dirs$hydro}/"),
       # sites
       sitefile=glue("D:/sealice_ensembling/data/farm_sites_2024-MarMay.csv"),
       sitefileEnd=glue("D:/sealice_ensembling/data/farm_sites_2024-MarMay.csv"),
       siteDensityPath=glue("D:/sealice_ensembling/data/lice_daily_2024-MarMay.csv"),
       # dynamics
       D_h=sim.i$D_h[.x],
       D_hVert=sim.i$D_hVert[.x],
       stepsPerStep=30,
       # biology
       fixDepth=sim.i$fixDepth[.x],
       startDepth=1,
       eggTemp_fn=sim.i$eggTemp_fn[.x],
       mortSal_fn=sim.i$mortSal_fn[.x],
       salinityThreshMin=sim.i$salinityThreshMin[.x],
       salinityThreshMax=sim.i$salinityThreshMax[.x],
       lightThreshCopepodid=sim.i$lightThreshCopepodid[.x],
       lightThreshNauplius=sim.i$lightThreshNauplius[.x],
       swimUpSpeedCopepodidMean=sim.i$swimUpSpeedMean[.x],
       swimUpSpeedCopepodidStd=abs(sim.i$swimUpSpeedMean[.x]/5),
       swimDownSpeedCopepodidMean=sim.i$swimDownSpeedMean[.x],
       swimDownSpeedCopepodidStd=abs(sim.i$swimDownSpeedMean[.x]/5),
       swimUpSpeedNaupliusMean=sim.i$swimUpSpeedMean[.x]/2,
       swimUpSpeedNaupliusStd=abs(sim.i$swimUpSpeedMean[.x]/10),
       swimDownSpeedNaupliusMean=sim.i$swimDownSpeedMean[.x]/2,
       swimDownSpeedNaupliusStd=abs(sim.i$swimDownSpeedMean[.x]/10),
       passiveSinkingIntercept=sim.i$swimDownSpeedMean[.x],
       passiveSinkingSlope=0,
       viableDegreeDays=sim.i$viableDegreeDays[.x],
       connectivityThresh=100,
       # recording
       verboseSetUp="true",
       recordConnectivity="false",
       recordPsteps="true",
       splitPsteps="false",
       pstepsInterval=1,
       pstepsMaxDepth=5,
       recordVertDistr="true",
       vertDistrInterval=1,
       vertDistrMax=30,
       recordElemActivity="false"))


# run simulations ---------------------------------------------------------

if(os=="windows") {
  plan(multisession, workers=parallel_sims)
} else {
  plan(multicore, workers=parallel_sims)
}

sim_sets <- split(sim_seq, rep(1:parallel_sims, length(sim_seq)/parallel_sims))
foreach(j=1:parallel_sims, .options.future=list(globals=structure(TRUE, add="sim.i"))) %dofuture% {
  for(i in sim_sets[[j]]) {
    setwd(dirs$proj)
    biotrackR::run_biotracker(
      jdk_path=dirs$jdk,
      jar_path=dirs$jar,
      f_properties=glue::glue("{dirs$out}/sim_{sim.i$i[i]}.properties"),
      sim_dir=glue::glue("{sim.i$outDir[i]}")
    )
  }
}
plan(sequential)
