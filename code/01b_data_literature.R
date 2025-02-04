# Project: Sealice IP Ensemble
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Estimate functions from literature


# This script fits functions for egg production ~ Temperature and larval 
# mortality ~ Salinity using data from the literature. Logistic models were
# selected as most suitable for each. Models were fit using Bayesian regression
# and the posterior mean was used for each parameter.


# setup -------------------------------------------------------------------
library(tidyverse)
library(brms)
theme_set(theme_bw())



# egg production ----------------------------------------------------------

# Scottish data + biotracker assume eggs/day/AF, NOT eggs/day/Gravid.
# egg_df is adjusted to assume 82.6% of adult females are gravid based on
# average development times in Toorians & Adams 2020 (35d: AF, 55-150d: Gravid)

# With data aggregated in Kragesteen 2023 (Table 1), Bayesian regression was
# used to fit 1) linear, 2) quadratic, and 3) logistic functions, where the 
# form of the quadratic function follows Stein 2005 and current Norwegian models

# . data ------------------------------------------------------------------
# Kragesteen 2023: Egg = f(temperature)
# typo in Table 1: Hamre 2019 has 28.9, not 28.6 for T=6
egg_df <- read_csv("data/lit/Kragesteen2023_Table1.csv") |>
  filter(temperature > 4) |>
  mutate(eggs_day_AG=if_else(temperature==6, 28.9, eggs_day_AG),
         eggs_day_AF=eggs_day_AG*(150-55)/(150-35))

# . linear ----------------------------------------------------------------
egg_temp_linear <- brm(
  eggs_day_AF ~ temperature, data=egg_df, cores=4,
  prior=c(prior(normal(0, 10), class=sigma)))
egg_temp_linear |> saveRDS("data/lit/fit/egg_temp_linear_brmsfit.rds")
as_draws_df(egg_temp_linear) |>
  rename(a=b_Intercept,
         b=b_temperature) |>
  select(a, b) |>
  write_csv("data/lit/fit/egg_temp_linear_post.csv")

# . quadratic -------------------------------------------------------------
egg_temp_quadratic <- brm(
  bf(eggs_day_AF ~ a * (b + temperature)^2, nl=T) +
    lf(a ~ 1, b ~ 1),
  prior=c(prior(normal(0.17, 3), nlpar="a"),
          prior(normal(4.28, 3), nlpar="b"),
          prior(normal(0, 3), class=sigma)),
  data=egg_df, cores=4)
egg_temp_quadratic |> saveRDS("data/lit/fit/egg_temp_quadratic_brmsfit.rds")
as_draws_df(egg_temp_quadratic) |>
  rename(a=b_a_Intercept, b=b_b_Intercept) |>
  select(a, b) |>
  write_csv("data/lit/fit/egg_temp_quadratic_post.csv")

# . logistic --------------------------------------------------------------
egg_temp_logistic <- brm(
  bf(eggs_day_AF ~ eggMax/(1+exp(-k*(temperature - tempMid))) + eggMin, nl=T) +
    lf(eggMax ~ 1, k~1, tempMid~1, eggMin~1),
  data=egg_df, cores=4, control=list(adapt_delta=0.99),
  prior=c(prior(normal(70, 5), nlpar="eggMax"),
          prior(normal(0, 1), nlpar="k"),
          prior(normal(10, 3), nlpar="tempMid"),
          prior(normal(10, 5), nlpar="eggMin", lb=0),
          prior(normal(0, 5), class=sigma)))
egg_temp_logistic |> saveRDS("data/lit/fit/egg_temp_logistic_brmsfit.rds")
as_draws_df(egg_temp_logistic) |>
  rename(eggMax=b_eggMax_Intercept,
         k=b_k_Intercept,
         tempMid=b_tempMid_Intercept,
         eggMin=b_eggMin_Intercept) |>
  select(eggMax, k, tempMid, eggMin) |>
  write_csv("data/lit/fit/egg_temp_logistic_post.csv")

# . visualization ---------------------------------------------------------
post_linear <- read_csv("data/lit/fit/egg_temp_linear_post.csv")
post_quadratic <- read_csv("data/lit/fit/egg_temp_quadratic_post.csv")
post_logistic <- read_csv("data/lit/fit/egg_temp_logistic_post.csv")

temp_seq <- seq(3, 20, by=0.1)
plot(NA, NA, xlim=c(3, 20), ylim=c(0, 100),
     xlab="Temperature (C)", ylab="Mean eggs per adult female per day")
for(i in 1:100) {
  lines(temp_seq, 
        post_linear$a[i] + post_linear$b[i] * temp_seq, 
        col=rgb(217, 95, 2, 25, maxColorValue=256))
  lines(temp_seq, 
        post_quadratic$a[i] * (post_quadratic$b[i] + temp_seq)^2, 
        col=rgb(117, 112, 179, 25, maxColorValue=256))
  lines(temp_seq, 
        post_logistic$eggMax[i] / 
          (1+exp(-post_logistic$k[i]*(temp_seq - post_logistic$tempMid[i]))) + 
          post_logistic$eggMin[i], 
        col=rgb(27, 158, 119, 25, maxColorValue=256))
}
points(egg_df$temperature, egg_df$eggs_day_AF, col="black", pch=19)
abline(h=28.2, lty=3)
legend("topleft", 
       c("Linear", "Quadratic", "Logistic", "Constant (28.2)"), 
       col=c(rgb(217,95,2, maxColorValue=256), 
             rgb(117, 112, 179, maxColorValue=256), 
             rgb(27, 158, 119, maxColorValue=256),
             "black"),
       lwd=1, 
       lty=c(1,1,1,3))

# . posterior mean --------------------------------------------------------
post_logistic |> 
  summarise(across(everything(), mean))



# mortality ---------------------------------------------------------------
# With data from Bricknell 2006 (Table 1), Bayesian regression was
# used to fit a logistic function

# . data ------------------------------------------------------------------
# Assume survival times follow an exponential distribution
# hazard rate = hourly mortality rate = ln(2)/LT50

# sal = 5, 9, 12 are listed as LT50: <1h
# "At 12 ppt and below the initial death rate was rapid, with all copepodids in
# 9 and 5 ppt dying within the first 2 h."
# Assume this translates to 5ppt = 0.95 and 9ppt = 0.9 
sal_df <- read_csv("data/lit/Bricknell2006_Table1.csv") |>
  mutate(mort_h=-log(0.5)/LT50_h) |>
  mutate(mort_h=case_when(salinity==5 ~ 0.95,
                          salinity==9 ~ 0.9,
                          .default=mort_h))

# . logistic --------------------------------------------------------------
mort_sal_logistic <- brm(
  bf(mort_h ~ mortMax/(1+exp(-k*(salinity - salMid))) + mortMin, nl=T) +
    lf(mortMax ~ 1, k~1, salMid~1, mortMin~1),
  data=sal_df, cores=4,
  prior=c(prior(normal(1, 0.1), nlpar="mortMax", lb=0, ub=1),
          prior(normal(0, 1), nlpar="k"),
          prior(normal(15, 4), nlpar="salMid"),
          prior(normal(0.01, 0.01), nlpar="mortMin", lb=0, ub=1),
          prior(normal(0, 0.4), class=sigma)))
mort_sal_logistic |> saveRDS("data/lit/fit/mort_sal_logistic_brmsfit.rds")
as_draws_df(mort_sal_logistic) |>
  rename(mortMax=b_mortMax_Intercept,
         k=b_k_Intercept,
         salMid=b_salMid_Intercept,
         mortMin=b_mortMin_Intercept) |>
  select(mortMax, k, salMid, mortMin) |>
  write_csv("data/lit/fit/mort_sal_logistic_post.csv")


# . visualization ---------------------------------------------------------
mort_post <- read_csv("data/lit/fit/mort_sal_logistic_post.csv")
sal_seq <- seq(5, 35, by=0.1)

plot(NA, NA, xlim=c(5, 35), ylim=c(0, 1),
     xlab="Salinity (psu)", ylab="Hourly larval mortality rate")
for(i in 1:300) {
  lines(sal_seq, 
        mort_post$mortMax[i] / 
          (1+exp(-mort_post$k[i]*(sal_seq - mort_post$salMid[i]))) + 
          mort_post$mortMin[i], 
        col=rgb(27, 158, 119, 25, maxColorValue=256))
}
points(sal_df$salinity, sal_df$mort_h, col="black", pch=19)
abline(h=0.01, lty=3)
legend("topright", 
       c("Logistic", "Constant (0.01)"), 
       col=c(rgb(27, 158, 119, maxColorValue=256),
             "black"),
       lwd=1, 
       lty=c(1,3))

# . posterior mean --------------------------------------------------------
mort_post |> 
  summarise(across(everything(), mean))


