#----------------------------------
# Purpose of the script: Checking modified flexTPC function
#
# What this script does:
# 1. Writes two models to visualise the flexTPC function, one with the inv.logit to create a smooth transition outside of tmin and tmax range, and one without the inv.logit to create a hard cut-off at tmin and tmax
# 2. Checks that the two models give similar results across a range of parameter values, and visualises any differences between the two models
# 3. When there are large differences between the models, try and visualise how big these are and why
#
# Author: Dr. Daniel Padfield
# Date Created: 10 April 2026
#
# Daniel Padfield, 2026
# This code is licensed under an MIT license and has been developed with the help of GitHub copilot, and potentially other LLMs. However, all output has been checked and modified by a human and we take full responsibility for any errors.
#
#----------------------------------
#
# Notes (potentially about software used:
#
# Generally there are not big differences between the two model formulations, but the inv.logit versions allow better convergence in a Bayesian framework.
#
#----------------------------------
# if librarian is not installed, install it
if (!requireNamespace('librarian', quietly = TRUE)) {
  install.packages('librarian')
}
# if BiocManager is not installed, install it
if (!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager')
}
# if Biobase is not installed, install it from Bioconductor
if (!requireNamespace('Biobase', quietly = TRUE)) {
  BiocManager::install('Biobase')
}
# load packages
librarian::shelf(tidyverse, patchwork)
## ---------------------------

# write custom functions

# we present two versions of the flextpc function.
# the first, flextpc_2024_mod, is the one we use in the manuscript. It has inv.logit terms that create a smoother transition outside of the range of tmin and tmax. This seems to help with convergence in brms, especially in hierarchical models.
# we use a value of 25 as it seems to give a steep transition, which allows the parameters of the model to retain their meaning, while still helping the fitting of complex models. If the value is too low, then the transition is gradual and the same model parameters give rise to different curves compared to the original model formulation.
# the second, flextpc_2024_mod2, is equivalent to the original model formulation, where rates outside of tmin and tmax are 0. This model still includes the transformations we used to get the parameter estimates on the same magnitude as each other, as this also helps convergence in brms.

# flextpc modified function without the inv.logit to create a smooth transition outside of tmin and tmax range
# the size of value (here 25) controls how steep the transition is.
flextpc_2024_mod <- function(temp, tmin, tmax, rmax, alpha, beta) {
  # transform parameters back
  tmax <- tmax / 100
  rmax <- log(rmax)
  alpha <- alpha / 10
  tmin <- tmin / 100

  # calculate kappa
  kappa <- (boot::inv.logit(alpha) * (1 - boot::inv.logit(alpha)) / beta^2)

  out <-
    boot::inv.logit(25 * (temp - (tmin * 100))) *
    boot::inv.logit(25 * ((tmax * 100) - temp)) *
    exp(rmax) *
    exp(
      kappa *
        (boot::inv.logit(alpha * 10) *
          log(pmax(temp - (tmin * 100), 10^-20) / boot::inv.logit(alpha * 10)) +
          (1 - boot::inv.logit(alpha * 10)) *
            log(
              pmax((tmax * 100) - temp, 10^-20) /
                (1 - boot::inv.logit(alpha * 10))
            ) -
          log((tmax * 100) - (tmin * 100)))
    )

  return(out)
}


# flextpc modified with a hard cut-off at tmin and tmax, so rates beyond these are zero
flextpc_2024_mod2 <- function(temp, tmin, tmax, rmax, alpha, beta) {
  # transform parameters back
  tmax <- tmax / 100
  rmax <- log(rmax)
  alpha <- alpha / 10
  tmin <- tmin / 100

  # calculate kappa
  kappa <- (boot::inv.logit(alpha) * (1 - boot::inv.logit(alpha)) / beta^2)

  out <-
    ifelse(temp <= (tmin * 100), 0, 1) *
    ifelse(temp >= (tmax * 100), 0, 1) *
    exp(rmax) *
    exp(
      kappa *
        (boot::inv.logit(alpha * 10) *
          log(pmax(temp - (tmin * 100), 10^-20) / boot::inv.logit(alpha * 10)) +
          (1 - boot::inv.logit(alpha * 10)) *
            log(
              pmax((tmax * 100) - temp, 10^-20) /
                (1 - boot::inv.logit(alpha * 10))
            ) -
          log((tmax * 100) - (tmin * 100)))
    )

  return(out)
}

# we can see what the inv.logit terms are doing. They are replacing the hard cut-off at tmin and tmax with a gradual transition, but it is essentially the same thing as making rates beyond tmin and tmax zero

# set tmin and tmax
tmin <- 10
tmax <- 30

# create a dataset to check the behaviour of the inv.logit terms
d_check_invlogit <- data.frame(temp = seq(0, 50, by = 0.1)) %>%
  mutate(
    rate = boot::inv.logit(25 * (temp - tmin)) *
      boot::inv.logit(25 * (tmax - temp)),
    rate2 = ifelse(temp <= tmin, 0, 1) * ifelse(temp >= tmax, 0, 1)
  )

# plot them
p1 <- ggplot(d_check_invlogit, aes(x = temp, y = rate)) +
  geom_line() +
  theme_bw(base_size = 14) +
  ggtitle("With inv.logit with a steepness of 25")
p2 <- ggplot(d_check_invlogit, aes(x = temp, y = rate2)) +
  geom_line() +
  theme_bw(base_size = 14) +
  ggtitle("Make rates zero beyond tmin and tmax")

p1 + p2

# run check
# we can see whether these two models differ to each other by passing them the same parameter values
d_check <- tibble(temp = seq(0, 50, by = 0.5)) %>%
  mutate(
    rate = flextpc_2024_mod(
      temp,
      tmin = 0,
      tmax = 35,
      rmax = exp(0.6),
      alpha = 1.4,
      beta = rlnorm(1, -1.5, 0)
    ),
    rate2 = flextpc_2024_mod2(
      temp,
      tmin = 0,
      tmax = 35,
      rmax = exp(0.6),
      alpha = 1.4,
      beta = rlnorm(1, -1.5, 0)
    )
  )

# check curves are similar, the lines should overlap
ggplot(d_check, aes(x = temp)) +
  geom_line(aes(y = rate), color = 'blue') +
  geom_line(aes(y = rate2), color = 'red') +
  theme_bw(base_size = 14)

# plot rates against each other, should be a 1:1 relationship if they are the same
ggplot(d_check, aes(rate, rate2)) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  geom_point(size = 3) +
  theme_bw(base_size = 14)

# we can simulate a bunch of parameter values that give random curve shapes to see how the two models differ.

# how many curves to estimate
n_curves <- 1000

# we will use the priors as used in 07_brms_tpc.R, but with slightly wider SDs. These are then back-transformed as we use raw values for these models.
alpha_vals <- rnorm(n_curves, 0.14, 0.1) * 10
beta_vals <- rlnorm(n_curves, -1.5, 0.5)
rmax_vals <- rnorm(n_curves, 1, 1.2) %>% exp()
tmin_vals <- rnorm(n_curves, 0, 0.07) * 100
tmax_vals <- rnorm(n_curves, 0.35, 0.2) * 100

# stack all these random values together and apply each function to them
d_check2 <- tibble(alpha_vals, beta_vals, rmax_vals, tmin_vals, tmax_vals) %>%
  mutate(temp = list(seq(min(tmin_vals) - 5, max(tmax_vals) + 5, by = 0.5))) %>%
  unnest(temp) %>%
  group_by(alpha_vals, beta_vals, rmax_vals, tmin_vals, tmax_vals) %>%
  mutate(
    rate = flextpc_2024_mod(
      temp,
      tmin = tmin_vals,
      tmax = tmax_vals,
      rmax = rmax_vals,
      alpha = alpha_vals,
      beta = beta_vals
    ),
    rate2 = flextpc_2024_mod2(
      temp,
      tmin = tmin_vals,
      tmax = tmax_vals,
      rmax = rmax_vals,
      alpha = alpha_vals,
      beta = beta_vals
    ),
    id = cur_group_id()
  ) %>%
  ungroup()

# check that there are not places where rate failed but rate2 didn't fail, or vice versa
d_check2 %>%
  filter((is.na(rate) & !is.na(rate2)) | (!is.na(rate) & is.na(rate2)))
# should be zero

# check that they both failed when one failed
d_check2 %>%
  filter(is.na(rate) & is.na(rate2)) %>%
  nrow()
# yes

# randomly sample 30 curves to plot
d_check2 %>%
  filter(id %in% sample(unique(id), 30)) %>%
  ggplot(., aes(x = temp)) +
  geom_line(aes(y = rate), color = 'blue') +
  geom_line(aes(y = rate2), color = 'red') +
  facet_wrap(~id, scales = "free_y") +
  theme_bw(base_size = 14)
# most of them look identical

# plot rates against each other, should be a 1:1 relationship if they are the same
d_check2 %>%
  ggplot(aes(rate, rate2)) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  geom_point(size = 3) +
  theme_bw(base_size = 14)
# looks really good

# calculate different between the two rates at each temperature
d_check2 <- d_check2 %>%
  mutate(diff = abs(rate - rate2))

# filter for when the max difference between two of the rates are more than 10% the value of rmax
# this will allow us to pinpoint any issues and should be visible
d_check3 <- d_check2 %>%
  group_by(id) %>%
  mutate(max_diff = max(diff, na.rm = TRUE)) %>%
  filter(max_diff > 0.1 * rmax_vals) %>%
  ungroup()
unique(d_check3$id) %>% length()

# plot some of these
d_check3 %>%
  filter(
    id %in% sample(unique(id), size = unique(d_check3$id) %>% length())
  ) %>%
  ggplot(., aes(x = temp)) +
  geom_line(aes(y = rate), color = 'blue') +
  geom_line(aes(y = rate2), color = 'red') +
  facet_wrap(~id, scales = "free_y") +
  theme_bw(base_size = 14)

# these look more or less equivalent
max(d_check3$max_diff / d_check3$rmax_vals)

# lets look at the really big differences
d_check3 %>%
  filter(max_diff > 0.2 * rmax_vals) %>%
  ggplot(., aes(x = temp)) +
  geom_line(aes(y = rate), color = 'blue') +
  geom_line(aes(y = rate2), color = 'red') +
  facet_wrap(~id, scales = "free_y") +
  theme_bw(base_size = 14)
# they still mainly the same
