#----------------------------------
# Purpose of the script: Fit a dose-response curve to kanamycin data using a Bayesian mixed effects model
#
# What this script does:
# 1. Reads in the data and does preprocessing
# 2. Set up model and grab priors from a single curve bayesnec of all the data
# 3. Set priors on population parameters and group-level effects
# 4. Sample prior to check priors are sensible
# 5. Fit full model to data, using ecxwb1 model from bayesnec
# 6. Calculates key traits from the model (EC50, EC95/MIC, NSEC, sub-MIC selection window)
# 7. Analyses these traits
# 8. Analyses estimates of selection
# 9. Makes Figures 2 and 3
#
# Author: Dr. Daniel Padfield
# Date Created: 15 December 2025
#
# Daniel Padfield, 2025
# This code is licensed under an MIT license and has been developed with the help of GitHub copilot, and potentially other LLMs. However, all output has been checked and modified by a human and we take full responsibility for any errors.
#
#----------------------------------
#
# Notes (potentially about software used:
#
#
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
librarian::shelf(
  posterior,
  brms,
  tidybayes,
  bayesplot,
  patchwork,
  bayesnec,
  future,
  monochromeR,
  flextable,
  ggarrow,
  tidyverse,
  cmdstanr
)
## ---------------------------

# set seed for reproducibility
set.seed(42)

# load in custom functions
source('scripts/00_functions.R')

# read in metric data
d_metrics <- read.csv(
  'data/processed/05_kanamycin_metrics_analysis.csv'
) %>%
  mutate(id3 = paste(id2, temp, sep = '_'))
head(d_metrics)
unique(d_metrics$classification)

# are there any growth rate values that are exactly the same
dupes <- d_metrics %>% count(gr) %>% filter(n > 1) %>% filter(gr > 0)

# the same file is in twice
d_metrics %>%
  filter(gr %in% dupes$gr) %>%
  select(file, serial_no) %>%
  distinct()

# remove 21ºC data because it only has 2 reps per strain because it was a mistake when collecting the data. This is unlikely to be enough to estimate the random effect.
d_metrics <- filter(d_metrics, temp != 21)

# going to use the beta family for the model - between 0 and 1
# this was recommended by Rebecca Fisher, an example is here: https://academic.oup.com/ieam/article/20/1/279/7725045 in case study 1
# divide through by the highest observed growth rate and add a small constant if any values are 0 or 1
# convert concentration to asinh scale to get more uniform distance between values
# asinh is inverse hyperbolic sine, similar to log but works with 0 values
# if any values are 0 or 1, add a small constant as beta family needs to be between 0 and 1.
d_metrics <- mutate(
  d_metrics,
  gr2 = gr / max(gr),
  gr2 = ifelse(gr2 == 1, gr2 - 10^-6, gr2),
  gr2 = ifelse(gr2 == 0, gr2 + 10^-6, gr2),
  conc = asinh(conc_µl),
  temp = as.character(temp)
)

# quick plot
ggplot(d_metrics, aes(asinh(conc_µl), gr2, col = id2)) +
  geom_point() +
  facet_wrap(~ id2 + temp, labeller = labeller(.multi_line = FALSE)) +
  theme_bw()

#--------------#
# Modelling ####
#--------------#

# test fit with just a single curve
test <- filter(d_metrics, id3 == 'wt_37')

# lets just fit a single model to the data to grab the default priors out
# using bayesnec out of the box, ecxwb1 is kind of arbitrarily chosen because it fits most curve shapes ok
# just run on a single chain as just testing
mod_test <- bnec(
  bnf(gr2 ~ crf(conc, model = "ecxwb1")),
  data = test,
  open_progress = TRUE,
  iter = 3000,
  chains = 1,
  init = '0'
)

mod_test$fit

# grab the priors out
bayesnec::pull_prior(mod_test)

# look at kanamycin concentrations
test$conc_µl %>% unique() %>% sort() %>% asinh()
test$conc_µl %>% unique() %>% sort() %>% asinh() %>% mean

# write brms formula
# use pretty basic random effect specification of combined temp and isolate
brms_form <- bf(
  gr2 ~ bot + (top - bot) * exp(-exp(exp(beta) * (conc - ec50))),
  brms::nlf(beta ~ betatrans * -10),
  brms::nlf(ec50 ~ ec50trans * 10),
  bot ~ 1 + (1 | id3),
  ec50trans ~ 1 + (1 | id3),
  top ~ 1 + (1 | id3),
  betatrans ~ 1 + (1 | id3),
  nl = TRUE
)

# set priors for full model
# use default beta, bot, top priors from bayesnec
# set ec50 prior based on mean of the concentrations on the asinh scale and a wide-ish sd. Limits set to upper and lower concentrations (on the asinh scale) respectively
priors <- c(
  prior(normal(0, 0.5), nlpar = "betatrans"),
  prior(beta(2, 5), nlpar = "bot", lb = 0, ub = 1),
  prior(beta(5, 2), nlpar = "top", lb = 0, ub = 1),
  prior(normal(0.275, 0.1), nlpar = 'ec50trans', lb = 0.001)
)

# set random priors
# wide-ish priors for bot, top, and ec50
# narrower on beta as not expecting this to vary as much
# default phi prior from bayesnec
priors_random <- random_priors <- c(
  set_prior("student_t(3,0,0.25)", class = "sd", nlpar = "bot"),
  set_prior("student_t(3,0,0.25)", class = "sd", nlpar = "top"),
  set_prior("student_t(5,0,0.25)", class = "sd", nlpar = "ec50trans"),
  set_prior("student_t(5,0,0.25)", class = "sd", nlpar = "betatrans"),
  set_prior("gamma(0.01, 0.01)", class = 'phi')
)

# combine priors
priors_all <- c(priors, priors_random)

table_priors <- priors_all %>%
  data.frame() %>%
  mutate(nlpar = ifelse(nlpar == '', class, nlpar)) %>%
  mutate(nlpar = ifelse(class == 'sd', paste(nlpar, 'sd', ' '), nlpar)) %>%
  select(parameter = nlpar, prior, lb, ub) %>%
  flextable() %>%
  set_header_labels(
    parameter = 'Parameter',
    prior = 'Prior',
    lb = 'Lower boundary',
    ub = 'Upper boundary'
  ) %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 14, part = 'all') %>%
  align(align = 'center', part = 'all') %>%
  autofit(add_w = 0.1) %>%
  fix_border_issues()

# save this out
save_as_image(table_priors, 'tables/table_priors_drc.png', webshot = 'webshot2')
save_as_docx(table_priors, path = 'tables/table_priors_drc.docx')

# set controls
controls <- list(adapt_delta = 0.99, step_size = 0.001, max_treedepth = 15)

# sample prior only
model_prior <- brm(
  brms_form,
  data = d_metrics,
  prior = priors_all,
  chains = 4,
  cores = 4,
  iter = 5000,
  family = 'beta',
  sample_prior = 'only',
  backend = 'cmdstan'
)
# lots of divergent transitions even for the prior, but this doesnt seem to impact the final model
# should really check this

# plot priors of population level effects
# grab draws, `b` usually determines fixed or population level effects
d_priors <- spread_draws(model_prior, `b_.*`, regex = TRUE) %>%
  # remove b_ and _Intercept
  rename_with(~ gsub("b_|_Intercept", "", .x), starts_with("b_")) %>%
  mutate(beta = betatrans * -10, ec50 = ec50trans * 10) %>%
  select(-betatrans, -ec50trans) %>%
  pivot_longer(
    cols = c(top, bot, beta, ec50),
    names_to = 'parameter',
    values_to = 'value'
  )

# plot the priors
ggplot(d_priors, aes(value)) +
  geom_histogram(fill = 'white', col = 'black') +
  facet_wrap(
    ~parameter,
    scales = 'free',
    labeller = labeller(parameter = letter_facets)
  ) +
  theme_bw(base_size = 18) +
  labs(x = 'Parameter value', y = 'Count') +
  theme(
    legend.position = c(0.85, 0.25),
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0)
  )
# these should look good. Do not move on if not

ggsave('plots/growth_drc_priors.png', width = 8, height = 6)

# check priors for random effects
d_priors_random <- spread_draws(model_prior, `sd_.*`, regex = TRUE) %>%
  # remove b_ and _Intercept
  rename_with(~ gsub("sd_id3__|_Intercept", "", .x), starts_with("sd_")) %>%
  mutate(beta = betatrans * 10, ec50 = ec50trans * 10) %>%
  select(-betatrans, -ec50trans) %>%
  pivot_longer(
    cols = c(top, bot, beta, ec50),
    names_to = 'parameter',
    values_to = 'value'
  )
# grab those for phi
d_priors_phi <- spread_draws(model_prior, `phi`, regex = TRUE) %>%
  rename(value = phi) %>%
  mutate(parameter = 'phi')

# combine them and plot
ggplot(bind_rows(d_priors_random, d_priors_phi), aes(value)) +
  geom_histogram(fill = 'white', col = 'black') +
  facet_wrap(
    ~parameter,
    scales = 'free',
    labeller = labeller(parameter = letter_facets)
  ) +
  theme_bw(base_size = 18) +
  labs(x = 'Parameter value', y = 'Count') +
  theme(
    legend.position = c(0.85, 0.25),
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0)
  )

# save plot out
ggsave('plots/growth_drc_priors_random.png', width = 10, height = 6)

# fit full model

# probably do not run this on a laptop as it'll take the best part of an hour
# read in model if you do not want to run it
model_full <- readRDS(
  'data/models/brms_drc_kanamycin_full.rds'
)

model_full <- brm(
  brms_form,
  data = d_metrics,
  prior = priors_all,
  chains = 4,
  cores = 4,
  iter = 20000,
  init = 0,
  family = 'beta',
  backend = 'cmdstan'
)

# save model out
saveRDS(model_full, 'data/models/brms_drc_kanamycin_full.rds')

# check model fit
model_full
pp_check(model_full, ndraws = 100)
plot(model_full)

# plot predictions
# quick method
conditions <- data.frame(id3 = unique(d_metrics$id3))
rownames(conditions) <- unique(d_metrics$id3)

preds_drc <- conditional_effects(
  model_full,
  conditions = conditions,
  re_formula = NULL,
  method = "posterior_epred",
  niter = 500
)

plot(preds_drc, ncol = 5, points = TRUE)

# more controlled method
# get predictions across 200 concentrations
preds_drc <- select(d_metrics, id3) %>%
  distinct() %>%
  group_by(id3) %>%
  do(data.frame(
    conc = seq(min(d_metrics$conc), max(d_metrics$conc), length.out = 200)
  )) %>%
  add_epred_draws(
    model_full,
    re_formula = NULL,
    ndraws = 1000
  ) %>%
  group_by(id3, conc) %>%
  median_hdci() %>%
  separate(id3, into = c('id2', 'temp'), sep = '_', remove = FALSE) %>%
  # back transform concentrations to the original scale
  mutate(
    conc_µl = sinh(conc),
    strain = factor(
      id2,
      levels = c("wt", 'nal.5', "R1", "RP4")
    ),
    temp2 = as.numeric(temp)
  )

d_metrics <- mutate(
  d_metrics,
  strain = factor(
    id2,
    levels = c("wt", 'nal.5', "R1", "RP4")
  ),
  temp2 = as.numeric(temp)
)

# plot predictions
ggplot() +
  geom_point(aes(conc_µl, gr2, group = temp, col = temp2), data = d_metrics) +
  geom_ribbon(
    aes(conc_µl, ymin = .lower, ymax = .upper, fill = temp2, group = temp),
    data = preds_drc,
    alpha = 0.3,
    show.legend = FALSE
  ) +
  geom_line(
    aes(conc_µl, .epred, col = temp2, group = temp),
    data = preds_drc,
    show.legend = FALSE
  ) +
  facet_wrap(~strain) +
  theme_bw(base_size = 14) +
  labs(x = 'Kanamycin concentration (µl)', y = 'Relative growth rate') +
  scale_color_viridis_c('Temperature (°C)') +
  scale_fill_viridis_c()

# from now on remove nal.5 from the dataset as we do not have plasmid data for this strain
# first calculate max growth in the dataset against which to calculate MIC
max_growth <- max(d_metrics$gr2)
d_metrics <- filter(d_metrics, id2 != 'nal.5')
preds_drc <- filter(preds_drc, id2 != 'nal.5')

#--------------------------------------------#
# look at posterior vs prior distribution ####
#--------------------------------------------#

# plot priors of population level effects
# grab draws, `b` usually determines fixed or population level effects
d_post <- spread_draws(model_full, `b_.*`, regex = TRUE) %>%
  # remove b_ and _Intercept
  rename_with(~ gsub("b_|_Intercept", "", .x), starts_with("b_")) %>%
  mutate(beta = betatrans * -10, ec50 = ec50trans * 10) %>%
  select(-betatrans, -ec50trans) %>%
  pivot_longer(
    cols = c(top, bot, beta, ec50),
    names_to = 'parameter',
    values_to = 'value'
  ) %>%
  mutate(type = 'posterior') %>%
  bind_rows(
    d_priors %>%
      mutate(type = 'prior')
  )

# plot the priors and posteriors
ggplot(
  d_post,
  aes(value, col = type, group = interaction(.chain, type))
) +
  geom_density(fill = NA, key_glyph = draw_key_rect) +
  facet_wrap(
    ~parameter,
    scales = 'free',
    labeller = labeller(parameter = letter_facets)
  ) +
  scale_color_manual(name = '', values = c('red', 'black')) +
  theme_bw(base_size = 18) +
  labs(x = 'Parameter value', y = 'Count') +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0)
  )


# save plot out
ggsave('plots/growth_drc_priors_posteriors.png', width = 10, height = 6)

# ec50 has moved a lot and is at the boundary. This makes sense when we think about the plasmids being resistant. Their EC50 must be at concentrations we have not measured.

# look at what parameters exist
tidybayes::get_variables(model_full)

# grab random/group level sd parameters
d_post_pop <- spread_draws(model_full, `sd_.*`, regex = TRUE) %>%
  # remove b_ and _Intercept
  rename_with(~ gsub("sd_id3__|_Intercept", "", .x), starts_with("sd_")) %>%
  mutate(beta = betatrans * 10, ec50 = ec50trans * 10) %>%
  select(-betatrans, -ec50trans) %>%
  pivot_longer(
    cols = c(top, bot, beta, ec50),
    names_to = 'parameter',
    values_to = 'value'
  ) %>%
  mutate(type = 'posterior') %>%
  bind_rows(
    mutate(d_priors_random, type = 'prior')
  )
d_post_phi <- spread_draws(model_full, `phi`, regex = TRUE) %>%
  rename(value = phi) %>%
  mutate(parameter = 'phi', type = 'posterior') %>%
  bind_rows(
    mutate(d_priors_phi, type = 'prior')
  )

# plot them
bind_rows(d_post_pop, d_post_phi) %>%
  ggplot(aes(value, col = type, group = interaction(.chain, type))) +
  geom_density(fill = NA, key_glyph = draw_key_rect) +
  facet_wrap(
    ~parameter,
    scales = 'free',
    labeller = labeller(parameter = letter_facets)
  ) +
  scale_color_manual(name = '', values = c('red', 'black')) +
  theme_bw(base_size = 18) +
  labs(x = 'Parameter value', y = 'Count') +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0),
    legend.position = c(0.85, 0.25)
  )
# looking good!

# save out plot
ggsave('plots/growth_drc_priors_posteriors_random.png', width = 10, height = 6)

#-------------------------------------------#
# calculate traits for each level of id3 ####
#-------------------------------------------#

# for each strain x temperature combination, we want to calculate the no significant effect concentration (NSEC), MIC, and EC50 as common traits which are used to quantify selection.
# NSEC is the concentration at which the growth rate is significantly different from the growth rate at 0 concentration
# EC50 is the concentration at which growth rate is reduced by 50% relative to the maximum growth rate for the strain x temperature combination.
# MIC is the minimum inhibitory concentration, which is the concentration at which no growth occurs. We use EC95 as a proxy for this, and specifically EC95 relative to the maximum growth rate across all the data as this makes the value comparable across all of our curves (e.g.) Our definition of zero growth on a plate reader should not change between strains and temperatures.
# sub-MIC selection window is the difference between MIC and NSEC

# we have to write a few functions to calculate these

# write a function to calculate NSEC
# returns the concentration of the rate closest to the estimated nsec rate
get_nsec <- function(x, y, nsec_rate) {
  conc_nsec <- x[which.min(abs(
    y - nsec_rate
  ))]
  return(conc_nsec)
}

# write function to get the rate at which ECx occurs
# fmt: skip
# z is the growth rate against which you want to calculate absolute ECx, likely to be 1, but could be whatever else you fancy
# you can choose whether you want to set a value for ECx (using type = absolute and giving z), or whether you want the relative ECx values which are set to the max growth rate in that data grouping
# y is the predicted rates
# ec_value is the x in ECx, between 1 and 100.
get_ecx_rate <- function(y, ec_value, type = c("relative", "absolute"), z = NULL) {
  # set type to relative if not specified
  if (missing(type)) {
    type = 'relative'
  }

  if (type == 'relative') {
    # find maximum rate
    y_max <- max(y)
    y_threshold <- y_max * (1 - (ec_value / 100)) %>% mean()
    return(y_threshold)
  }

  if (type == 'absolute') {
    y_threshold <- z * (1 - (ec_value / 100)) %>% mean()
    return(y_threshold)
  }
}

# write function to get the ECx concentration
# this finds the concentration at which the rate is closest to the ECx rate, which is calculated using the get_ecx_rate function.
# fmt: skip
get_ecx <- function(x, y, ec_value, type = c('relative', 'absolute'), z = NULL) {
  if (missing(z)) {
    z = NULL
  }
  if (missing(type)) {
    type = 'relative'
  }

  # find the concentration closest to the threshold rate
  conc_ecx <- x[which.min(abs(
    y - get_ecx_rate(y, ec_value, type = type, z = z)
  ))]

  return(conc_ecx)
}

# calculate NSEC, EC50, EC95/MIC, and sub-MIC window for each id3
d_traits_preds <- select(d_metrics, id3) %>%
  filter(!str_detect(id3, 'nal.5')) %>%
  distinct() %>%
  group_by(id3) %>%
  # create new high resolution dataframe of concentrations
  do(data.frame(
    conc = seq(min(d_metrics$conc), max(d_metrics$conc), length.out = 1000)
  )) %>%
  # add in predicted values, re_formula = NULL to include group level/random effects
  add_epred_draws(
    model_full,
    re_formula = NULL,
    ndraws = 1000
  ) %>%
  # for each id3, calculate the 0.01 quantile of the predicted rate at 0 concentration. This is against which we calculate the NSEC. This is equivalent to a p value of 0.01
  group_by(id3) %>%
  # calculate rates
  mutate(
    rate_nsec = quantile(.epred[conc == min(conc)], 0.01)
  )

d_traits <- d_traits_preds %>%
  group_by(id3, .draw, rate_nsec) %>%
  summarise(
    no_drugrate = .epred[conc == min(conc)],
    conc_nsec = get_nsec(conc, .epred, mean(rate_nsec)),
    conc_mic = get_ecx(
      conc,
      .epred,
      95,
      type = 'absolute',
      z = max(d_metrics$gr2)
    ),
    rate_mic = get_ecx_rate(
      .epred,
      95,
      type = 'absolute',
      z = max(d_metrics$gr2)
    ),
    conc_ec50 = get_ecx(conc, .epred, 50, type = 'relative'),
    rate_ec50 = get_ecx_rate(.epred, 50, type = 'relative'),
    .groups = 'drop'
  ) %>%
  mutate(
    conc_submic_window = sinh(conc_mic) - sinh(conc_nsec),
    rate_submic_window = 0,
  ) %>%
  # reshape dataset
  pivot_longer(
    cols = matches(
      "^(rate|conc)_(nsec|ec50|mic|submic_window)"
    ),
    names_to = c(".value", "trait"),
    names_pattern = "(rate|conc)_(nsec|ec50|mic|submic_window)"
  )

# create copy of this
d_traits_mod <- d_traits

d_traits <- group_by(d_traits, id3, trait) %>%
  # calculate 95% credible intervals
  median_hdci() %>%
  separate(id3, into = c('id2', 'temp'), sep = '_', remove = FALSE) %>%
  mutate(
    conc_µl = ifelse(
      !trait %in% c('submic_window'),
      sinh(conc),
      conc
    ),
    temp2 = as.numeric(temp),
    strain = factor(
      id2,
      levels = c("wt", "R1", "RP4")
    )
  )

# make figure 2 to demonstrate all curve fits and NSEC, EC50, and MIC values

# set colours
cols <- c("#94CBEC", "#C26A77", "#2E2585")
names(cols) <- c('wt', 'R1', 'RP4')

# add in ºC to temperature labels

letter_facets_mod <- function(string) {
  paste(letter_facets(string), '°C')
}

# refactor id2
d_metrics <- mutate(
  d_metrics,
  strain = factor(
    id2,
    levels = c("wt", "R1", "RP4")
  ),
  temp2 = as.numeric(temp)
)

# make plot
ggplot() +
  geom_point(
    aes(conc_µl, gr2, group = temp, col = strain),
    data = d_metrics
  ) +
  geom_ribbon(
    aes(conc_µl, ymin = .lower, ymax = .upper, fill = strain, group = strain),
    data = preds_drc,
    alpha = 0.3,
    show.legend = FALSE
  ) +
  geom_line(
    aes(conc_µl, .epred, col = strain, group = strain),
    data = preds_drc,
    show.legend = FALSE
  ) +
  facet_wrap(~temp, labeller = labeller(temp = letter_facets_mod)) +
  theme_bw(base_size = 14) +
  labs(x = 'Kanamycin concentration (µg/mL)', y = 'Relative growth rate') +
  scale_color_manual(
    values = cols,
    labels = c(
      'wt' = 'wild-type',
      'R1' = 'wild-type + R1',
      'RP4' = 'wild-type + RP4'
    )
  ) +
  scale_fill_manual(
    values = cols,
    labels = c(
      'wt' = 'wild-type',
      'R1' = 'wild-type + R1',
      'RP4' = 'wild-type + RP4'
    )
  ) +
  geom_point(
    aes(conc_µl, rate, shape = trait),
    data = filter(
      d_traits,
      strain == 'wt' & trait %in% c('mic', 'nsec', 'ec50')
    ),
    col = 'black',
    fill = 'white',
    size = 3,
    show.legend = FALSE
  ) +
  scale_shape_manual(values = c(21, 22, 24)) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0),
    legend.position = 'bottom',
    legend.title = element_blank()
  ) +
  NULL

# save out plot
ggsave(
  'plots/figure_2.png',
  width = 10,
  height = 8
)

#---------------------------------------------------#
# how do the traits differ between temperatures? ####
#---------------------------------------------------#

# filter this just for the wild type
d_traits_wt <- filter(d_traits, id2 == 'wt') %>%
  mutate(
    conc.lower = ifelse(
      !trait %in% c('submic_window'),
      sinh(conc.lower),
      conc.lower
    ),
    conc.upper = ifelse(
      !trait %in% c('submic_window'),
      sinh(conc.upper),
      conc.upper
    )
  ) %>%
  mutate(before_topt = ifelse(temp2 < 40.5, 'before', 'after'))

autocontrast <- aes(colour = after_scale(contrast(fill)))

# look at relationship between growth in drug free environment vs NSEC, EC50, and MIC
p_gr_mic <- filter(d_traits_wt, trait %in% c('mic')) %>%
  ggplot(aes(fill = temp2)) +
  geom_linerange(
    aes(
      x = no_drugrate,
      ymin = conc.lower,
      ymax = conc.upper
    ),
    linewidth = 0.4
  ) +
  geom_linerange(
    aes(
      xmin = no_drugrate.lower,
      xmax = no_drugrate.upper,
      y = conc_µl
    ),
    linewidth = 0.4
  ) +
  geom_arrow_chain(
    aes(x = no_drugrate, y = conc_µl, group = 1),
    linewidth = 0.8,
    data = . %>% arrange(temp2),
    resect_head = 3
  ) +
  geom_point(
    aes(
      x = no_drugrate,
      y = conc_µl,
      fill = temp2
    ),
    size = 7,
    shape = 21,
    show.legend = FALSE
  ) +
  geom_text(
    aes(x = no_drugrate, y = conc_µl, label = temp2, !!!autocontrast),
    size = pts(8)
  ) +
  theme_bw(base_size = 14) +
  labs(
    x = 'Relative growth with no kanamycin',
    y = 'MIC (µg/mL)',
    title = '(e) MIC vs. growth rate'
  ) +
  scale_fill_viridis_c('Temperature (°C)') +
  ylim(c(0, 28)) +
  xlim(c(0, 0.9)) +
  theme(plot.title = element_text(size = rel(0.8)))

p_gr_mic

# analyse how MIC, NSEC, and sub-MIC selection window change across temperatures

filter(d_traits_wt, trait %in% c('mic', 'nsec', 'submic_window')) %>%
  ggplot(aes(temp2, conc_µl, col = trait)) +
  geom_point(aes(col = trait), size = 3, position = position_dodge(1)) +
  geom_linerange(
    aes(
      ymin = conc.lower,
      ymax = conc.upper,
      group = interaction(strain, trait)
    ),
    position = position_dodge(width = 1),
    linewidth = 0.4
  ) +
  theme_bw(base_size = 14) +
  facet_wrap(~trait) +
  labs(x = 'Temperature (°C)', y = 'Kanamycin trait value')

# want to fit a quadratic model to these to see how they change with temperature, and see whether they differ between traits.
# interested in (a) position of any peak (b) fold change of each trait across temperatures, and (c) smallest and biggest differences between traits at any given temperature

# to propagate uncertainty, we will use the 500 draws of traits we have and run a model on each of those. We will then process those to get the distribution of the parameters we are interested in.
# this can neatly be done using brm_multiple

# set number of draws
num_draws <- 500

head(d_traits_mod)

# create list of datasets, 1 for each draw
d_traits_mod <- d_traits_mod %>%
  separate(id3, into = c('id2', 'temp'), sep = '_', remove = FALSE) %>%
  mutate(
    conc_µl = ifelse(
      !trait %in% c('submic_window', 'hillcoef'),
      sinh(conc),
      conc
    ),
    temp2 = as.numeric(temp),
    strain = factor(
      id2,
      levels = c("wt", "R1", "RP4")
    )
  ) %>%
  filter(id2 == 'wt') %>%
  dplyr::select(.draw, temp2, trait, conc_µl)

min(d_traits_mod$conc_µl)

d_traits_list <- filter(
  d_traits_mod,
  trait %in% c('mic', 'nsec', 'ec50', 'submic_window')
) %>%
  filter(.draw %in% sample(unique(.draw), num_draws)) %>%
  group_by(.draw) %>%
  group_split()

# fit a non-truncated model to all the data
# fixed effects of temp, trait
fit_test <- brm(
  conc_µl ~ poly(temp2, 2) * trait,
  data = filter(
    d_traits_mod,
    trait %in% c('mic', 'nsec', 'ec50', 'submic_window')
  ),
  iter = 5000,
  warmup = 3500,
  chains = 4,
  cores = 4
)

# set a prior for sigma based on this
# this helps the model to converge
prior <- c(
  prior(normal(3.15, 2), class = "sigma")
)

# create flextable of custom prior
table_priors2 <- prior %>%
  dplyr::select(parameter = class, prior, lb, ub) %>%
  flextable() %>%
  set_header_labels(
    parameter = 'Parameter',
    prior = 'Prior',
    lb = 'Lower boundary',
    ub = 'Upper boundary'
  ) %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 14, part = 'all') %>%
  align(align = 'center', part = 'all') %>%
  autofit(add_w = 0.1) %>%
  fix_border_issues()

save_as_docx(table_priors2, path = 'tables/table_priors_drc_traits.docx')
save_as_image(
  table_priors2,
  'tables/table_priors_drc_traits.png',
  webshot = 'webshot2'
)

# create temporary dataset as needed
d_traits_temp <- d_traits_list[1:10]

# fit truncated normal model

# this takes a long time
# read it in here
fit_trait <- readRDS('data/models/brms_kana_traits.rds')

# set up parallel processing
plan(multisession, workers = 2)

fit_trait <- brm_multiple(
  conc_µl | trunc(lb = 0) ~ poly(temp2, 2) * trait * id2,
  data = d_traits_list,
  prior = prior,
  iter = 5000,
  warmup = 3500,
  chains = 4
)

# check convergence of individual models
d_draws <- posterior::as_draws_df(fit_trait)

# add column for each model
# this number is based on each chain comes through each time
d_draws <- mutate(d_draws, model = rep(1:500, each = 4 * 1500))

# nest the datasets by model and run convergence diagnostics on each model
d_converge <- d_draws %>%
  group_by(model) %>%
  nest() %>%
  mutate(
    diag = map(
      data,
      ~ posterior::summarise_draws(
        .x,
        posterior::default_convergence_measures()
      )
    )
  ) %>%
  select(model, diag) %>%
  unnest(diag)

# look at whether any models have rhats > 1.01
d_converge %>%
  filter(rhat > 1.1)
# how many models have Rhats over 1.01
d_converge %>%
  filter(rhat >= 1.02) %>%
  pull(model) %>%
  unique() %>%
  length()
# generally pretty good, will leave them in because they do contribute to the overall posterior distribution

# save out model
saveRDS(fit_trait, 'data/models/brms_kana_traits.rds')

# set up a new dataset and create predictions
preds_traits <- expand.grid(
  temp2 = seq(
    min(d_traits_mod$temp2),
    max(d_traits_mod$temp2),
    length.out = 100
  ),
  trait = c('nsec', 'ec50', 'mic', 'submic_window')
) %>%
  add_epred_draws(
    fit_trait,
    re_formula = NULL,
    ndraws = 1000
  )

# get summary of these predictions
preds_traits2 <- preds_traits %>%
  group_by(temp2, trait) %>%
  median_hdci()

# calculate summary stats for each trait
trait_pred_sum <- preds_traits %>%
  filter(trait %in% c('nsec', 'ec50', 'mic', 'submic_window')) %>%
  group_by(trait, .draw) %>%
  summarise(
    # where is the trait highest (indicator of tolerance)
    max_temp = temp2[which.max(.epred)],
    # what is the value of the trait at this temperature
    max_epred = max(.epred),
    # where is the trait lowest (indicator of sensitivity)
    min_temp = temp2[which.min(.epred)],
    # what is the value of the trait at this temperature
    min_epred = min(.epred),
    # what is the fold change between the min and max values of the trait across temperatures
    fold_change = max(.epred) / min(.epred),
    .groups = 'drop'
  ) %>%
  ungroup()

#estimates for traits
trait_pred_sum2 <- trait_pred_sum %>%
  pivot_longer(
    cols = contains('_'),
    names_to = 'parameter',
    values_to = 'value'
  ) %>%
  group_by(trait, parameter) %>%
  median_hdci(value)

# check that the identified traits are positioned where we expect them to be.
ggplot() +
  geom_line(aes(temp2, .epred), filter(preds_traits, .draw < 5)) +
  geom_point(aes(max_temp, max_epred), filter(trait_pred_sum, .draw < 5)) +
  geom_point(aes(min_temp, min_epred), filter(trait_pred_sum, .draw < 5)) +
  facet_wrap(~ trait + .draw, labeller = labeller(.multi_line = FALSE))
# yes they are!

# make Figure 3

# reorder traits for plotting
preds_traits2$trait <- factor(
  preds_traits2$trait,
  levels = c('nsec', 'ec50', 'mic', 'submic_window')
)
d_traits_wt$trait <- factor(
  d_traits_wt$trait,
  levels = c('nsec', 'ec50', 'mic', 'submic_window')
)

# uppercase all the trait labels for plotting
levels(preds_traits2$trait) <- toupper(levels(preds_traits2$trait))
levels(d_traits_wt$trait) <- toupper(levels(d_traits_wt$trait))

# change submic_window to sub-MIC window
levels(preds_traits2$trait) <- gsub(
  'SUBMIC_WINDOW',
  'sub-MIC window',
  levels(preds_traits2$trait)
)
levels(d_traits_wt$trait) <- gsub(
  'SUBMIC_WINDOW',
  'sub-MIC window',
  levels(d_traits_wt$trait)
)
unique(d_traits_wt$trait)
unique(preds_traits2$trait)

# make a plot
p_traits <- ggplot(d_traits_wt %>% filter(!is.na(trait))) +
  geom_ribbon(
    aes(temp2, ymin = .lower, ymax = .upper),
    preds_traits2,
    alpha = 0.2
  ) +
  geom_line(aes(temp2, .epred), preds_traits2) +
  geom_linerange(
    aes(
      x = temp2,
      ymin = conc.lower,
      ymax = conc.upper,
      group = interaction(strain, trait)
    ),
    position = position_dodge(width = 1),
    linewidth = 0.4
  ) +
  geom_point(
    aes(temp2, conc_µl, shape = trait),
    size = 3,
    col = 'black',
    fill = 'white',
    show.legend = FALSE
  ) +
  scale_shape_manual(values = c(21, 22, 24, 21)) +
  theme_bw(base_size = 18) +
  facet_wrap(~trait, labeller = labeller(trait = letter_facets)) +
  labs(x = 'Temperature (°C)', y = 'Kanamycin concentration (µg / mL)') +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0)
  )

p_traits

#--------------------------------#
# relative fitness of strains ####
#--------------------------------#

# we can calculate the selection rate constant (SRC) and selection coefficient (SC) of susceptible vs resistant strains
# this is a metric of selection by taking into account the impact of carrying the plasmid and the impact of antibiotics on the susceptible strain.

# we can choose several concentrations to look at. We will choose 0, 1, 5, 10, 15, and 20 µg/mL as these are concentrations that span the NSEC to the MIC across the entire range of temperatures.

# calculate SRC and SC of R1 and RP4 compared to wt across temperatures and concentrations
d_sel <- dplyr::select(d_metrics, id3) %>%
  distinct() %>%
  filter(grepl('wt|R1|RP4', id3)) %>%
  group_by(id3) %>%
  do(data.frame(
    conc = asinh(c(0, 2, 5, 10, 15))
  )) %>%
  add_epred_draws(
    model_full,
    re_formula = NULL,
    ndraws = 1000
  ) %>%
  separate(id3, into = c('id2', 'temp'), sep = '_') %>%
  filter(id2 %in% c('wt', 'R1', 'RP4')) %>%
  pivot_wider(
    id_cols = c(.draw, temp, conc),
    names_from = id2,
    values_from = .epred
  ) %>%
  mutate(
    sel_R1 = (wt - R1) / R1, # SC
    sel_RP4 = (wt - RP4) / RP4, # SC
    sel_R1_2 = wt - R1,
    sel_RP4_2 = wt - RP4
  ) %>%
  dplyr::select(
    .draw,
    temp,
    conc,
    sel_R1,
    sel_RP4,
    sel_R1_2,
    sel_RP4_2,
  ) %>%
  pivot_longer(
    cols = starts_with('sel'),
    names_to = 'plasmid',
    values_to = 'sel',
    names_prefix = 'sel_'
  ) %>%
  mutate(
    temp = as.numeric(temp),
    conc_µl = sinh(conc),
    conc2 = as.character(conc_µl)
  )

# create a summary of these
d_sel_summary <- d_sel %>%
  group_by(temp, conc, plasmid, conc2, conc_µl) %>%
  median_hdci(sel) %>%
  ungroup()

# quick plot
ggplot(d_sel_summary, aes(temp, sel, col = conc2)) +
  geom_hline(yintercept = 0, lty = 2, col = 'grey50') +
  #geom_linerange(aes(ymin = .lower, ymax = .upper), linewidth = 0.4) +
  geom_point(size = 3, position = position_jitter(width = 1)) +
  facet_wrap(~plasmid, scales = 'free_y') +
  theme_bw(base_size = 14) +
  labs(x = 'Temperature (°C)', y = 'Selection coefficient')

# lets look at the concentrations we have and whether they are right
unique(d_sel_summary$conc) %>% sinh() %>% sort()

# are they discrete in conc2?
unique(d_sel_summary$conc2)

# make plot of how selection coefficient changes across temperatures
# dont show all temperatures all the time
# plot selection coefficient

# add in costs of resistance from TPCs

# look at colour palette
cols <- monochromeR::generate_palette(
  "red",
  modification = "go_both_ways",
  n_colours = 5,
  view_palette = TRUE
)
# add in a light grey for when there is 0 antibiotic
cols <- c('light grey', cols[c(1, 2, 4, 5)])
names(cols) <- as.character(c(0, 2, 5, 10, 15))

# Selection rate constant: analysis on density and frequency dependent selection coefficient
# keep 500 draws to fit the model on
draws_to_keep <- sample(unique(d_sel$.draw), 500)

d_src_list <- d_sel %>%
  filter(plasmid %in% c('R1_2', 'RP4_2')) %>%
  filter(.draw %in% draws_to_keep) %>%
  group_by(.draw) %>%
  group_split()

d_src_list_2 <- d_sel %>%
  filter(plasmid %in% c('R1_2', 'RP4_2')) %>%
  filter(.draw %in% draws_to_keep) %>%
  group_by(.draw) %>%
  group_split()

d_src_temp <- d_src_list_2[1:10]

# this takes a long time. Is fitting 500 brms models

# you can read the model object in here:
fit_src <- readRDS('data/models/brms_kana_src.rds')

# set up parallel processing
plan(multisession, workers = 2)

# fit model on four R instances
fit_src <- brm_multiple(
  sel ~ poly(temp, 2) * conc2 * plasmid,
  data = d_src_list_2,
  iter = 3000,
  chains = 4,
  init = 0
)

# check convergence of individual models
d_draws <- posterior::as_draws_df(fit_src)

# add column for each model
# this number is based on each chain comes through each time
d_draws <- mutate(d_draws, model = rep(1:500, each = 4 * 1500))

# nest the datasets by model and run convergence diagnostics on each model
d_converge <- d_draws %>%
  group_by(model) %>%
  nest() %>%
  mutate(
    diag = map(
      data,
      ~ posterior::summarise_draws(
        .x,
        posterior::default_convergence_measures()
      )
    )
  ) %>%
  dplyr::select(model, diag) %>%
  unnest(diag)

# how many models have Rhats over 1.01
d_converge %>%
  filter(rhat >= 1.02) %>%
  pull(model) %>%
  unique() %>%
  length()
# 0 out of 500
# generally pretty good, will leave them in because they do contribute to the overall posterior distribution

# save out model
saveRDS(fit_src, 'data/models/brms_kana_src.rds')

# create predictions
preds_src <- expand.grid(
  temp = seq(min(d_sel$temp), max(d_sel$temp), length.out = 100),
  plasmid = c('R1_2', 'RP4_2'),
  conc2 = unique(d_sel$conc2)
) %>%
  add_epred_draws(
    fit_src,
    re_formula = NULL,
    ndraws = 1000
  )

# get summary of these predictions
preds_src2 <- preds_src %>%
  filter(plasmid %in% c('R1_2', 'RP4_2')) %>%
  group_by(temp, conc2, plasmid) %>%
  median_hdci()

# summary traits for these models:
# find first temperature where CI fully excludes zero (lower bound > 0)
# not sure how you did this but looks right.
src_threshold <- preds_src2 %>%
  group_by(conc2, plasmid) %>%
  summarise(
    {
      candidates <- which(.upper < 0)

      if (length(candidates) == 0) {
        tibble(
          threshold_temp = NA_real_,
          threshold_epred = NA_real_,
          threshold_lower = NA_real_,
          threshold_upper = NA_real_
        )
      } else {
        # Among temperatures whose lower CI > 0,
        # select the lowest temperature
        i <- candidates[which.min(temp[candidates])]

        tibble(
          threshold_temp = temp[i],
          threshold_epred = .epred[i],
          threshold_lower = .lower[i],
          threshold_upper = .upper[i]
        )
      }
    },
    .groups = "drop"
  )
src_threshold

filter(src_threshold, conc2 == 10 & plasmid == 'R1_2')

# make a plot of this point and see if it makes sense
ggplot() +
  geom_hline(yintercept = 0, lty = 2, col = 'grey50') +
  geom_ribbon(
    aes(temp, fill = conc2, ymin = .lower, ymax = .upper),
    alpha = 0.2,
    show.legend = FALSE,
    preds_src2
  ) +
  facet_wrap(~ plasmid + conc2, labeller = labeller(.multi_line = FALSE)) +
  geom_point(
    aes(threshold_temp, threshold_epred),
    src_threshold,
    col = 'black',
    size = 3
  ) +
  geom_point(
    aes(threshold_temp, threshold_lower),
    src_threshold,
    col = 'black',
    size = 2
  ) +
  theme_bw()

# try different summaries
src_pred_estim <- preds_src %>%
  filter(conc2 %in% c('0', '2', '5', '10', '15')) %>%
  group_by(conc2, plasmid, .draw) %>%
  summarise(
    # maximum epred and its temperature
    min_epred = .epred[which.min(.epred)],
    min_temp = temp[which.min(.epred)]
  ) %>%
  ungroup()
src_pred_estim

# estimates for traits
src_pred_sum <- src_pred_estim %>%
  group_by(plasmid, conc2) %>%
  median_hdci()

filter(
  preds_src2,
  abs(temp - 31) == min(abs(temp - 31)) & conc2 == 10 & plasmid == 'R1_2'
)

# make plot
ggplot() +
  geom_hline(yintercept = 0, lty = 2, col = 'grey50') +
  geom_ribbon(
    aes(temp, fill = conc2, ymin = .lower, ymax = .upper),
    alpha = 0.2,
    show.legend = FALSE,
    preds_src2
  ) +
  facet_wrap(~ plasmid + conc2, labeller = labeller(.multi_line = FALSE)) +
  geom_vline(
    aes(xintercept = threshold_temp),
    src_threshold,
    col = 'black'
  ) +
  geom_vline(aes(xintercept = min_temp), src_pred_sum, col = 'black', lty = 2) +
  geom_vline(
    aes(xintercept = min_temp.lower),
    src_pred_sum,
    col = 'black',
    lty = 2
  ) +
  theme_bw()

# create figure
ggplot(preds_src2) +
  geom_hline(yintercept = 0, lty = 2, col = 'grey50') +
  geom_ribbon(
    aes(temp, fill = conc2, ymin = .lower, ymax = .upper),
    alpha = 0.2,
    show.legend = FALSE
  ) +
  geom_line(aes(temp, .epred, col = conc2), show.legend = FALSE) +
  geom_linerange(
    aes(temp, ymin = .lower, ymax = .upper),
    filter(d_sel_summary, plasmid %in% c('R1_2', 'RP4_2'))
  ) +
  geom_point(
    aes(temp, sel, fill = conc2),
    col = 'black',
    size = 3,
    shape = 21,
    filter(d_sel_summary, plasmid %in% c('R1_2', 'RP4_2'))
  ) +
  facet_wrap(~plasmid, labeller = labeller(plasmid = letter_facets)) +
  theme_bw(base_size = 14) +
  labs(
    x = 'Temperature (ºC)',
    y = 'Selection rate constant'
  ) +
  scale_color_manual(values = cols) +
  scale_fill_manual(
    values = cols,
    labels = c(
      '0' = '0µg/mL',
      '2' = '2µg/mL',
      '5' = '5µg/mL',
      '10' = '10µg/mL',
      '15' = '15µg/mL'
    )
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0),
    legend.position = 'bottom',
    legend.title = element_blank()
  )

# use just R1 for Figure 3
p_src_r1 <- filter(preds_src2, str_detect(plasmid, 'R1')) %>%
  ggplot() +
  geom_hline(yintercept = 0, lty = 2, col = 'grey50') +
  geom_ribbon(
    aes(temp, fill = conc2, ymin = .lower, ymax = .upper),
    alpha = 0.2,
    show.legend = FALSE
  ) +
  geom_line(aes(temp, .epred, col = conc2), show.legend = FALSE) +
  geom_linerange(
    aes(temp, ymin = .lower, ymax = .upper),
    filter(d_sel_summary, plasmid %in% c('R1_2'))
  ) +
  geom_point(
    aes(temp, sel, fill = conc2),
    col = 'black',
    size = 3,
    shape = 21,
    filter(d_sel_summary, plasmid %in% c('R1_2'))
  ) +
  theme_bw(base_size = 14) +
  labs(
    x = 'Temperature (ºC)',
    y = 'Selection rate constant',
    title = '(f) R1 selection rate constant'
  ) +
  scale_color_manual(values = cols) +
  scale_fill_manual(
    values = cols,
    labels = c(
      '0' = '0µg/mL',
      '2' = '2µg/mL',
      '5' = '5µg/mL',
      '10' = '10µg/mL',
      '15' = '15µg/mL'
    )
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0),
    legend.position = 'bottom',
    legend.title = element_blank(),
    plot.title = element_text(size = rel(0.8))
  ) +
  guides(fill = guide_legend(nrow = 2, ncol = 3, byrow = TRUE))

# make plot for RP4
filter(preds_src2, str_detect(plasmid, 'RP4')) %>%
  ggplot() +
  geom_hline(yintercept = 0, lty = 2, col = 'grey50') +
  geom_ribbon(
    aes(temp, fill = conc2, ymin = .lower, ymax = .upper),
    alpha = 0.2,
    show.legend = FALSE
  ) +
  geom_line(aes(temp, .epred, col = conc2), show.legend = FALSE) +
  geom_linerange(
    aes(temp, ymin = .lower, ymax = .upper),
    filter(d_sel_summary, plasmid %in% c('RP4_2'))
  ) +
  geom_point(
    aes(temp, sel, fill = conc2),
    col = 'black',
    size = 3,
    shape = 21,
    filter(d_sel_summary, plasmid %in% c('RP4_2'))
  ) +
  theme_bw(base_size = 14) +
  labs(
    x = 'Temperature (ºC)',
    y = 'Selection rate constant',
    title = 'RP4 selection rate constant'
  ) +
  scale_color_manual(values = cols) +
  scale_fill_manual(
    values = cols,
    labels = c(
      '0' = '0µg/mL',
      '2' = '2µg/mL',
      '5' = '5µg/mL',
      '10' = '10µg/mL',
      '15' = '15µg/mL'
    )
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0),
    legend.position = 'bottom',
    legend.title = element_blank()
  )

ggsave('plots/RP4_src.png', height = 6, width = 8)

#----------------------------------------------------------------#
# analyse how selection coefficient changes with temperatures ####
#----------------------------------------------------------------#

# keep 500 draws to fit the model on
draws_to_keep <- sample(unique(d_sel$.draw), 500)

d_sc_list <- d_sel %>%
  filter(plasmid %in% c('R1', 'RP4')) %>%
  filter(.draw %in% draws_to_keep) %>%
  group_by(.draw) %>%
  group_split()

d_sc_temp <- d_sc_list[1:10]

# this takes a long time. Is fitting 500 brms models

# you can read the model object in here:
fit_sc <- readRDS('data/models/brms_kana_sc.rds')

# set up parallel processing
plan(multisession, workers = 2)

# fit model on four R instances
fit_sc <- brm_multiple(
  sel ~ poly(temp, 2) * conc2 * plasmid,
  data = d_sc_list,
  iter = 3000,
  chains = 4,
  init = 0
)

# check convergence of individual models
d_draws <- posterior::as_draws_df(fit_sc)

# add column for each model
# this number is based on each chain comes through each time
d_draws <- mutate(d_draws, model = rep(1:500, each = 4 * 1500))

# nest the datasets by model and run convergence diagnostics on each model
d_converge <- d_draws %>%
  group_by(model) %>%
  nest() %>%
  mutate(
    diag = map(
      data,
      ~ posterior::summarise_draws(
        .x,
        posterior::default_convergence_measures()
      )
    )
  ) %>%
  dplyr::select(model, diag) %>%
  unnest(diag)

# how many models have Rhats over 1.01
d_converge %>%
  filter(rhat >= 1.02) %>%
  pull(model) %>%
  unique() %>%
  length()
# 0 out of 500
# generally pretty good, will leave them in because they do contribute to the overall posterior distribution

# save out model
saveRDS(fit_sc, 'data/models/brms_kana_sc.rds')

# create predictions
preds_sc <- expand.grid(
  temp = seq(min(d_sel$temp), max(d_sel$temp), length.out = 100),
  plasmid = c('R1', 'RP4'),
  conc2 = unique(d_sel$conc2)
) %>%
  add_epred_draws(
    fit_sc,
    re_formula = NULL,
    ndraws = 1000
  )

# get summary of these predictions
preds_sc2 <- preds_sc %>%
  filter(plasmid %in% c('R1', 'RP4')) %>%
  group_by(temp, conc2, plasmid) %>%
  median_hdci()

# try different summaries
sc_pred_estim <- preds_sc %>%
  filter(conc2 %in% c('0', '2', '5', '10', '15')) %>%
  group_by(conc2, plasmid, .draw) %>%
  summarise(
    # maximum absolute epred and its temperature
    max_abs_epred = .epred[which.max(abs(.epred))],
    max_abs_temp = temp[which.max(abs(.epred))],
    # epred closest to zero and its temperature
    min_abs_epred = .epred[which.min(abs(.epred))],
    min_abs_temp = temp[which.min(abs(.epred))]
  ) %>%
  ungroup()
sc_pred_estim

# estimates for traits
sc_preds_sum <- sc_pred_estim %>%
  pivot_longer(
    cols = contains('_'),
    names_to = 'parameter',
    values_to = 'value'
  ) %>%
  group_by(plasmid, conc2, parameter) %>%
  median_hdci(value)
sc_preds_sum

# test for no overlap limits, e.g. where do we see selection for resistance (CI fully exclude zero.)
# we only care about concentration == 5, bc 10 and 15 no overlap with 0
thresholds <- preds_sc2 %>%
  group_by(conc2, plasmid) %>%
  filter(!all(.upper < 0)) %>% # remove groups where CI excludes 0 at every temperature
  arrange(temp, .by_group = TRUE) %>%
  summarise(
    # peak = row where .epred is closest to zero
    peak_epred = .epred[which.min(abs(.epred))],
    peak_temp = temp[which.min(abs(.epred))],

    # low side: highest temp BELOW peak where .upper first becomes negative
    low_temp_threshold = {
      candidates <- temp[temp < temp[which.min(abs(.epred))] & .upper < 0]
      if (length(candidates) == 0) NA_real_ else max(candidates)
    },

    # high side: lowest temp ABOVE peak where .upper first becomes negative
    high_temp_threshold = {
      candidates <- temp[temp > temp[which.min(abs(.epred))] & .upper < 0]
      if (length(candidates) == 0) NA_real_ else min(candidates)
    }
  ) %>%
  left_join(
    preds_sc2 %>%
      select(conc2, plasmid, temp, .epred) %>%
      rename(low_epred = .epred),
    by = c('conc2', 'plasmid', 'low_temp_threshold' = 'temp')
  ) %>%
  left_join(
    preds_sc2 %>%
      select(conc2, plasmid, temp, .epred) %>%
      rename(high_epred = .epred),
    by = c('conc2', 'plasmid', 'high_temp_threshold' = 'temp')
  ) %>%
  ungroup()

# create figure
ggplot(preds_sc2) +
  geom_hline(yintercept = 0, lty = 2, col = 'grey50') +
  geom_ribbon(
    aes(temp, fill = conc2, ymin = .lower, ymax = .upper),
    alpha = 0.2,
    show.legend = FALSE
  ) +
  geom_line(aes(temp, .epred, col = conc2), show.legend = FALSE) +
  geom_linerange(
    aes(temp, ymin = .lower, ymax = .upper),
    filter(d_sel_summary, plasmid %in% c('R1', 'RP4'))
  ) +
  geom_point(
    aes(temp, sel, fill = conc2),
    col = 'black',
    size = 3,
    shape = 21,
    filter(d_sel_summary, plasmid %in% c('R1', 'RP4'))
  ) +
  facet_wrap(~plasmid, labeller = labeller(plasmid = letter_facets)) +
  theme_bw(base_size = 14) +
  labs(
    x = 'Temperature (ºC)',
    y = 'Selection coefficient'
  ) +
  scale_color_manual(values = cols) +
  scale_fill_manual(
    values = cols,
    labels = c(
      '0' = '0µg/mL',
      '2' = '2µg/mL',
      '5' = '5µg/mL',
      '10' = '10µg/mL',
      '15' = '15µg/mL'
    )
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0),
    legend.position = 'bottom',
    legend.title = element_blank()
  )

# save out plot
ggsave(
  'plots/figure_selection_coefficient.png',
  width = 7,
  height = 4
)

#------------------#
# make figure 3 ####
#------------------#

p_traits /
  (p_gr_mic + p_src_r1) +
  plot_layout(heights = c(0.67, 0.33)) &
  theme(text = element_text(size = 14), plot.title = element_text(size = 14))

ggsave(
  'plots/figure_3.png',
  width = 8,
  height = 10
)
