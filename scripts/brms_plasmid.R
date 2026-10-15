#----------------------------------
# Purpose of the script: Fit thermal performance curves to R1 and RP4 plasmid transfer data
#
# What this script does:
# 1. Reads in conjugation data
# 2. Estimates conjugation rate using t/r.d.time and compares it to the Approximate Simonsen method
# 3. Fits thermal performance curves to t/r.d.time gamma estimates for R1 and RP4 separately
#
# Author: Dr. Daniel Padfield
# Date Created: 08 April 2026
#
# Daniel Padfield, 2026
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
# fmt: skip
librarian::shelf(
  readxl,
  JSHuisman/conjugator,
  paletteer,
  patchwork,
  brms,
  tidybayes,
  ggdist,  
  bayesplot,
  tidyverse,
  flextable
)
## ---------------------------

# set seed for reproducibility
set.seed(42)

# load in custom functions
source('scripts/00_functions.R')

# read in datset containing endpoint densities from conjugation assays
conjugation <- read_xlsx(
  'data/raw/conjugation/conjugation.xlsx',
  sheet = 1
)

# remove rep 7 for RP4 because some of the transconjugants were uncountable
conjugation <- conjugation %>%
  filter(!(rep == 7 & strain == "RP4"))

# transform variables
conjugation <- conjugation %>%
  mutate(
    dilution = as.numeric(dilution),
    temp = as.numeric(temp),
    rep = as.character(rep),
  )


unique(conjugation$strain)
# keep only R1 and RP4. E.g. remove controls

conjugation <- conjugation %>%
  filter(strain %in% c("R1", "RP4"))

# estimate gamma for plasmid transfer rate
conjugation_wide <- conjugation %>%
  dplyr::select(strain, plate, dilution, rep, time_h, temp, density, date) %>%
  # separate datasets for each selective plate
  group_by(strain, rep, time_h, temp, date) %>%
  # create donor, recipient, transconjugant columns
  summarise(
    trans = density[plate == "T"],
    recip = density[plate == "R"],
    donor_plus_trans = density[plate == "DT"],
    .groups = "drop"
  ) %>%
  # calculate donor density and conjugation rate: gamma
  mutate(
    donor = donor_plus_trans - trans,
    gamma = trans / (recip * donor * time_h)
  )
conjugation_wide

# quick plot
ggplot(conjugation_wide, aes(temp, gamma)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~strain, scale = 'free_y')

#------------------------------------------------#
# estimate conjugation rate using ASM method  ####
#------------------------------------------------#

# lets see how our estimates change if we estimate transfer rate accounting for growth rates of donor, recipient and transconjugant. We shall do this using the Bayesian model from brms_tpc.R.R

# load in model
model_growth <- readRDS(
  'data/models/brms_flextpc.rds'
)

# define stan function not included in R to allow for R to be able to process the model
fmax <- pmax

# get growth rate estimates for the donor (d), recipient (r), and transconjugant (t) strains at each temperature we measured plasmid transfer
temp_grid <- unique(conjugation_wide$temp)

# the ids the brms model to keep which are our d, r, and t
ids_keep <- c("wt + R1", "nal.5", "wt + RP4", "nal.5 + R1", 'nal.5 + RP4')

# create new data for predictions of growth rates for our d, r, and t
newdata <- expand.grid(
  temp = temp_grid,
  id2 = ids_keep
) %>%
  add_epred_draws(model_growth, re_formula = NULL, ndraws = 1000)

# get credible intervals for these
pred_summary <- newdata %>%
  group_by(temp, id2) %>%
  median_qi()
# the median estimate is what we will use for the estimate of ASM

# make plot of predictions
ggplot(pred_summary, aes(temp, .epred)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(
    x = "Temperature (°C)",
    y = "Estimated growth rate",
  ) +
  facet_wrap(~id2)

# create data frame as needed by JHuisman/conjugator

# D.t, R.t, T.t are the final population sizes
# D.0, R.0, T.0 are the initial population sizes. Initial population sizes are guestimates and set to 100,000, but importantly should be the same between donor and recipient.
# psi.D, psi.R, psi.T are the growth rates we have just estimated
# t is time of measurement
# id_cols are ID columns

#  create this from our existing data
d_conjugator <- select(
  conjugation_wide,
  strain,
  rep,
  date,
  temp,
  D.t = donor,
  R.t = recip,
  T.t = trans,
  t = time_h
) %>%
  mutate(D.0 = 100000, R.0 = 100000, T.0 = 0)

# create growth rate dataset

# recipient growth rates are the same for both plasmids
d_growth_r <- filter(pred_summary, id2 == "nal.5") %>%
  select(temp, .epred) %>%
  rename(psi.R = .epred)
# duplicate it and add plasmid column
d_growth_r <- bind_rows(
  d_growth_r %>% mutate(strain = "R1"),
  d_growth_r %>% mutate(strain = "RP4")
)

# donor and transconjugant growth rates are different
d_growth_dt <- select(pred_summary, temp, id2, .epred) %>%
  filter(id2 %in% c("wt + R1", "wt + RP4", "nal.5 + R1", 'nal.5 + RP4')) %>%
  mutate(
    strain = gsub(".*?\\+\\s*", "", id2) %>% trimws(),
    type = gsub("\\s*\\+.*", "", id2) %>% trimws()
  ) %>%
  select(-id2) %>%
  pivot_wider(names_from = type, values_from = .epred) %>%
  rename(psi.D = wt, psi.T = nal.5)

d_growth <- left_join(d_growth_dt, d_growth_r, by = c("temp", "strain"))

# join them together
d_conjugator <- left_join(d_conjugator, d_growth, by = c("temp", "strain"))

# Use JSHuisman/conjugator to estimate conjugation rates
huisman_conjrate <- estimate_conj_rate(
  d_conjugator,
  method = c("ASM", "T_DR"),
  id_cols = c("rep", "date", "temp", 'strain')
)
# a few negative rates, but we will make those 0 for now. They return as NAs

# replace NAs with 0
huisman_conjrate <- huisman_conjrate %>%
  mutate(estimate = ifelse(is.na(estimate), 0, estimate))

# plot ASM and T_DR estimates against each other.
# should correlate strongly, but not expected to be 1:1
pivot_wider(huisman_conjrate, names_from = method, values_from = estimate) %>%
  ggplot(aes(x = log10(T_DR), y = log10(ASM), color = rep)) +
  geom_point(size = 3, show.legend = FALSE) +
  theme_bw(base_size = 18) +
  labs(
    # make this t/dr delta time
    x = expression(paste(
      "log10 Conjugation rate (",
      t ~ "/" ~ dr ~ Delta ~ time ~ ')'
    )),
    y = "log Conjugation rate\n(Approximate Simonsen Method)"
  ) +
  scale_colour_paletteer_d("ggthemes::Classic_Color_Blind") +
  facet_wrap(~strain, labeller = labeller(strain = letter_facets)) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0)
  )
# ok these methods are very similar, use the T-DR method for the rest of the analysis.

# save this out
ggsave(
  'plots/conjugation_rate_comparison.png',
  width = 10,
  height = 5
)

#######################################
# fitting TPCs to conjugation data ####
#######################################

# the rates are super small, this could make fitting tricky
# check priors are working and they give a sensible curve
# rmaxtrans = exp(rmax) to force positive values. log(10^-11) = -25, meaning our prior needs to be around that. To get it on the same magnitude as the other parameters, we make rmaxtrans = exp(rmax * 100)

# do quick plot
ggplot(conjugation_wide, aes(temp, gamma)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~strain, scale = 'free_y')
# R1 looks like a sensible TPC, RP4 is not doing much of anything.

# we will analyse them separately

# save out processed dataset to be used for plotting in a separate script
select(conjugation_wide, strain, rep, temp, gamma) %>%
  write.csv('data/processed/conjugation_rates.csv', row.names = FALSE)

# write out the model
# set random effect of rep on every parameter
# except beta because we could not get a consistent fit
# fmt: skip
brms_form_random <- bf(
  gamma ~ inv_logit(25 * (temp - tmin)) *
    inv_logit(25 * (tmax - temp)) *
    rmax *
    exp(
      kappa *
        (inv_logit(alpha) *
          log(fmax(temp - tmin, 10^-20) / inv_logit(alpha)) +
          (1 - inv_logit(alpha)) *
            log(fmax(tmax - temp, 10^-20) / (1 - inv_logit(alpha))) -
          log(tmax - tmin))
    ),
  # this will put some of the parameters on similar scales
  # tmaxtrans is tmax/100
  brms::nlf(tmax ~ tmaxtrans * 100),
  brms::nlf(tmin ~ tmintrans * 100),
  # rmaxtrans is the exponential transformation of rmax, forcing it to be positive
  brms::nlf(rmax ~ exp(rmaxtrans * 100)),
  # alphatrans is alpha/10 to put it on the same scale as other parameters
  brms::nlf(alpha ~ alphatrans * 10),
  # set kappa
  brms::nlf(kappa ~ (inv_logit(alpha) * (1 - inv_logit(alpha)) / beta^2)),
  alphatrans ~ 1 + (1 | rep),
  beta ~ 1 + (1 | rep),
  tmintrans ~ 1 + (1 | rep),
  tmaxtrans ~ 1 + (1 | rep),
  rmaxtrans ~ 1 + (1 | rep),
  nl = TRUE
)

# set priors
# for setting priors, https://ben18785.shinyapps.io/distribution-zoo/ is a great website
# this helped me choose the priors, especially for beta
# probably want a wider alpha prior
priors <- c(
  prior(normal(0.12, 0.05), nlpar = "alphatrans"),
  prior(lognormal(-1.5, 0.3), nlpar = "beta", lb = 0),
  prior(normal(-0.25, 0.005), nlpar = "rmaxtrans"),
  prior(normal(0.15, 0.025), nlpar = "tmintrans", ub = 0.25),
  prior(normal(0.35, 0.05), nlpar = "tmaxtrans", lb = 0.25)
)

# set some priors on the random effects
random_priors <- c(
  set_prior("student_t(5,0,0.1)", class = "sd", nlpar = "alphatrans"),
  set_prior("student_t(10,0,0.025)", class = "sd", nlpar = "beta"),
  set_prior("student_t(10,0,0.005)", class = "sd", nlpar = "rmaxtrans"),
  set_prior("student_t(5,0,0.1)", class = "sd", nlpar = "tmaxtrans"),
  set_prior("student_t(5,0,0.1)", class = "sd", nlpar = "tmintrans"),
  set_prior("student_t(3,0,1e-12)", class = "sigma")
)
# combine priors
priors_random <- c(priors, random_priors)

table_priors <- priors_random %>%
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
save_as_image(
  table_priors,
  'tables/table_priors_plasmid.png',
  webshot = 'webshot2'
)
save_as_docx(table_priors, path = 'tables/table_priors_plasmid.docx')

# fit model

# check if there are any NAs
filter(conjugation_wide, is.na(gamma))

# sample just from the priors
model_check <- brm(
  brms_form_random,
  data = conjugation_wide,
  prior = priors_random,
  chains = 4,
  cores = 4,
  iter = 3000,
  init = 0,
  backend = 'cmdstan',
  sample_prior = 'only'
)

# plot priors
# grab draws, `b` usually determines fixed or population level effects
d_priors <- spread_draws(model_check, `b_.*`, regex = TRUE) %>%
  # remove b_ and _Intercept
  rename_with(~ gsub("b_|_Intercept", "", .x), starts_with("b_")) %>%
  # back transform parameters
  mutate(
    alpha = boot::inv.logit(alphatrans * 10),
    tmax = tmaxtrans * 100,
    rmax = exp(rmaxtrans * 100),
    tmin = tmintrans * 100
  ) %>%
  select(-c(alphatrans, tmaxtrans, rmaxtrans, tmintrans)) %>%
  pivot_longer(
    cols = c(alpha, beta, rmax, tmax, tmin),
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
  ) +
  theme(axis.text.x = element_text(size = 8))

ggsave(
  'plots/plasmid_tpc_priors.png',
  width = 8,
  height = 5
)

# plot priors for random effects
# grab draws, `sd_` usually determines random effects
d_priors_random <- spread_draws(model_check, `sd_.*`, regex = TRUE) %>%
  # remove b_ and _Intercept
  rename_with(~ gsub("sd_rep__|_Intercept", "", .x), starts_with("sd_")) %>%
  # back transform parameters
  mutate(
    tmax = tmaxtrans * 100,
    tmin = tmintrans * 100
  ) %>%
  select(-c(tmintrans, tmaxtrans)) %>%
  pivot_longer(
    cols = c(rmaxtrans, tmax, tmin, alphatrans, beta),
    names_to = 'parameter',
    values_to = 'value'
  )
# grab prior on sigma too
d_prior_sigma <- spread_draws(model_check, `sigma`, regex = TRUE) %>%
  rename(value = sigma) %>%
  mutate(parameter = 'sigma')

# combine and plot
d_priors_random <- bind_rows(d_priors_random, d_prior_sigma)

ggplot(d_priors_random, aes(value)) +
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
    strip.text = element_text(hjust = 0),
    axis.text.x = element_text(size = 8)
  )

ggsave(
  'plots/plasmid_tpc_priors_random.png',
  width = 8,
  height = 5
)

controls <- list(adapt_delta = 0.99, max_treedepth = 12)

# fit model
model_r1 <- brm(
  brms_form_random,
  data = filter(conjugation_wide, strain == "R1"),
  prior = priors_random,
  chains = 4,
  cores = 4,
  iter = 5000,
  init = 0,
  backend = 'cmdstan',
  #control = controls
)

summary(model_r1)
plot(model_r1)

# model for RP4
model_rp4 <- brm(
  brms_form_random,
  data = filter(conjugation_wide, strain == "RP4"),
  prior = priors_random,
  chains = 4,
  cores = 4,
  iter = 5000,
  init = 0,
  backend = 'cmdstan',
  #control = controls
)

# quick model checks
summary(model_rp4)
plot(model_rp4)

plot(conditional_effects(model_rp4), points = TRUE, re_formula = NULL)

# save out models
saveRDS(model_r1, 'data/models/brms_flextpc_r1.rds')
saveRDS(model_rp4, 'data/models/brms_flextpc_rp4.rds')

#------------------------------------#
# do posterior probability checks ####
#------------------------------------#

# for fixed effects
# look at how much the model has changed based on the priors
d_prior_post <- spread_draws(model_r1, `b_.*`, regex = TRUE) %>%
  # remove b_ and _Intercept
  rename_with(~ gsub("b_|_Intercept", "", .x), starts_with("b_")) %>%
  # back transform parameters
  mutate(
    alpha = boot::inv.logit(alphatrans * 10),
    tmax = tmaxtrans * 100,
    rmax = exp(rmaxtrans * 100),
    tmin = tmintrans * 100,
    type = 'posterior'
  ) %>%
  pivot_longer(
    cols = c(alpha, beta, rmax, tmax, tmin),
    names_to = 'parameter',
    values_to = 'value'
  ) %>%
  bind_rows(
    d_priors %>%
      mutate(type = 'prior')
  )

ggplot(
  d_prior_post,
  aes(value, col = type, group = interaction(.chain, type))
) +
  geom_density(fill = NA, key_glyph = draw_key_rect) +
  facet_wrap(
    ~parameter,
    scales = 'free',
    labeller = labeller(parameter = letter_facets)
  ) +
  theme_bw(base_size = 18) +
  scale_color_manual(name = '', values = c('#FF000080', '#00000080')) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = 'Parameter value', y = 'Density') +
  theme(
    legend.position = c(0.85, 0.25),
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0),
    axis.text.x = element_text(size = 8)
  )

# save out this plot
ggsave(
  'plots/r1_tpc_priors_posteriors.png',
  width = 8,
  height = 5
)

# for the random effects
# do these need to backtransformed?
d_post_pop <- spread_draws(model_r1, `sd_.*|sigma`, regex = TRUE) %>%
  # remove b_ and _Intercept
  rename_with(~ gsub("sd_rep__|_Intercept", "", .x), starts_with("sd_")) %>%
  # back transform parameters
  mutate(
    tmax = tmaxtrans * 100,
    tmin = tmintrans * 100,
    type = 'posterior'
  ) %>%
  select(-c(tmintrans, tmaxtrans)) %>%
  pivot_longer(
    cols = c(alphatrans, tmin, tmax, rmaxtrans, beta, sigma),
    names_to = 'parameter',
    values_to = 'value'
  ) %>%
  mutate(type = 'posterior') %>%
  bind_rows(
    d_priors_random %>%
      mutate(type = 'prior')
  )

ggplot(
  d_post_pop,
  aes(value, col = type, group = interaction(.chain, type))
) +
  geom_density(fill = NA, key_glyph = draw_key_rect) +
  facet_wrap(
    ~parameter,
    scales = 'free',
    labeller = labeller(parameter = letter_facets)
  ) +
  theme_bw(base_size = 18) +
  scale_color_manual(name = '', values = c('#FF000080', '#00000080')) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = 'Parameter value', y = 'Density') +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0),
    legend.position = 'bottom',
    axis.text = element_text(size = 7)
  )

ggsave(
  'plots/r1_tpc_priors_posteriors_random.png',
  width = 8,
  height = 5
)

# do the same with rp4
# for fixed effects
# look at how much the model has changed based on the priors
d_prior_post <- spread_draws(model_rp4, `b_.*`, regex = TRUE) %>%
  # remove b_ and _Intercept
  rename_with(~ gsub("b_|_Intercept", "", .x), starts_with("b_")) %>%
  # back transform parameters
  mutate(
    alpha = boot::inv.logit(alphatrans * 10),
    tmax = tmaxtrans * 100,
    rmax = exp(rmaxtrans * 100),
    tmin = tmintrans * 100,
    type = 'posterior'
  ) %>%
  pivot_longer(
    cols = c(alpha, beta, rmax, tmax, tmin),
    names_to = 'parameter',
    values_to = 'value'
  ) %>%
  bind_rows(
    d_priors %>%
      mutate(type = 'prior')
  )

ggplot(
  d_prior_post,
  aes(value, col = type, group = interaction(.chain, type))
) +
  geom_density(fill = NA, key_glyph = draw_key_rect) +
  facet_wrap(
    ~parameter,
    scales = 'free',
    labeller = labeller(parameter = letter_facets)
  ) +
  theme_bw(base_size = 18) +
  scale_color_manual(name = '', values = c('#FF000080', '#00000080')) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = 'Parameter value', y = 'Density') +
  theme(
    legend.position = c(0.85, 0.25),
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0),
    axis.text.x = element_text(size = 8)
  )

# save out this plot
ggsave(
  'plots/rp4_tpc_priors_posteriors.png',
  width = 8,
  height = 5
)

# for the random effects

# do it for model_rp4
d_post_pop <- spread_draws(model_rp4, `sd_.*|sigma`, regex = TRUE) %>%
  # remove b_ and _Intercept
  rename_with(~ gsub("sd_rep__|_Intercept", "", .x), starts_with("sd_")) %>%
  # back transform parameters
  mutate(
    tmax = tmaxtrans * 100,
    tmin = tmintrans * 100,
    type = 'posterior'
  ) %>%
  select(-c(tmintrans, tmaxtrans)) %>%
  pivot_longer(
    cols = c(alphatrans, tmin, tmax, rmaxtrans, beta, sigma),
    names_to = 'parameter',
    values_to = 'value'
  ) %>%
  mutate(type = 'posterior') %>%
  bind_rows(d_priors_random %>% mutate(type = 'prior'))

ggplot(
  d_post_pop,
  aes(value, col = type, group = interaction(.chain, type))
) +
  geom_density(fill = NA, key_glyph = draw_key_rect) +
  facet_wrap(
    ~parameter,
    scales = 'free',
    labeller = labeller(parameter = letter_facets)
  ) +
  theme_bw(base_size = 18) +
  scale_color_manual(name = '', values = c('#FF000080', '#00000080')) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = 'Parameter value', y = 'Density') +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0),
    legend.position = 'bottom',
    axis.text = element_text(size = 7)
  )

ggsave(
  'plots/rp4_tpc_priors_posteriors_random.png',
  width = 8,
  height = 5
)
