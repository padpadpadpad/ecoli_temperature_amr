#----------------------------------
# Purpose of the script: fitting Thermal Performance Curves (TPCs) to wild type, nalR, and plasmid-bearing strains of E. coli
#
# What this script does:
# 1. Reads in growth rate datasets
# 2. Wrangles and combines datasets
# 3. Defines model
# 4. Sets priors and samples to check priors
# 5. Fits model to test data
# 6. Fits full model and check model fit
# 7. Saves out model and data for plotting in create_figure_1.R
#
# Author: Dr. Daniel Padfield
# Date Created: 03 December 2025
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
  tidyverse,
  brms,
  tidybayes,
  bayesplot,
  patchwork,
  job,
  flextable,
  stan - dev / cmdstanr
)
## ---------------------------

# set seed for reproducibility
set.seed(42)

# load custom functions
source('scripts/00_functions.R')

# set fmax to pmax
# this is needed for some of the functions in brms and other Bayesian R packages to work
fmax <- pmax

# read in tpc data across antibiotics
d_metrics <- read.csv(
  'data/processed/05_kanamycin_metrics_analysis.csv'
)

# filter to not have any antibiotics
d_metrics <- filter(d_metrics, conc_µl == 0)

# select key columns
# convert id2 to show which strain has the plasmid
d_metrics <- select(d_metrics, id2, temp, gr) %>%
  mutate(
    id2 = case_when(
      id2 == 'R1' ~ 'wt + R1',
      id2 == 'RP4' ~ 'wt + RP4',
      TRUE ~ id2
    )
  )

unique(d_metrics$id2)

# read in tpc data for nal.5 with R1 and RP4 in the absence of antibiotics
# these are modelled simply for estimating plasmid transfer rates
d_metrics_nal <- read.csv(
  'data/processed/04_Trans_gr_metrics.csv'
)

# rename nalR mutants when they carry plasmids
d_metrics_nal <- select(d_metrics_nal, id2 = id, temp, gr) %>%
  filter(id2 %in% c('R1', 'RP4')) %>%
  mutate(
    id2 = case_when(
      id2 == 'R1' ~ 'nal.5 + R1',
      id2 == 'RP4' ~ 'nal.5 + RP4',
      TRUE ~ id2
    )
  )

unique(d_metrics_nal$id2)

# bind the two datasets together
d_metrics <- bind_rows(d_metrics, d_metrics_nal)

nrow(d_metrics)

# how many points do we have for each curve?
# we have a lot of replicates for the nal5 strains containing the plasmid
d_metrics %>%
  group_by(id2) %>%
  summarise(n = n())

# save out dataset for plotting in figure 1
write.csv(
  d_metrics,
  'data/processed/temperature_growth_data.csv',
  row.names = FALSE
)

# quick plot of the data
ggplot(d_metrics, aes(temp, gr, col = id2)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~id2)

#--------------------------------------#
# fitting flextpc on a single curve ####
#--------------------------------------#

# write out brms model
# we can then fit this model using transformed parameters
# inv_logit allows a smooth transition between when the rate is either 0 or positive. This helps with convergence rather than having a hard cutoff between 0 and 1. This does not make a difference to the parameters of the model, as can be seen in 00_flextpc_check.R
# alpha is inv_logit to ensure it is between 0 and 1. alphatrans is alpha/10 to put it on the same scale as other parameters
# rmax is exp(rmaxtrans) to ensure it is positive
# beta is on the original scale, but we set a lognormal prior on it to ensure it is positive
# tmaxtrans and tmintrans are tmax/100 and tmin/100 to put them on the same scale as other parameters
# fmt: skip
brms_form <- bf(
  gr ~ inv_logit(25 * (temp - tmin)) *
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
  brms::nlf(rmax ~ exp(rmaxtrans)),
  # alphatrans is alpha/10 to put it on the same scale as other parameters
  brms::nlf(alpha ~ alphatrans * 10),
  # set kappa
  brms::nlf(kappa ~ (inv_logit(alpha) * (1 - inv_logit(alpha)) / beta^2)),
  alphatrans + beta + tmintrans + tmaxtrans + rmaxtrans ~ 1,
  nl = TRUE
)

# write function so that we can plot a set of parameters
# this function takes values on the original scale
# fmt: skip
flextpc_2024_mod <- function(temp, tmin, tmax, rmax, alpha, beta) {
  
  # transform parameters back
  tmax <- tmax / 100
  rmax <- log(rmax)
  alpha <- alpha / 10
  tmin <- tmin / 100
  
  # calculate kappa
    kappa <- (boot::inv.logit(alpha) * (1 - boot::inv.logit(alpha)) / beta^2)

    out <- 
    exp(rmax) *
    exp(kappa * (boot::inv.logit(alpha * 10) *
        log(pmax(temp - (tmin * 100), 10^-20) / boot::inv.logit(alpha * 10)) + (1 - boot::inv.logit(alpha * 10)) *
        log(pmax((tmax * 100) - temp, 10^-20) / (1 - boot::inv.logit(alpha * 10))) -
        log((tmax * 100) - (tmin * 100)))
    )
  
  return(out)

}

# set priors
# for setting priors, https://ben18785.shinyapps.io/distribution-zoo/ is a great website
# this helped me choose the priors, especially for beta
priors <- c(
  prior(normal(0.14, 0.04), nlpar = "alphatrans"),
  prior(lognormal(-1.5, 0.3), nlpar = "beta", lb = 0),
  prior(normal(0.6, 0.25), nlpar = "rmaxtrans"),
  prior(normal(0.0, 0.025), nlpar = "tmintrans"),
  prior(normal(0.35, 0.05), nlpar = "tmaxtrans", lb = 0.2),
  set_prior("student_t(3,0,0.2)", class = "sigma")
)

# make a table of the custom priors
table_priors <- priors %>%
  data.frame() %>%
  mutate(nlpar = ifelse(nlpar == '', class, nlpar)) %>%
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

table_priors

# save this out
save_as_image(table_priors, 'tables/table_priors_tpc.png', webshot = 'webshot2')
save_as_docx(table_priors, path = 'tables/table_priors_tpc.docx')

# we can see what these priors are like as a fit
d_check <- tibble(temp = 1:50) %>%
  mutate(
    rate = flextpc_2024_mod(
      temp,
      tmin = 0,
      tmax = 35,
      rmax = exp(0.6),
      alpha = 1.4,
      beta = rlnorm(1, -1.5, 0)
    )
  )

# plot
ggplot(d_check, aes(x = temp, y = rate)) +
  geom_line() +
  theme_bw()

# Looks great

# sample only from the priors first
# on a test data set
test <- filter(d_metrics, id2 == 'wt')

model_check <- brm(
  brms_form,
  data = test,
  prior = priors,
  chains = 4,
  cores = 1,
  iter = 3000,
  sample_prior = 'only',
  init = '0',
  backend = 'cmdstan'
)

# grab draws, `b` usually determines fixed or population level effects
d_priors <- spread_draws(model_check, `b_.*`, regex = TRUE) %>%
  # remove b_ and _Intercept
  rename_with(~ gsub("b_|_Intercept", "", .x), starts_with("b_")) %>%
  # back transform parameters
  mutate(
    alpha = boot::inv.logit(alphatrans * 10),
    tmax = tmaxtrans * 100,
    rmax = exp(rmaxtrans),
    tmin = tmintrans * 100
  ) %>%
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
  )

ggsave(
  'plots/growth_tpc_priors.png',
  width = 8,
  height = 5
)

# seems all reasonable
# fit the model on the test dataset
model1 <- brm(
  brms_form,
  data = test,
  prior = priors,
  chains = 4,
  cores = 2,
  iter = 3000,
  init = '0',
  backend = 'cmdstan'
)

# we can plot the model
plot(conditional_effects(model1), points = TRUE)[[1]]

# we can look at some model checks
pp_check(model1, ndraws = 100)

# pairs plot
pairs(model1)
# beta and tmin are quite correlated, as are alpha and tmax

# parameter distributions
plot(model1)

# look at posterior vs prior distribution
# check available parameters
variables(model1)

# look at how much the model has changed based on the priors
d_prior_post <- spread_draws(model1, `b_.*`, regex = TRUE) %>%
  # remove b_ and _Intercept
  rename_with(~ gsub("b_|_Intercept", "", .x), starts_with("b_")) %>%
  # back transform parameters
  mutate(
    alpha = boot::inv.logit(alphatrans * 10),
    tmax = tmaxtrans * 100,
    rmax = exp(rmaxtrans),
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

# compare prior to posterior
ggplot(
  d_prior_post,
  aes(value, col = type, group = interaction(.chain, type))
) +
  geom_density(fill = NA, key_glyph = draw_key_rect) +
  facet_wrap(~parameter, scales = 'free') +
  theme_bw() +
  scale_color_manual(name = '', values = c('red', 'black'))
# this is fun, you can see how your parameters change from prior to posterior
# all looks good, prior on tmin has not changed much, but thats because we do not have much information on them

#---------------#
# full model ####
#---------------#

# scaling up to all the curves without antibiotics
# we fit a model with a fixed effect of id2 (strain) on all parameters
# 0 + id2 means that we fit a separate parameter for each id2

# write the full model
# fmt: skip
brms_form_v2 <- bf(
  gr ~ inv_logit(25 * (temp - tmin)) *
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
  brms::nlf(rmax ~ exp(rmaxtrans)),
  # alphatrans is alpha/10 to put it on the same scale as other parameters
  brms::nlf(alpha ~ alphatrans * 10),
  # set kappa
  brms::nlf(kappa ~ (inv_logit(alpha) * (1 - inv_logit(alpha)) / beta^2)),
  alphatrans ~ 0 + id2,
  beta ~ 0 + id2,
  tmintrans ~ 0 + id2,
  tmaxtrans ~ 0 + id2,
  rmaxtrans ~ 0 + id2,
  nl = TRUE
)

# set controls
controls <- list(adapt_delta = 0.99, step_size = 0.001, max_treedepth = 10)

# if you do not want to run the model, can load it in
model <- readRDS('data/models/brms_flextpc.rds')

# run the full model
model <- brm(
  brms_form_v2,
  data = d_metrics,
  prior = priors,
  chains = 4,
  cores = 2,
  iter = 3000,
  control = controls,
  init = 0,
  backend = 'cmdstan'
)

saveRDS(model, 'data/models/brms_flextpc.rds')
# all rhats are similar
summary(model)

plot(model)

# plot predictions of the model
conditions <- data.frame(id2 = unique(d_metrics$id2))
rownames(conditions) <- unique(d_metrics$id2)

preds_bact_tpc <- conditional_effects(
  model,
  conditions = conditions,
  method = "posterior_epred"
)

# this looks pretty good
plot(preds_bact_tpc, ncol = 2, points = TRUE, plot = FALSE)[[1]]

# plot predictions
preds_bact_tpc <- select(d_metrics, id2) %>%
  distinct() %>%
  group_by(id2) %>%
  do(data.frame(
    temp = seq(min(d_metrics$temp), max(d_metrics$temp), length.out = 200)
  )) %>%
  add_epred_draws(
    model,
    ndraws = 1000
  ) %>%
  group_by(id2, temp) %>%
  median_hdci()

# for plotting, we only want TPCs for wild-type and wild-type + R1.
d_plot <- filter(d_metrics, id2 %in% c('wt', 'wt + R1', 'wt + RP4'))
d_preds_bact_tpc <- filter(
  preds_bact_tpc,
  id2 %in% c('wt', 'wt + R1', 'wt + RP4')
)

# plot predictions
ggplot() +
  geom_point(
    aes(temp, gr, col = id2),
    data = d_plot
  ) +
  geom_ribbon(
    aes(
      temp,
      ymin = .lower,
      ymax = .upper,
      fill = id2
    ),
    data = d_preds_bact_tpc,
    alpha = 0.3,
    show.legend = FALSE
  ) +
  geom_line(
    aes(temp, .epred, col = id2),
    data = d_preds_bact_tpc,
    show.legend = FALSE
  ) +
  facet_wrap(~id2) +
  theme_bw(base_size = 14) +
  labs(x = 'Temperature (ºC)', y = 'Growth rate')
# this is good, but the final plot is made in figure_1.R

# check posterior vs prior distributions for full model
# as we used 0 + id2, there is no intercept

# look at how much the model has changed based on the priors
d_prior_post <- gather_draws(model, `b_.*`, regex = TRUE) %>%
  # remove b_ and _Intercept
  mutate(.variable = gsub("b_|_intercept", "", .variable)) %>%
  separate(.variable, into = c('parameter', 'id2'), sep = '_') %>%
  pivot_wider(names_from = parameter, values_from = .value) %>%
  # back transform parameters
  mutate(
    alpha = boot::inv.logit(alphatrans * 10),
    tmax = tmaxtrans * 100,
    rmax = exp(rmaxtrans),
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
      mutate(type = 'prior', id2 = 'prior')
  )

unique(d_prior_post$id2)

# compare prior to posterior
ggplot(
  d_prior_post,
  aes(value, group = interaction(.chain, type, id2))
) +
  geom_density(aes(color = type), fill = NA, key_glyph = draw_key_rect) +
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
    strip.text = element_text(hjust = 0)
  )

# save out this plot
ggsave(
  'plots/growth_tpc_priors_posteriors.png',
  width = 8,
  height = 5
)
