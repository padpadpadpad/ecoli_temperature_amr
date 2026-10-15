#----------------------------------
# Purpose of the script: Fit beta model to invasion data.
#
# What this script does:
# 1. Reads in invasion data
# 2. Fits beta regression model to proportion of resistant bacteria
# 3. Fits negative binomial to count data
#
# Author: Dr. Daniel Padfield
# Date Created: 25 April 2026
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
librarian::shelf(
  readxl,
  brms,
  marginaleffects,
  tidybayes,
  patchwork,
  tidyverse,
  janitor,
  performance,
  paletteer
)
## ---------------------------

# load in functions
source('scripts/00_functions.R')

# custom function
# letter facets but allow starting at a different letter.
# and add a ºC to the end
letter_facets_mod <- function(string, start_letter) {
  # find number to start at based on start_letter
  start_num <- which(letters == start_letter)
  len <- length(string)
  string <- paste(
    "(",
    letters[start_num:(start_num + len - 1)],
    ") ",
    string,
    'ºC',
    sep = ""
  )
  return(string)
}

# this link is great: https://www.andrewheiss.com/blog/2021/11/08/beta-regression-guide/#explore-the-data

# set seed so we always get the same model and output for calculating credible intervals
set.seed(79)

# read in data
invasion <- read_xlsx('data/raw/invasion/invasion.xlsx', sheet = 1)

# change some of the columns to the correct format
invasion <- invasion %>%
  mutate(
    rep = as.character(rep),
    conc = as.character(conc),
    temp = as.character(temp),
    density = as.numeric(density)
  )

# create dataframe with numbers of resistant, susceptible, and transconjugants
# do some processing
invasion_2 <- invasion %>%
  dplyr::select(strain, plate, conc, temp, dil, density, rep) %>%
  # separate datasets for each selective plate
  group_by(strain, temp, conc, rep, plate) %>%
  mutate(n = 1:n()) %>%
  ungroup() %>%
  pivot_wider(
    id_cols = c(strain, temp, conc, rep, n),
    names_from = plate,
    values_from = density
  ) %>%
  mutate(
    D = D - T,
    R = R - T,
    total_res = T + D,
    total_d = R + D + T,
    prop_T = T / total_d,
    prop_D = D / total_d,
    prop_R = R / total_d,
    prop_resistant = (T + D) / total_d,
    log_total = log10(total_d),
    log_res = log10(total_res)
  )

# check NAs
filter(invasion_2, is.na(prop_resistant))
# these are NA because...

# get rid of NAs
invasion_2 <- invasion_2 %>%
  filter(!is.na(total_d))

# quick plot
ggplot(invasion_2, aes(x = conc, y = prop_resistant, color = conc)) +
  geom_point() +
  geom_boxplot() +
  facet_wrap(
    strain ~ temp,
    scales = "free_y",
    labeller = labeller(.multi_line = FALSE)
  )

ggplot(invasion_2, aes(x = conc, y = total_d, color = conc)) +
  geom_point() +
  geom_boxplot() +
  facet_wrap(
    rep ~ temp,
    scales = "free_y",
    labeller = labeller(.multi_line = FALSE)
  )
# missing 21ºC, 15µg/mL for one rep of RP4, everything else is there

#-------------------------#
# statistical analysis ####
#-------------------------#

# we can do a beta regression for proportion of resistant bacteria using brms
# if there are any 0s or 1s, make them very close to 0 or 1
# modelling proportion change takes into account the impact of temperature and antibiotics on both susceptible and resistant bacteria.

summary(invasion_2$prop_resistant)
# minimum prop resistant is 0.008, when we put in ~0.1
# so there is a cost to being resistant

# create a column for conc2 that has ug/mL on the end of conc and orders it by size of conc
invasion_2 <- mutate(invasion_2, conc2 = paste(conc, "µg/mL", sep = "")) %>%
  mutate(conc2 = fct_reorder(conc2, as.numeric(conc)))

# lets make a histogram to see proportion of resistance across kanamycin concentration

panel_labeller <- function(labels) {
  combined <- do.call(paste, c(labels, sep = ", "))
  list(letter_facets(combined))
}

ggplot(invasion_2, aes(x = prop_resistant)) +
  geom_histogram(bins = 30, fill = 'grey', color = 'black') +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "red") +
  theme_bw() +
  facet_wrap(strain ~ conc2, labeller = panel_labeller, ncol = 5) +
  labs(x = "Proportion of resistant bacteria", y = "Count") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0)
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = c("0", "0.25", "0.5", "0.75", "1")
  )
# kind of makes sense. As you increase the antibiotic concentration, there is a higher proportion of resistant bacteria.

# save this out
ggsave(
  'plots/invasion_histogram.png',
  width = 8,
  height = 4
)

# check whether plasmid transfer might be costly. Do we see lower proportions of resistant bacteria where conjugation is happening and selection isnt very strong. Signal isn't great but some suggestion this is true. Could we look at this from the model output?
group_by(invasion_2, strain, conc) %>%
  summarise(
    mean_prop_resistant = mean(prop_resistant),
    median_prop_resistant = median(prop_resistant),
    .groups = "drop"
  )

filter(invasion_2, strain == 'R1') %>%
  group_by(conc, temp) %>%
  summarise(
    mean_prop_resistant = mean(prop_resistant),
    median_prop_resistant = median(prop_resistant),
    .groups = "drop"
  )

# set up some colours for plotting
my_colours <- c('light grey', paletteer::paletteer_c("ggthemes::Red", 4))
names(my_colours) <- c("0", "2", "3.5", "5", "15")

# set some priors, recommended by Andrew Heiss in his post at the top of the script
priors <- c(
  set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
  set_prior("normal(0, 2.5)", class = "b")
)

table_invasion_priors <- priors %>%
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
  table_invasion_priors,
  'tables/table_priors_invasion.png',
  webshot = 'webshot2'
)
save_as_docx(table_invasion_priors, path = 'tables/table_priors_invasion.docx')

# run model
# has a parameter for the mean, and a parameter for the precision (phi)
mod_beta <- brm(
  bf(
    prop_resistant ~ strain * conc * temp + (1 | rep),
    phi ~ strain * conc * temp + (1 | rep)
  ),
  data = invasion_2,
  iter = 10000,
  family = Beta(),
  prior = priors,
  backend = 'cmdstanr',
  cores = 4
)

# save out model
saveRDS(mod_beta, 'data/models/beta_invasion.rds')
mod_beta <- readRDS('data/models/beta_invasion.rds')
# check the fit of the model to the data
pp_check(mod_beta, ndraws = 100)
# really good

# look at the model summary
summary(mod_beta)

performance::r2_bayes(mod_beta)

# lets get some predictions

# we can look at differences between concentrations within each temperature
# using marginal effects
mod_contrasts <- mod_beta %>%
  avg_comparisons(
    variables = list(conc = 'reference'),
    by = c("strain", "temp"),
    type = "response",
    re_formula = NULL
  )

mod_contrasts

# we can grab these out as MCMC draws using marginaleffects too
mod_contrasts_draws <- mod_beta %>%
  avg_comparisons(
    variables = list(conc = 'reference'),
    by = c("strain", "temp"),
    type = "response",
    re_formula = NULL
  ) %>%
  posterior_draws() %>%
  janitor::clean_names() %>%
  mutate(
    conc = as.character(parse_number(contrast)),
    contrast2 = paste("0 - ", conc, sep = "")
  )

# estimate per rep random effects in the contrasts
mod_contrasts_draws_re <- mod_beta %>%
  avg_comparisons(
    variables = list(conc = 'reference'),
    by = c("strain", "temp", 'rep'),
    type = "response",
    re_formula = NULL
  ) %>%
  posterior_draws() %>%
  janitor::clean_names() %>%
  mutate(
    conc = as.character(parse_number(contrast)),
    contrast2 = paste("0 - ", conc, sep = "")
  )

# can use marginal effects to look at key differences and estimate means
# COMMENT do not want to predict over levels we do not have, so create new dataset
preds_beta <- mod_beta %>%
  epred_draws(
    newdata = select(invasion_2, strain, conc, temp, rep) %>% distinct(),
    re_formula = NULL
  ) %>%
  janitor::clean_names() %>%
  group_by(strain, conc, temp, draw) %>%
  summarise(epred = mean(epred), .groups = "drop")

preds_re_beta <- mod_beta %>%
  epred_draws(
    newdata = select(invasion_2, strain, conc, temp, rep) %>% distinct(),
    re_formula = NULL
  ) %>%
  janitor::clean_names()

preds_beta %>%
  group_by(strain, conc, temp) %>%
  median_qi(epred)

# use invasion data
d_trans <- filter(invasion, plate == "T") %>%
  mutate(
    dil = as.numeric(dil),
    count = as.numeric(count),
    count = replace_na(count, 0)
  )

# look at the max and min dilution
d_trans$dil %>% max()
d_trans$dil %>% min()
# smallest dilution is 10^0, so normalise to 10^0

# normalise counts to 10^-0 to get consistent counts across dilutions
d_trans <- mutate(d_trans, count2 = count * 10^(dil - 0))

# changed this from count2 to count
ggplot(d_trans, aes(x = conc, y = count2, color = conc)) +
  geom_point() +
  geom_boxplot() +
  facet_wrap(
    strain ~ temp,
    scales = "free_y",
    labeller = labeller(.multi_line = FALSE)
  )

filter(d_trans, count2 == 0)

# fit a negative binomial model
# have a parameter on the dispersion parameter as well and fit a negative binomial
# what are the effects for the dispersion parameter?
# will try:
# 1. the same as the predictors for mean
# 2. just replicate representing day (e.g. the overnight)
# 3. same as predictors for mean but no interactions

# read in models
mod_negbinom <- readRDS('data/models/negbinom_invasion.rds')
mod_negbinom2 <- readRDS('data/models/negbinom_invasion2.rds')
mod_negbinom3 <- readRDS('data/models/negbinom_invasion3.rds')

mod_negbinom <- brm(
  bf(
    count2 ~ strain * conc * temp + (1 | rep),
    shape ~ strain * conc * temp + (1 | rep)
  ),
  data = d_trans,
  iter = 5000,
  family = negbinomial(),
  backend = 'cmdstanr',
  cores = 4,
  save_pars = save_pars(all = TRUE)
)

mod_negbinom <- add_criterion(
  mod_negbinom,
  criterion = "loo",
  save_psis = TRUE,
  reloo = TRUE,
  cores = 4,
  backend = "cmdstanr",
  silent = 1
)

mod_negbinom2 <- brm(
  bf(
    count2 ~ strain * conc * temp + (1 | rep),
    shape ~ (1 | rep)
  ),
  data = d_trans,
  iter = 5000,
  family = negbinomial(),
  backend = 'cmdstanr',
  cores = 4,
  save_pars = save_pars(all = TRUE)
)

mod_negbinom2 <- add_criterion(
  mod_negbinom2,
  criterion = "loo",
  save_psis = TRUE,
  reloo = TRUE,
  cores = 4,
  backend = "cmdstanr"
)

mod_negbinom3 <- brm(
  bf(
    count2 ~ strain + conc + temp + (1 | rep),
    shape ~ strain + conc + temp + (1 | rep)
  ),
  data = d_trans,
  iter = 5000,
  family = negbinomial(),
  backend = 'cmdstanr',
  cores = 4,
  save_pars = save_pars(all = TRUE)
)

mod_negbinom3 <- add_criterion(
  mod_negbinom3,
  criterion = "loo",
  save_psis = TRUE,
  reloo = TRUE,
  cores = 4,
  backend = "cmdstanr",
  overwrite = TRUE
)

mod_negbinom3$criteria$loo
mod_negbinom2$criteria$loo
mod_negbinom$criteria$loo

loo_compare(mod_negbinom, mod_negbinom2, mod_negbinom3)
# so the most complicated shape is the best, even though its ppcheck sometimes predicts some very high values

pp_check(mod_negbinom, ndraws = 100, type = 'ecdf_overlay')

summary(mod_negbinom)
# everything has converged

# plot number of transconjugants
preds_negbinom <- mod_negbinom %>%
  epred_draws(
    newdata = select(invasion_2, strain, conc, temp, rep) %>% distinct(),
    re_formula = NULL
  ) %>%
  janitor::clean_names() %>%
  group_by(strain, conc, temp, draw) %>%
  summarise(epred = mean(epred), .groups = "drop")

preds_re_negbinom <- mod_negbinom %>%
  epred_draws(
    newdata = select(invasion_2, strain, conc, temp, rep) %>% distinct(),
    re_formula = NULL
  ) %>%
  janitor::clean_names()

# panel 1: proportion resistant
make_p1 <- function(strain_name, start_letter = "a") {
  ggplot(
    preds_beta %>% filter(strain == strain_name),
    aes(x = fct_reorder(conc, as.numeric(conc)), y = epred, fill = conc)
  ) +
    stat_eye(
      aes(fill = conc),
      .width = c(0.8, 0.95),
      point_interval = "median_hdci",
      point_size = 4,
      shape = 21
    ) +
    stat_pointinterval(
      aes(x = fct_reorder(conc, as.numeric(conc)), y = epred, group = rep),
      position = position_dodge(width = 0.5),
      point_size = 1,
      .width = c(0.8, 0.95),
      point_interval = "median_hdci",
      data = preds_re_beta %>% filter(strain == strain_name),
      alpha = 0.5
    ) +
    geom_dots(
      aes(x = fct_reorder(conc, as.numeric(conc)), y = prop_resistant),
      data = invasion_2 %>% filter(strain == strain_name),
      side = 'both',
      layout = 'hex',
      dotsize = 1,
      col = 'black',
      stackratio = 0.7,
      binwidth = 0.02
    ) +
    facet_wrap(
      ~temp,
      labeller = labeller(temp = function(labels) {
        letter_facets_mod(labels, start_letter = start_letter)
      })
    ) +
    labs(
      x = "Kanamycin concentration (µg/mL)",
      y = "Proportion of resistant bacteria"
    ) +
    theme_bw(base_size = 14) +
    scale_fill_manual(values = my_colours) +
    ylim(c(0, 1)) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(hjust = 0)
    )
}

# panel 2: contrasts
make_p2 <- function(strain_name, start_letter = "d") {
  ggplot(
    mod_contrasts_draws %>% filter(strain == strain_name),
    aes(
      x = forcats::fct_reorder(contrast2, as.numeric(conc)),
      y = draw,
      fill = conc
    )
  ) +
    stat_eye(
      .width = c(0.8, 0.95),
      point_interval = "median_hdci",
      point_size = 3
    ) +
    stat_pointinterval(
      aes(
        x = forcats::fct_reorder(contrast2, as.numeric(conc)),
        y = draw,
        group = rep
      ),
      position = position_dodge(width = 0.5),
      point_size = 1,
      .width = c(0.8, 0.95),
      point_interval = "median_hdci",
      data = mod_contrasts_draws_re %>% filter(strain == strain_name),
      alpha = 0.5
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(
      ~temp,
      labeller = labeller(temp = function(labels) {
        letter_facets_mod(labels, start_letter = start_letter)
      })
    ) +
    labs(x = "Contrast", y = "Difference in\nresistance proportion") +
    theme_bw(base_size = 14) +
    scale_fill_manual(values = my_colours) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(hjust = 0)
    )
}

# panel 3: transconjugants
make_p3 <- function(strain_name, start_letter = "g") {
  ggplot(
    preds_negbinom %>% filter(strain == strain_name),
    aes(x = fct_reorder(conc, as.numeric(conc)), y = epred, fill = conc)
  ) +
    stat_eye(
      aes(fill = conc),
      .width = c(0.8, 0.95),
      point_interval = "median_hdci",
      point_size = 4,
      shape = 21
    ) +
    stat_pointinterval(
      aes(x = fct_reorder(conc, as.numeric(conc)), y = epred, group = rep),
      position = position_dodge(width = 0.5),
      point_size = 1,
      .width = c(0.8, 0.95),
      point_interval = "median_hdci",
      data = preds_re_negbinom %>% filter(strain == strain_name),
      alpha = 0.5
    ) +
    geom_dots(
      aes(x = fct_reorder(conc, as.numeric(conc)), y = count2),
      data = d_trans %>% filter(strain == strain_name),
      side = 'both',
      layout = 'hex',
      dotsize = 10,
      col = 'black',
      stackratio = 0.7,
      binwidth = 125
    ) +
    facet_wrap(
      ~temp,
      labeller = labeller(temp = function(labels) {
        letter_facets_mod(labels, start_letter = start_letter)
      })
    ) +
    labs(
      x = "Kanamycin concentration (µg/mL)",
      y = "Number of transconjugants (per 10µL)"
    ) +
    ylim(0, 50000) +
    theme_bw(base_size = 14) +
    scale_fill_manual(values = my_colours) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(hjust = 0)
    )
}
# R1
p1_r1 <- make_p1('R1', 'a')
p2_r1 <- make_p2('R1', 'd')
p3_r1 <- make_p3('R1', 'g')
p1_r1 / p2_r1 / p3_r1 + plot_layout(heights = c(0.33, 0.33, 0.33))

#figure 4
ggsave(
  'plots/figure_4.png',
  width = 9,
  height = 12
)

# RP4
p1_rp4 <- make_p1('RP4', 'a')
p2_rp4 <- make_p2('RP4', 'd')
p3_rp4 <- make_p3('RP4', 'g')
p1_rp4 / p2_rp4 / p3_rp4 + plot_layout(heights = c(0.33, 0.33, 0.33))

ggsave(
  'plots/invasion_rp4.png',
  width = 9,
  height = 12
)
