#----------------------------------
# Purpose of the script: Create Figure 1
#
# What this script does:
# 1. Loads in the data and models for growth rates and plasmid transfer rates
# 2. Creates predictions
# 3. Creates and saves out Figure 1
#
# Author: Dr. Daniel Padfield
# Date Created: 15 April 2026
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
librarian::shelf(tidyverse, brms, tidybayes, patchwork)
## ---------------------------

fmax <- pmax

# read in data for plotting growth rates
d_growth <- read.csv('data/processed/temperature_growth_data.csv')
model_growth <- readRDS(
  'data/models/brms_flextpc.rds'
)

# read in data for plotting plasmid transfer
d_conjugation <- read.csv('data/processed/conjugation_rates.csv')
model_r1 <- readRDS(
  'data/models/brms_flextpc_r1.rds'
)
model_rp4 <- readRDS(
  'data/models/brms_flextpc_rp4.rds'
)

# set colours
cols <- c("#94CBEC", "#C26A77", "#2E2585")
names(cols) <- c('wt', 'wt + R1', 'wt + RP4')

# create TPC graph
preds_growth <- select(d_growth, id2) %>%
  distinct() %>%
  group_by(id2) %>%
  do(data.frame(
    temp = seq(min(d_growth$temp), max(d_growth$temp), length.out = 200)
  )) %>%
  add_epred_draws(
    model_growth,
    ndraws = 1000
  ) %>%
  group_by(id2, temp) %>%
  median_hdci

# for plotting, we only want TPCs for wild-type and wild-type + R1.
d_growth <- filter(d_growth, id2 %in% c('wt', 'wt + R1', 'wt + RP4'))
preds_growth <- filter(preds_growth, id2 %in% c('wt', 'wt + R1', 'wt + RP4'))

# plot predictions
p1 <- ggplot() +
  geom_point(
    aes(temp, gr, col = id2),
    data = d_growth,
    size = 1.5,
    alpha = 0.5
  ) +
  geom_ribbon(
    aes(
      temp,
      ymin = .lower,
      ymax = .upper,
      fill = id2
    ),
    data = preds_growth,
    alpha = 0.3,
    show.legend = FALSE
  ) +
  geom_line(
    aes(temp, .epred, col = id2),
    data = preds_growth,
    show.legend = FALSE,
    linewidth = 1.2
  ) +
  theme_bw(base_size = 14) +
  labs(
    x = 'Temperature (ºC)',
    y = expression("Growth rate (hr"^-1 ~ ")"),
    title = expression(paste('(a) ', italic("E. coli")))
  ) +
  scale_color_manual(
    values = cols,
    labels = c(
      'wt' = 'wild-type',
      'wt + R1' = 'wild-type + R1',
      'wt + RP4' = 'wild-type + RP4'
    )
  ) +
  scale_fill_manual(values = cols) +
  theme(legend.position = c(0.25, 0.8), legend.title = element_blank()) +
  scale_y_continuous(limits = c(0.1, 2.8), n.breaks = 5) +
  scale_x_continuous(limits = c(10, 50), n.breaks = 5) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

p1

# create plasmid transfer graph
# get predictions of the model, marginalising over the random effects
preds_pop_r1 <- select(d_conjugation, strain, rep) %>%
  filter(strain == 'R1') %>%
  distinct() %>%
  group_by(rep) %>%
  do(data.frame(
    temp = seq(
      min(d_conjugation$temp),
      max(d_conjugation$temp),
      length.out = 200
    )
  )) %>%
  add_epred_draws(
    model_r1,
    re_formula = NULL,
    ndraws = 1000
  ) %>%
  group_by(temp, .draw) %>%
  summarise(.epred = mean(.epred)) %>%
  group_by(temp) %>%
  tidybayes::median_hdci() %>%
  ungroup()

# calculate topt from where the rate is highest based on the predictions
# this is a better way than using parameter draws because of the non-linear shapes of each individual curve
select(d_conjugation, strain, rep) %>%
  filter(strain == 'R1') %>%
  distinct() %>%
  group_by(rep) %>%
  do(data.frame(
    temp = seq(
      min(d_conjugation$temp),
      max(d_conjugation$temp),
      length.out = 400
    )
  )) %>%
  add_epred_draws(
    model_r1,
    re_formula = NULL,
    ndraws = 1000
  ) %>%
  group_by(temp, .draw) %>%
  summarise(.epred = mean(.epred)) %>%
  group_by(.draw) %>%
  filter(.epred == max(.epred)) %>%
  ungroup() %>%
  median_hdci(temp)

# get predictions including random effects
preds_group_r1 <- select(d_conjugation, strain, rep) %>%
  filter(strain == 'R1') %>%
  distinct() %>%
  group_by(rep) %>%
  do(data.frame(
    temp = seq(
      min(d_conjugation$temp),
      max(d_conjugation$temp),
      length.out = 200
    )
  )) %>%
  add_epred_draws(
    model_r1,
    re_formula = NULL,
    ndraws = 1000
  ) %>%
  group_by(temp, rep) %>%
  tidybayes::median_hdci() %>%
  ungroup()

# make plots
p2 <- ggplot() +
  geom_ribbon(
    aes(temp, ymin = .lower, ymax = .upper, group = rep),
    preds_group_r1,
    alpha = 0.025,
    fill = cols['wt + R1']
  ) +
  geom_line(
    aes(temp, .epred, group = rep),
    preds_group_r1,
    linewidth = 0.5,
    alpha = 0.3,
    col = cols['wt + R1']
  ) +
  geom_ribbon(
    aes(temp, ymin = .lower, ymax = .upper),
    preds_pop_r1,
    alpha = 0.5,
    fill = cols['wt + R1']
  ) +
  geom_line(
    aes(temp, .epred),
    preds_pop_r1,
    linewidth = 1,
    color = cols['wt + R1']
  ) +
  theme_bw(base_size = 14) +
  labs(
    x = "Temperature (°C)",
    y = expression("Plasmid transfer rate (mL" ~ cell^-1 ~ hr^-1 ~ ")"),
    title = '(b) R1 transfer'
  ) +
  geom_point(
    aes(temp, gamma),
    filter(d_conjugation, strain == 'R1'),
    col = cols['wt + R1'],
    alpha = 0.5
  ) +
  scale_x_continuous(limits = c(10, 50), n.breaks = 5) +
  scale_y_continuous(limits = c(0, 1.5e-11), n.breaks = 6)

p2

# for rp4
preds_pop_rp4 <- select(d_conjugation, strain, rep) %>%
  filter(strain == 'RP4') %>%
  distinct() %>%
  group_by(rep) %>%
  do(data.frame(
    temp = seq(
      min(d_conjugation$temp),
      max(d_conjugation$temp),
      length.out = 200
    )
  )) %>%
  add_epred_draws(
    model_rp4,
    re_formula = NULL,
    ndraws = 1000
  ) %>%
  group_by(temp, .draw) %>%
  summarise(.epred = mean(.epred)) %>%
  group_by(temp) %>%
  tidybayes::median_hdci() %>%
  ungroup()

preds_group_rp4 <- select(d_conjugation, strain, rep) %>%
  filter(strain == 'RP4') %>%
  distinct() %>%
  group_by(rep) %>%
  do(data.frame(
    temp = seq(
      min(d_conjugation$temp),
      max(d_conjugation$temp),
      length.out = 200
    )
  )) %>%
  add_epred_draws(
    model_rp4,
    re_formula = NULL,
    ndraws = 1000
  ) %>%
  group_by(temp, rep) %>%
  tidybayes::median_hdci() %>%
  ungroup()

p3 <- ggplot() +
  geom_ribbon(
    aes(temp, ymin = .lower, ymax = .upper, group = rep),
    preds_group_rp4,
    alpha = 0.025,
    fill = cols['wt + RP4']
  ) +
  geom_line(
    aes(temp, .epred, group = rep),
    preds_group_rp4,
    linewidth = 0.5,
    alpha = 0.3,
    col = cols['wt + RP4']
  ) +
  geom_ribbon(
    aes(temp, ymin = .lower, ymax = .upper),
    preds_pop_rp4,
    alpha = 0.5,
    fill = cols['wt + RP4']
  ) +
  geom_line(
    aes(temp, .epred),
    preds_pop_rp4,
    linewidth = 1,
    color = cols['wt + RP4']
  ) +
  theme_bw(base_size = 14) +
  labs(
    x = "Temperature (°C)",
    y = expression("Plasmid transfer rate (mL" ~ cell^-1 ~ hr^-1 ~ ")"),
    title = '(c) RP4 transfer'
  ) +
  geom_point(
    aes(temp, gamma),
    filter(d_conjugation, strain == 'RP4'),
    col = cols['wt + RP4'],
    alpha = 0.5
  ) +
  scale_x_continuous(limits = c(10, 50), n.breaks = 5) +
  scale_y_continuous(limits = c(0, 1.5e-11), n.breaks = 6)

p3


p1 + p2 + p3

ggsave('plots/figure_1.png', width = 13, height = 4)

# create summary table of key TPC parameters

# grab parameters
# marginalise over random effect for calculating parameters. e.g. calculate random and population level parameters, add them together, and then average those before estimate credible intervals
# this does not work well for topt because of the different shapes and heights of the curves that make the topt estimate poor. This is done better by estimating it from predictions from the underlying data

# grab intercepts
d_r1_intercepts <- spread_draws(model_r1, `b_.*`, regex = TRUE) %>%
  # remove b_ and _Intercept
  rename_with(~ gsub("b_|_Intercept", "", .x), starts_with("b_")) %>%
  pivot_longer(
    cols = c(alphatrans, tmaxtrans, rmaxtrans, tmintrans, beta),
    names_to = 'parameter',
    values_to = 'intercept'
  )

# grab random effects
d_r1_random <- gather_draws(model_r1, `r_rep__.*`, regex = TRUE) %>%
  # grab everything between __ and [
  mutate(
    parameter = gsub("r_rep__|\\[.*\\]", "", .variable),
    # grab everything between [ and ,
    rep = gsub(".*\\[(.*),.*\\]", "\\1", .variable)
  ) %>%
  ungroup() %>%
  select(-.variable) %>%
  rename(random_effect = .value)

d_r1_summary <- left_join(d_r1_random, d_r1_intercepts) %>%
  mutate(value = intercept + random_effect) %>%
  select(-intercept, -random_effect) %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  mutate(
    alpha = boot::inv.logit(alphatrans * 10),
    tmax = tmaxtrans * 100,
    rmax = exp(rmaxtrans * 100),
    tmin = tmintrans * 100,
    # calculate Topt
    topt = alpha * tmax + (1 - alpha) * tmin,
    thermal_tolerance = tmax - tmin
  ) %>%
  select(-c(alphatrans, tmaxtrans, rmaxtrans, tmintrans)) %>%
  pivot_longer(
    cols = c(alpha, tmax, rmax, tmin, topt, thermal_tolerance, beta),
    names_to = 'parameter',
    values_to = 'value'
  ) %>%
  group_by(.chain, .iteration, .draw, parameter) %>%
  summarise(value = mean(value)) %>%
  group_by(parameter) %>%
  median_hdci(value)
# lol why don't these match with Figure 1b

d_r1_summary

# alternative topt derived from predictions
r1_topt <- model_r1$data %>%
  select(rep) %>%
  distinct() %>%
  group_by(rep) %>%
  do(data.frame(
    temp = seq(
      min(model_r1$data$temp),
      max(model_r1$data$temp),
      length.out = 500
    )
  )) %>%
  add_epred_draws(
    model_r1,
    re_formula = NULL,
    ndraws = 1000
  ) %>%
  group_by(temp, .draw) %>%
  summarise(.epred = mean(.epred)) %>%
  group_by(.draw) %>%
  filter(.epred == max(.epred)) %>%
  ungroup() %>%
  median_hdci(value = temp)

# plasmid RP4

# grab intercepts
d_rp4_intercepts <- spread_draws(model_rp4, `b_.*`, regex = TRUE) %>%
  # remove b_ and _Intercept
  rename_with(~ gsub("b_|_Intercept", "", .x), starts_with("b_")) %>%
  pivot_longer(
    cols = c(alphatrans, tmaxtrans, rmaxtrans, tmintrans, beta),
    names_to = 'parameter',
    values_to = 'intercept'
  )

# grab random effects
d_rp4_random <- gather_draws(model_rp4, `r_rep__.*`, regex = TRUE) %>%
  # grab everything between __ and [
  mutate(
    parameter = gsub("r_rep__|\\[.*\\]", "", .variable),
    # grab everything between [ and ,
    rep = gsub(".*\\[(.*),.*\\]", "\\1", .variable)
  ) %>%
  ungroup() %>%
  select(-.variable) %>%
  rename(random_effect = .value)

d_rp4_summary <- left_join(d_rp4_random, d_rp4_intercepts) %>%
  mutate(value = intercept + random_effect) %>%
  select(-intercept, -random_effect) %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  mutate(
    alpha = boot::inv.logit(alphatrans * 10),
    tmax = tmaxtrans * 100,
    rmax = exp(rmaxtrans * 100),
    tmin = tmintrans * 100,
    # calculate Topt
    topt = alpha * tmax + (1 - alpha) * tmin,
    thermal_tolerance = tmax - tmin
  ) %>%
  select(-c(alphatrans, tmaxtrans, rmaxtrans, tmintrans)) %>%
  pivot_longer(
    cols = c(alpha, tmax, rmax, tmin, topt, thermal_tolerance, beta),
    names_to = 'parameter',
    values_to = 'value'
  ) %>%
  group_by(.chain, .iteration, .draw, parameter) %>%
  summarise(value = mean(value)) %>%
  group_by(parameter) %>%
  median_hdci(value)

d_rp4_summary

# alternative topt derived from predictions
rp4_topt <- model_rp4$data %>%
  select(rep) %>%
  distinct() %>%
  group_by(rep) %>%
  do(data.frame(
    temp = seq(
      min(model_rp4$data$temp),
      max(model_rp4$data$temp),
      length.out = 500
    )
  )) %>%
  add_epred_draws(
    model_rp4,
    re_formula = NULL,
    ndraws = 1000
  ) %>%
  group_by(temp, .draw) %>%
  summarise(.epred = mean(.epred)) %>%
  group_by(.draw) %>%
  filter(.epred == max(.epred)) %>%
  ungroup() %>%
  median_hdci(value = temp)

# summary stats of TPCs

# grab parameters
d_tpc <- gather_draws(model_growth, `b_.*`, regex = TRUE) %>%
  # remove b_ and _Intercept
  mutate(.variable = gsub("b_|_intercept", "", .variable)) %>%
  separate(.variable, into = c('parameter', 'id2'), sep = '_') %>%
  mutate(id2 = gsub('id2', '', id2)) %>%
  pivot_wider(names_from = parameter, values_from = .value) %>%
  # back transform parameters
  mutate(
    alpha = boot::inv.logit(alphatrans * 10),
    tmax = tmaxtrans * 100,
    rmax = exp(rmaxtrans),
    tmin = tmintrans * 100,
    # calculate Topt
    topt = alpha * tmax + (1 - alpha) * tmin,
    thermal_tolerance = tmax - tmin
  )

unique(d_tpc$id2)

# get summary of parameters
d_tpc_summary <- d_tpc %>%
  pivot_longer(
    cols = c(alpha, tmax, rmax, tmin, topt, thermal_tolerance, beta),
    names_to = 'parameter',
    values_to = 'value'
  ) %>%
  group_by(id2, parameter) %>%
  median_hdci(value)

d_tpc_summary %>%
  filter(parameter == 'thermal_tolerance') %>%
  filter(str_detect(id2, 'wt'))

# make a table of parameters for all the TPCs.
table_bact <- d_tpc_summary %>%
  filter(str_detect(id2, 'wt')) %>%
  mutate(
    id2 = case_when(
      id2 == 'wtPR1' ~ 'wt + R1',
      id2 == 'wtPRP4' ~ 'wt + RP4',
      TRUE ~ id2
    )
  ) %>%
  select(id2, parameter, value, .lower, .upper)

table_R1 <- d_r1_summary %>%
  filter(parameter != 'topt') %>%
  bind_rows(r1_topt %>% mutate(parameter = 'topt')) %>%
  mutate(id2 = 'R1 transfer') %>%
  select(id2, parameter, value, .lower, .upper)

table_RP4 <- d_rp4_summary %>%
  filter(parameter != 'topt') %>%
  bind_rows(rp4_topt %>% mutate(parameter = 'topt')) %>%
  mutate(id2 = 'RP4 transfer') %>%
  select(id2, parameter, value, .lower, .upper)

# Custom formatting function
format_num <- function(x, threshold = 0.001) {
  ifelse(abs(x) < threshold, sprintf("%.1e", x), round(x, 2))
}

table1 <- bind_rows(table_bact, table_R1, table_RP4) %>%
  mutate(
    parameter = case_when(
      parameter == 'thermal_tolerance' ~ 'thermal tolerance',
      TRUE ~ parameter
    ),
    across(where(is.numeric), \(x) format_num(x, threshold = 0.001))
  ) %>%
  flextable(,
    col_keys = c(
      'id2',
      'parameter',
      'value',
      '.lower',
      '.upper'
    )
  ) %>%
  set_header_labels(
    id2 = 'Thermal performance',
    parameter = 'TPC parameter',
    value = 'Median estimate',
    .lower = 'Lower 95% CI',
    .upper = 'Upper 95% CI'
  ) %>%
  merge_v(~id2) %>%
  valign(valign = 'top', j = 1, part = 'body') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 14, part = 'all') %>%
  align(align = 'center', part = 'all') %>%
  hline(i = c(7, 14, 21, 28), border = fp_border_default()) %>%
  fix_border_issues()

# save this out
save_as_image(table1, 'tables/table_tpc.png', webshot = 'webshot2')
save_as_docx(table1, path = 'tables/table_tpc.docx')
