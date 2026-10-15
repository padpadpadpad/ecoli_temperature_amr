# ---------------------------
# Purpose of script: Calculate growth rate metrics from processed data
#
# What this script does:
# 1. Read in processed data from 03_process_od.R
# 2. Remove noise at low ODs. It does this in two steps. If you have 50% or more low ODs, we leave all the data in, if you have 50% or less low ODs, we remove the low ODs. This is because sometimes the low ODs are noise, but sometimes they are real (e.g. no growth ever). The characterisation exists in the column threshold_density
# 3. Calculate derivative for each well
# 4. Calculate growth metrics and save them out
# 5. Visualise growth metrics on growth curves
#
# Author: Dr. Daniel Padfield
#
# Date Created:  2025-02-20
#
# Copyright (c) Daniel Padfield, 2025
# This code is licensed under an MIT license and has been developed with the help of GitHub copilot, and potentially other LLMs. However, all output has been checked and modified by a human and we take full responsibility for any errors.
#
# ---------------------------
#
# Notes:
#
# ---------------------------

# if librarian is not installed, install it
if (!requireNamespace("librarian", quietly = TRUE)) {
  install.packages("librarian")
}

# if BiocManager is not installed, install it
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# if Biobase is not installed, install it from Bioconductor
if (!requireNamespace("Biobase", quietly = TRUE)) {
  BiocManager::install("Biobase")
}

# load packages
# fmt: skip
librarian::shelf(tidyverse, growthrates, gcplyr)

## ---------------------------

#---------------------#
# things to change ####
#---------------------#
# The processing scripts (01_read_od - 05_clean_data) needs to be run separately
# for the dose response growth rates (DR sheets) and the transconjugant growth rates (T_gr).
# The necessary changes will be coded in the script with a #comment

# name of files to read in
file_names <- list.files(
  'data/processed',
  full.names = TRUE,
  #pattern = '03_Trans_gr'
  pattern = '03_Kanamycin_DR'
)

# read in data
d <- file_names %>% map_df(read.csv)

# split file identifier if needed, might have been done in 03_process_od.R

# check the number of unique combinations there should be
select(d, file, well, id, serial_no) %>%
  distinct() %>%
  nrow()

# set up your grouping variables to get each individual well for each plate
# this bit is interactive
groupings <- c("file", "well", "serial_no", "id", 'temp', 'run')

# output name
#output <- 'Trans_gr_metrics.csv'
output <- 'Kanamycin_DR_gr_metrics.csv'
output <- paste('04_', output, sep = '')

# window_width - set to calculate a rolling regression over the course of an hour
# this measurement window code might need to be changed?
measurements_every_X_min <- 10
window_width <- 60 / measurements_every_X_min

# if window_width is even, subtract 1
if (window_width %% 2 == 0) {
  window_width <- window_width - 1
}

rm(measurements_every_X_min)


#-----------------#
# wrangle data ####
#-----------------#

# This needs to be treated with care, think some wells have no measurements over this threshold, but that means they did not grow. This is important data!
# There should be 2 modelling steps then, if max density is close to the limit, leave the points in, if not, then remove them

# set threshold
threshold <- 0.0090

# for each well, calculate the proportion of points over the filtering threshold
d_check <- group_by(d, file, well) %>%
  summarise(
    n_over = sum(od_cor > threshold, na.rm = TRUE),
    n_total = n(),
    prop_over = n_over / n_total,
    .groups = 'drop'
  )

# do a plot
ggplot(d_check, aes(x = prop_over)) +
  geom_histogram(binwidth = 0.05) +
  theme_bw() +
  labs(x = 'Proportion of points over threshold', y = 'Count') +
  NULL

# lets say at least 50% of points should be over the threshold, and if not then we need to treat the growth rate estimates with caution - this is arbitrary
d_check <- mutate(
  d_check,
  over_threshold = ifelse(prop_over >= 0.5, TRUE, FALSE)
)

# run two filtering steps
# if you have 50% of your points above the threshold, then filter out the below it.
d_filt1 <- left_join(d, select(d_check, file, well, over_threshold)) %>%
  filter((od_cor > threshold & over_threshold == TRUE)) %>%
  select(-over_threshold) %>%
  mutate(threshold_density = 'mostly_above')
# if you have less than 50% of your points above the threshold, then do no filtering.
d_filt2 <- left_join(d, select(d_check, file, well, over_threshold)) %>%
  filter(over_threshold == FALSE) %>%
  select(-over_threshold) %>%
  mutate(threshold_density = 'mostly_below')

# recombine the datasets
d_filt <- bind_rows(d_filt1, d_filt2)

# check it has worked
group_by(d_filt, file, well) %>%
  summarise(
    n_over = sum(od_cor > threshold, na.rm = TRUE),
    n_total = n(),
    prop_over = n_over / n_total,
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = prop_over)) +
  geom_histogram(binwidth = 0.05) +
  theme_bw() +
  labs(x = 'Proportion of points over threshold', y = 'Count')

# quick plot of data ####
d_filt %>%
  ggplot(aes(
    x = measurement_time_hr,
    y = od_cor,
    group = interaction(well)
  )) +
  geom_line(alpha = 0.3) +
  facet_wrap(~file) +
  theme_bw()

# calculate the expected number of rows in the output - 1 per group
n_expected <- select(d_filt, all_of(groupings)) %>%
  distinct() %>%
  nrow()
n_expected

# look at distribution of measurement times
group_by(d_filt, file) %>%
  summarise(
    max_time = max(measurement_time_hr, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  group_by(max_time) %>%
  tally()

#----------------------------------------#
# calculate model-free growth metrics ####
#----------------------------------------#

# this uses gcplyr: https://mikeblazanin.github.io/gcplyr/index.html

# calculate per-capita derivative - this is the slope across a window of time
# a key parameter here is window_width_n, which is the number of time points to use in the calculation

# add threshold_density to groupings
groupings <- c(groupings, 'threshold_density')

# this group by is the important step here, needs to contain all the variables to get a single well
d_filt <- d_filt %>%
  group_by(across(all_of(groupings))) %>%
  mutate(
    deriv_percap = calc_deriv(
      x = measurement_time_hr,
      y = od_cor,
      percapita = TRUE,
      blank = 0,
      window_width_n = window_width,
      trans_y = "log"
    )
  ) %>%
  ungroup()

# calculate growth metrics
d_sum <- d_filt %>%
  group_by(across(all_of(groupings))) %>%
  summarise(
    min_dens = first_minima(od_cor, return = 'y'),
    lag_time = lag_time(
      y = od_cor,
      x = measurement_time_hr,
      deriv = deriv_percap,
      blank = 0,
      y0 = min_dens
    ),
    gr = max_gc(deriv_percap, na.rm = TRUE),
    gr_time = extr_val(measurement_time_hr, which_max_gc(deriv_percap)),
    gr_dens = extr_val(od_cor, which_max_gc(deriv_percap)),
    max_dens = max_gc(od_cor, na.rm = TRUE),
    max_time = extr_val(measurement_time_hr, which_max_gc(od_cor)),
    doub_time = doubling_time(y = gr),
    final_time = max(measurement_time_hr, na.rm = TRUE),
    final_dens = od_cor[which.max(measurement_time_hr)],
    auc = auc(y = od_cor, x = measurement_time_hr),
    auc_48 = tryCatch(
      auc(y = od_cor, x = measurement_time_hr, xlim = c(0, 48)),
      error = function(e) Inf
    ),
    .groups = 'drop'
  )

# run a check to see if the number of rows is as expected
if (nrow(d_sum) == n_expected) {
  message("✅ Dataframe has the expected number of rows: ", n_expected)
} else {
  message("⚠️ Warning: Expected ", n_expected, " rows, but found ", nrow(d_sum))
}

# save this out ####
write.csv(
  d_sum,
  file.path('data/processed', output),
  row.names = FALSE
)

hist(d_sum$auc)

#----------------------#
# Visualise metrics ####
#----------------------#

# make a plot to visualise every plate
files <- unique(d$file)

# if metrics_plots does not exist, create it
if (!dir.exists('plots/metric_plots')) {
  dir.create('plots/metric_plots', recursive = TRUE)
}

# open a pdf
pdf(
  file.path(
    'plots/metric_plots',
    paste('check_gcplyr_log_Kanamycin_DR_gr', '.pdf', sep = '')
  ),
  #paste('check_gcplyr_log_Trans_gr', '.pdf', sep = '')),
  width = 10,
  height = 6.5
)

for (i in 1:length(files)) {
  temp_od <- filter(d, file == files[i]) |>
    mutate(column = str_extract(well, "[A-Z]+"), row = parse_number(well))

  temp_params <- filter(d_sum, file == files[i]) %>%
    mutate(column = str_extract(well, "[A-Z]+"), row = parse_number(well))

  temp_id <- select(temp_od, column, row, id) %>%
    distinct() %>%
    mutate(x = 0, y = 1)

  # fmt: skip
  temp_plot <- temp_od %>%
    ggplot(aes(x = measurement_time_hr, y = log(od_cor))) +
    geom_point(size = 0.1) +
    geom_vline(aes(xintercept = lag_time), temp_params, col = 'red') +
    geom_point(aes(gr_time, log(gr_dens)), temp_params, col = 'red') +
    geom_point(aes(max_time, log(max_dens)), temp_params, col = 'blue') +
    geom_label(aes(x = 0, y = 1.25, label = id), temp_id, hjust = 0, vjust = 0.8, label.size = NA) +
    facet_grid(column~row, switch = 'y') +
    theme_bw() +
    labs(title = paste('File:', files[i], sep = ' '),
         x = 'Time (hr)',
         y = 'log Absorbance (OD600)') +
    ylim(c(min(log(temp_od$od_cor)[is.finite(log(temp_od$od_cor))]), 1.25)) +
    NULL

  print(temp_plot)
}

dev.off()

# do the same plot but unlogged
# open a pdf
pdf(
  file.path(
    'plots/metric_plots',
    #paste('check_gcplyr_nolog_Trans_gr', '.pdf', sep = '')
    paste('check_gcplyr_nolog_Kanamycin_DR_gr', '.pdf', sep = '')
  ),
  width = 10,
  height = 6.5
)

for (i in 1:length(files)) {
  temp_od <- filter(d, file == files[i]) |>
    mutate(column = str_extract(well, "[A-Z]+"), row = parse_number(well))

  temp_params <- filter(d_sum, file == files[i]) %>%
    mutate(column = str_extract(well, "[A-Z]+"), row = parse_number(well))

  temp_id <- select(temp_od, column, row, id) %>%
    distinct() %>%
    mutate(x = 0, y = 1)

  # fmt: skip
  temp_plot <- temp_od %>%
    ggplot(aes(x = measurement_time_hr, y = od_cor)) +
    geom_point(size = 0.1) +
    geom_vline(aes(xintercept = lag_time), temp_params, col = 'red') +
    geom_point(aes(gr_time, gr_dens), temp_params, col = 'red') +
    geom_point(aes(max_time, max_dens), temp_params, col = 'blue') +
    geom_label(aes(x = 0, y = 1.25, label = id), temp_id, hjust = 0, vjust = 0.8, label.size = NA) +
    facet_grid(column~row, switch = 'y') +
    theme_bw() +
    labs(title = paste('File:', files[i], sep = ' '),
         x = 'Time (hr)',
         y = 'log Absorbance (OD600)') +
    ylim(c(min(temp_od$od_cor[is.finite(temp_od$od_cor)]), 1.25)) +
    NULL

  print(temp_plot)
}

dev.off()
