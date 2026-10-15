# ---------------------------
# Purpose of script: Create a visualisations for raw OD readings from the Multiskan Skyhigh
#
# What this script does:
# 1. Loads in data
# 2. Creates a plot for each file, with all wells combined
#
# Author: Dr. Daniel Padfield
#
# Date Created:  2025-01-06
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
librarian::shelf(tidyverse)

## ---------------------------

#---------------------#
# things to change ####
#---------------------#
# The processing scripts (01_read_od - 05_clean_data) needs to be run separately
# for the dose response growth rates (DR sheets) and the transconjugant growth rates (T_gr).
# The necessary changes will be coded in the script with a #comment

# set run, this will be used to label the plots
# should be the same as output from 01_process_od.R

#input <- 'Trans_gr.csv'
input <- 'Kanamycin_DR_gr.csv'
input <- paste('01_', input, sep = '')

input_no_ext <- tools::file_path_sans_ext(input)

# change this as needed
d_od <- read.csv(file.path('data/processed', input))

#---------------------#

# plot each well across each combination of factors ####

# create an ID column for each plate
d_od <- d_od %>%
  mutate(plate_id = paste(file, serial_no, sep = '_'))

plate_ids <- unique(d_od$plate_id)

# create a director if plots/first_look_plots does not exist
if (!dir.exists('plots/first_look_plots')) {
  dir.create('plots/first_look_plots', recursive = TRUE)
}

# open a pdf
pdf(
  file.path(
    'plots/first_look_plots',
    paste(input_no_ext, '.pdf', sep = '')
  ),
  width = 10,
  height = 6.5
)

for (i in 1:length(plate_ids)) {
  temp_od <- d_od |>
    filter(plate_id == plate_ids[i]) |>
    mutate(column = str_extract(well, "[A-Z]+"), row = parse_number(well))

  temp_id <- select(temp_od, column, row, id) %>%
    distinct() %>%
    mutate(x = 0, y = 1)

  temp_plot <- temp_od %>%
    ggplot(aes(x = measurement_time_hr, y = raw_absorbance)) +
    geom_line() +
    geom_text(aes(x = 0, y = 1, label = id), temp_id, hjust = 0) +
    facet_grid(column ~ row, switch = 'y') +
    theme_bw() +
    labs(
      title = paste('Plate ID:', plate_ids[i]),
      x = 'Time (hr)',
      y = 'Absorbance (OD600)'
    ) +
    ylim(c(0, 1.25))

  print(temp_plot)
}

dev.off()

# remove temporary items in plot
rm(list = c('temp_od', 'temp_plot'))
