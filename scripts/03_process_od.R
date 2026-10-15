# ---------------------------
# Purpose of script: process OD data
#
# What this script does:
# 1. Read in all data
# 2. Correct for blanks
# 3. Remove blanks
# 4. Filter out biologically implausible data
# 5. Save out file
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
librarian::shelf(tidyverse)

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
  #pattern = '01_Trans'
  pattern = '01_Kanamycin_DR'
)

# read in data
d <- file_names %>% map_df(read.csv)

# check the number of unique combinations there should be
select(d, file, well, id, serial_no) %>%
  distinct() %>%
  nrow()

# split file identifier
# what do the bits between the underscores in your file mean?
set_file_id <- c('project', 'temp', 'run')
d <- separate(
  d,
  file,
  into = c('project', 'temp', 'run'),
  sep = '_',
  remove = FALSE
)

# set up your grouping variables to get each individual well for each plate
# this bit is interactive
groupings <- c("file", "run", "well", "serial_no", "id", set_file_id)

# set blank method: options are: well_specific or blank_median
blank_method <- "blank_median"

# identity of blanks: what is the unique ID you give to your blank wells
blank_id <- 'B'

# output name
#output <- 'Trans_gr.csv'
output <- 'Kanamycin_DR_gr.csv'
output <- paste('03_', output, sep = '')

#---------------------#

#--------------------------------------------------------------#
# remove samples because of biologically implausible signal ####
#--------------------------------------------------------------#

# from the plots in first_look_plots, you can remove some wells that have biologically implausible signal (e.g. massive spikes) that are not going to be easy to model and will give poor estimates
# if some of the control/blank wells are wrong you can remove them too

# create unique ID of file and well
d_filt <- d %>%
  mutate(plate_well_id = paste(file, well, sep = '_'))

# remove uninoculated rows ('X' in plate layout)
d_filt <- d_filt %>%
  filter(id != 'X')

#---------------------#
# calculate blanks ####
#---------------------#

# method 1: create a well-specific blank
if (blank_method == "well_specific") {
  # take the first three time points
  # average per grouping
  d_blank <- d_filt %>%
    group_by(across(all_of(groupings))) %>%
    filter(measurement_time_hr %in% sort(measurement_time_hr)[1:3]) %>%
    summarise(
      ave_blank = mean(raw_absorbance),
      sd_blank = sd(raw_absorbance),
      .groups = 'drop'
    )
}

# method 2 Calculate blank from average of all blank wells in each plate

# name outside wells
# fmt: skip

# do blank subtraction
if (blank_method == "blank_median") {
  d_blank <- d_filt %>%
    filter(substr(id, 1, 1) == blank_id) %>%
    group_by(across(all_of(groupings[!groupings %in% c('well', 'id')]))) %>%
    summarise(
      ave_blank = median(raw_absorbance),
      sd_blank = sd(raw_absorbance),
      .groups = 'drop'
    )
}

# visualise average blank
# if the SDs are very big, you might want to check why!
hist(d_blank$ave_blank)
hist(d_blank$sd_blank)

#-----------------#
# wrangle data ####
#-----------------#

# add the blank to the data
d_filt <- left_join(d_filt, d_blank)

# make corrected absorbance
d_filt <- mutate(d_filt, od_cor = raw_absorbance - ave_blank)

# set negative values to the minimum value for that well
d_filt <- group_by(d_filt, file, well) %>%
  mutate(
    od_cor = ifelse(od_cor <= 0, min(od_cor[od_cor > 0], na.rm = TRUE), od_cor)
  ) %>%
  ungroup()

# remove all wells that are blanks - need to edit this code maybe
d_filt <- filter(d_filt, substr(id, 1, 1) != blank_id)

#-----------------------------#
# optional broad filtering ####
#-----------------------------#

# In 01_kanamycin_selection.pdf, the three 19ºC plates that run for >72 hours have a very weird kink in the data after >72 hours. This cannot be explained by growth dynamics so we will remove these. It is very odd that is so repeatable.
d_filt <- filter(d_filt, !(temp == '19' & measurement_time_hr > 72))

#------------------#
# save out data ####
#------------------#

write.csv(d_filt, file.path('data', 'processed', output))
