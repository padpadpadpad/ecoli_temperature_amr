# ---------------------------
# Purpose of script: Process files from Multiskan Skyhigh
#
# What this script does:
# 1. Lists files ready to be merged, this could be changed to improve file selection (e.g. only do "run1")
# 2. Calculates measurement time to add to the second files
# 3. Removes outer wells
# 4. Adds in correct master plate
# 5. Save out
#
# Author: Dr. Daniel Padfield
#
# Date Created:  2024-12-16
#
# Copyright (c) Daniel Padfield, 2024
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
librarian::shelf(tidyverse, readxl)

## ---------------------------

# load in extra functions
source('scripts/00_functions.R')

#---------------------#
# things to change ####
#---------------------#
# The processing scripts (01_read_od - 05_clean_data) needs to be run separately
# for the dose response growth rates (DR sheets) and the transconjugant growth rates (T_gr).
# The necessary changes will be coded in the script with a #comment

# list files in
# change this to only be the correct files
files <- list.files(
  'data/raw/growth_rates/',
  full.names = TRUE,
  #pattern = 'Tgr.*\\.xlsx'
  #for dose response # the below and unhash the above
  pattern = '^DR'
)

# set output name
#output <- 'Trans_gr.csv'
output <- 'Kanamycin_DR_gr.csv'
output <- paste('01_', output, sep = '')

#-----------------#
# read in data ####
#-----------------#

# remove file if there is a $ in the name - means they are temporary
files <- files[!grepl('\\$', files)]

# bind together all the OD data using get_od() - can see how it works in 00_functions.R
d_od <- map_df(files, get_od) %>%
  mutate(well = gsub(" ", "", well))

# bind all info data using get_info() - can see how it works in 00_functions.R
d_info <- map_df(files, get_info)

# bind all the runlog data - using get_runlog() - can see how it works in 00_functions.R
d_runlog <- map_df(files, possibly(get_runlog))
# check which files are not in d_runlog
no_runlog <- files[
  !tools::file_path_sans_ext(basename(files)) %in% unique(d_runlog$file)
]

# read in master plates
# d_master_plates
d_layout <- map_df(files, get_plate_layout)

#--------------------#
# wrangle od data ####
#--------------------#

# bind well data with od data
d_od <- left_join(d_od, d_layout)

# grab serial number from the info file
d_serial <- filter(d_info, str_detect(param, 'Serial number')) %>%
  select(file, serial_no = value)

# bind with od
d_od <- left_join(d_od, d_serial)

# split the file name of d_od to its useful parts
d_od <- d_od %>%
  select(
    file,
    serial_no,
    measurement_time_s,
    wavelength_nm,
    well,
    raw_absorbance,
    id
  )

# convert time to minutes and hours
d_od <- mutate(
  d_od,
  measurement_time_min = measurement_time_s / 60,
  measurement_time_hr = measurement_time_min / 60
)

# save out data
write.csv(
  d_od,
  file.path('data/processed', output),
  row.names = FALSE
)
