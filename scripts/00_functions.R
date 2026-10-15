# ---------------------------
# Purpose of script: Contains functions to help load in excel files from the Multiskan Skyhigh
#
# What this script does:
# 1. Function to load in Raw data spreadsheet
# 2. Function to load in General info spreadsheet
# 3. Function to load in Run info spreadsheet
# 4. Function to letter facets in a ggplot2 panel
# 5. Function to convert size of text in ggplot2 to pts
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
# All functions add the filename as a variable
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
# custom functions ####
#---------------------#

# function to parse od data from Sheet 2
get_od <- function(file) {
  # read in temporary file
  temp <- read_excel(file, sheet = 2, skip = 4) %>%
    # clean file names up
    janitor::clean_names()

  # grab the start time from the file
  start_time <- read_excel(file, sheet = 2, range = "A1:A2") %>%
    pull(1) %>%
    mdy_hms(tz = 'GMT')

  # calculate time of each measurement
  temp$time <- start_time
  temp$time <- temp$time + lubridate::seconds(temp$measurement_time_s)

  # add file name without extension
  temp$file <- tools::file_path_sans_ext(basename(file))

  # rearrange dataset
  temp <- temp %>%
    select(file, measurement_time_s, time, everything())

  return(temp)
}

# function to parse general information from Sheet 4
get_info <- function(file) {
  # read in temporary file
  temp <- read_excel(
    file,
    sheet = 3,
    range = "C3:D10",
    col_names = c('param', 'value')
  )
  temp2 <- read_excel(
    file,
    sheet = 3,
    range = "C17:D23",
    col_names = c('param', 'value')
  )

  temp <- bind_rows(temp, temp2) %>%
    mutate(file = tools::file_path_sans_ext(basename(file))) %>%
    select(file, everything())

  return(temp)
}

# function to grab the data from the Run log (Sheet 5)
get_runlog <- function(file) {
  # read in temporary file
  temp <- suppressWarnings(read_excel(
    file,
    sheet = 5,
    .name_repair = 'minimal'
  )) %>%
    janitor::clean_names() %>%
    select(2:3) %>%
    filter(!is.na(1))
  colnames(temp) <- c('temperature', 'time')

  # add file name without extension
  temp$file <- tools::file_path_sans_ext(basename(file))
  temp <- select(temp, file, everything())

  return(temp)
}

# function to grab the plate layout from the Run log (Sheet 5)
get_plate_layout <- function(file) {
  # read in temporary file
  temp <- suppressWarnings(readxl::read_excel(
    file,
    sheet = 4,
    .name_repair = 'minimal'
  )) %>%
    rename(well_letter = 1) %>%
    pivot_longer(
      cols = 2:ncol(.),
      values_to = 'id',
      names_to = 'well_number'
    ) %>%
    unite('well', well_letter, well_number, sep = '')

  # add file name without extension
  temp$file <- tools::file_path_sans_ext(basename(file))
  temp <- select(temp, file, everything())

  return(temp)
}

# create function to add letter to facet labels
letter_facets <- function(string) {
  len <- length(string)
  string <- paste("(", letters[1:len], ") ", string, sep = "")
  return(string)
}

# setting pts in a ggplot
pts <- function(x) {
  as.numeric(grid::convertX(grid::unit(x, "points"), "mm"))
}

contrast <- function(colour) {
  out <- rep("black", length(colour))
  light <- farver::get_channel(colour, "l", space = "hcl")
  out[light < 60] <- "white"
  out
}
