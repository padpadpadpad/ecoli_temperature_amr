#----------------------------------
# Purpose of the script: Tries to identify wells that have dried up and have that are de novo evolved
#
# What this script does:
# 1. Reads in growth metrics
# 2. Cleans data using a series of rules discussed below
#
# 3. Saves ou data
#
# Author: Dr. Daniel Padfield
# Date Created: 01 August 2025
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
librarian::shelf(tidyverse, ggbeeswarm, janitor, easystats, readxl)
## ---------------------------

# extra functions

# write function to return vector identifying outliers
# can change method used for scoring outliers, and pass extra arguments to performance::check_outliers
get_outliers <- function(x, method = 'zscore_robust', ...) {
  temp <- performance::check_outliers(
    x,
    method = 'zscore_robust',
    ...
  )
  return(
    as.data.frame(temp) %>%
      janitor::clean_names() %>%
      pull(outlier)
  )
}

# 1. check numbers in each grouping (temp x row/concentration x id) to look at number of replicates. Should be 4 for each combination, but could be more
# 2. If the value of max density is 0, classify as gcplyr has failed
# 3. Where wells were mostly below threshold density, set gr, lag time and doub time to 0 as they are not reliable
# 4. Try identified dried wells by looking at where OD drops by more than 40% between the max od and the max density is greater than 0.075. Dried wells are harder to identify when the max density is low
# 5. Try and identify de novo evolution in nal.5 and wt wells by looking for outliers in final density within each temp x row/concentration grouping. If you are an outlier with a high final density, classify as de novo evolution. If you are an outlier with a low final density, classify as failed.
# 6. If you are an outlier, but still no greater than 0.15 from the median final density, classify as fine
# 7. For all strains, if you are lower than 0.5 the median final density, classify as failed
# 8. Using our expertise of the system, we also manually coded what we thought were dried up wells and de novo evolution. We compared the programmatic classifications to these and two people manually checked any disagreements. For it to be failed, both people had to agree it was failed.
# 9. Update disagreements based on manual checking

# read in data
d_metrics <- read.csv(
  'data/processed/04_Kanamycin_DR_gr_metrics.csv'
) %>%
  mutate(n = 1:n())
d_processed <- read.csv(
  'data/processed/03_Kanamycin_DR_gr.csv'
)
d_dried <- read_excel(
  'data/processed/wells_to_discard.xlsx'
)

# remove wells with known problems ####

# Zoltan forgot to inoculate one row
d_metrics <- mutate(
  d_metrics,
  classification = ifelse(
    (file == 'DR1907_19_1' & substr(well, 1, 1) == 'H'),
    'failed',
    'fine'
  )
)

# check for NAs in growth
filter(d_metrics, is.na(gr)) %>%
  nrow()
filter(d_metrics, is.infinite(gr)) %>% nrow()
filter(d_metrics, is.infinite(max_dens)) %>%
  nrow()

filter(d_metrics, is.infinite(max_dens)) %>%
  filter(classification == 'fine')

# replace NAs and infinites with 0s
d_metrics <- mutate(
  d_metrics,
  across(where(is.numeric), ~ replace_na(., 0)),
  across(where(is.numeric), ~ replace(., is.infinite(.), 0))
)

# classify well where max dens is 0 as they have failed
filter(d_metrics, max_dens == 0) %>%
  nrow()
d_metrics <- mutate(
  d_metrics,
  classification = ifelse(
    (max_dens == 0 & classification == 'fine'),
    'gcplyr_failed',
    classification
  )
)

# where wells were mostly below threshold_density, set some of the metrics to 0 as they are not reliable.
d_metrics <- mutate(
  d_metrics,
  across(c(gr, lag_time, doub_time), function(x) {
    ifelse(threshold_density == 'mostly_below', 0, x)
  })
)

# Check total number of reps for each temp x row x id combination, should be 4 each
d_metrics %>%
  mutate(row = parse_number(well)) %>%
  group_by(temp, row, id) %>%
  tally() %>%
  filter(n < 4) %>%
  group_by(temp) %>%
  tally()

# count number of each type in classification
d_metrics %>%
  group_by(classification) %>%
  tally() %>%
  arrange(desc(n))

#-------------------------------------------------#
# try and filter out where wells have dried up ####
#-------------------------------------------------#

# our hypothesis is that they reduce their OD, this could be because they've dried up, or because some of the cells have died
# we think that the larger drops are because of drying up

# we will label wells where OD drops by more than 40% between the max od and the max density is greater than 0.075
# drying is harder to identify when the max density is low
d_metrics <- mutate(
  d_metrics,
  classification = ifelse(
    ((final_dens / max_dens) < 0.6 & max_dens > 0.075),
    'dried',
    classification
  )
)

# check number in each classification
d_metrics %>%
  group_by(classification) %>%
  tally() %>%
  arrange(desc(n))

#-------------------------------------------------------#
# filter wells that have undergone de novo evolution ####
#-------------------------------------------------------#

# some of the wells have had de novo evolution to occur. This typically occurs at high antibiotic concentrations, and only in the wt or nal.5 wells
# it is characterised by an increase in the lag time and an increase in the maximum and final density of the wells in which it occurs, but this may be tricky to characterise
# there is no perfect method to this

# keep only nal.5 or wt cells that are fine
d2 <- filter(d_metrics, str_detect(id, 'nal.5|wt')) %>%
  filter(classification == 'fine') %>%
  mutate(column = str_extract(well, "[A-Z]+"), row = parse_number(well))
unique(d2$id)

# how many have already lost wells?
d2 %>%
  group_by(temp, row) %>%
  tally() %>%
  filter(n < 6) %>%
  arrange(n)

# look at distribution of lag times, growth rate and max densities at each antibiotic concentration/column
# this will help us see what metric to use for identifying de novo evolution
select(d2, row, id, temp, lag_time, gr, max_dens, final_dens) %>%
  pivot_longer(
    cols = c(lag_time, gr, max_dens, final_dens),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(!is.infinite(value)) %>%
  ggplot(aes(x = row, y = value, col = id)) +
  geom_beeswarm() +
  theme_bw() +
  facet_grid(metric ~ temp, scales = "free_y")
# looks like final or max density might be the best metric to use

# check for "outliers" in the distribution using: https://easystats.github.io/performance/articles/check_outliers.

# calculate outliers for all groupings
# the threshold here can change based on how strict you want to be in removing outliers
# remove any point that is 0.5 or further from median final density, likely a failed well
# any point less than 0.5 of the median final density has probably not failed and is not de novo evolution, make those fine
# any point within 0.15 of the median final density is also fine, positive or negative
# these are kind of arbitrary but we're using our "domain knowledge"
d2 <- d2 %>%
  group_by(temp, row) %>%
  mutate(
    outlier = get_outliers(final_dens),
    classification = ifelse(outlier == 1, 'denovo', 'fine'),
    classification = ifelse(
      (final_dens < median(final_dens, na.rm = TRUE) &
        classification == 'denovo'),
      'failed',
      classification
    ),
    classification = ifelse(
      (final_dens < median(final_dens, na.rm = TRUE) &
        abs(final_dens - median(final_dens)) < 0.5),
      'fine',
      classification
    ),
    classification = ifelse(
      abs(final_dens - median(final_dens)) < 0.15,
      'fine',
      classification
    )
  ) %>%
  ungroup()

d2 %>%
  group_by(classification) %>%
  tally()

filter(d2) %>%
  ggplot(
    .,
    aes(
      x = forcats::fct_reorder(as.character(row), row),
      y = final_dens,
      col = classification
    )
  ) +
  geom_beeswarm() +
  facet_wrap(~temp)

# combine back with other data
d_metrics <- filter(d_metrics, !n %in% d2$n) %>%
  bind_rows(., select(d2, -c(column, row, outlier)))

# set the same classifications (for low final density) for the plasmid carriers
d3 <- filter(d_metrics, str_detect(id, 'R1|RP4')) %>%
  mutate(column = str_extract(well, "[A-Z]+"), row = parse_number(well)) %>%
  filter(classification == 'fine') %>%
  group_by(temp, row) %>%
  mutate(
    classification = ifelse(
      (final_dens < median(final_dens, na.rm = TRUE) &
        abs(final_dens - median(final_dens)) > 0.5),
      'failed',
      classification
    )
  ) %>%
  ungroup()

group_by(d3, classification) %>%
  tally()

# add these metrics back in
d_metrics <- d_metrics %>%
  filter(!n %in% d3$n) %>%
  bind_rows(., select(d3, -c(column:row)))

# check numbers
d_metrics %>%
  group_by(classification) %>%
  tally()

#----------------------------------------#
# check with Zoltan's classifications ####
#----------------------------------------#

# clean up data in dried column
d_dried <- janitor::clean_names(d_dried) %>%
  mutate(
    file = paste(plate_id, temp, run, sep = '_'),
    id = paste(file, well, sep = '_')
  ) %>%
  filter(parse_number(well) != 12) %>%
  filter(id %in% paste(d_metrics$file, d_metrics$well, sep = '_'))

# if well is in d_dried classify it as failed by Zoli
# if classification is fine and zoli told us it failed, classify it as zoli_failed
# if classification is not fine, and zoli told us it was fine, classify as zoli_fine
d_metrics <- mutate(
  d_metrics,
  id2 = paste(file, well, sep = '_'),
  classification2 = classification,
  zoli_remove = ifelse(
    id2 %in% d_dried$id,
    'failed',
    'fine'
  ),
  classification = case_when(
    (zoli_remove == 'failed' & classification == 'fine') ~ 'zoli_failed',
    (zoli_remove == 'fine' & classification != 'fine') ~ 'zoli_fine',
    TRUE ~ classification
  )
)

# look at the numbers of each classification
d_metrics %>%
  group_by(classification) %>%
  tally() %>%
  arrange(desc(n))

filter(d_metrics, classification == 'zoli_fine')

# create a plot of all wells to check their classification, split up different strains/plasmids

# set colours for classification
cols <- RColorBrewer::brewer.pal(
  unique(d_metrics$classification) %>% length(),
  'Set1'
)
names(cols) <- unique(d_metrics$classification)

filter(d_metrics) %>%
  mutate(
    row = parse_number(well),
    id = str_extract(id, 'nal.5|wt|R1|RP4')
  ) %>%
  ggplot(
    .,
    aes(
      x = forcats::fct_reorder(as.character(row), row),
      y = final_dens,
      col = classification
    )
  ) +
  geom_beeswarm(size = 0.75) +
  facet_grid(id ~ temp) +
  theme_bw() +
  scale_color_manual(values = cols) +
  labs(y = 'Final Density', x = 'Row (antibiotic concentration')

ggsave(
  file.path('plots', 'check_final_density.pdf'),
  width = 10,
  height = 6.5
)

# add into processed data
d_processed <- d_processed %>%
  left_join(
    select(d_metrics, file, well, classification),
    by = c("file", "well")
  )

#---------------------------------------------------------#
# make a plot of each plate to check how cleaning went ####
#---------------------------------------------------------#

# create an ID column for each plate
d_processed <- d_processed %>%
  mutate(plate_id = paste(file, serial_no, sep = '_'))


plate_ids <- unique(d_processed$plate_id)

# open a pdf
pdf(
  file.path(
    'plots',
    paste('check_outliers', '.pdf', sep = '')
  ),
  width = 10,
  height = 6.5
)

for (i in 1:length(plate_ids)) {
  temp_od <- d_processed |>
    filter(plate_id == plate_ids[i]) |>
    mutate(column = str_extract(well, "[A-Z]+"), row = parse_number(well))

  temp_id <- select(temp_od, column, row, id) %>%
    distinct() %>%
    mutate(x = 0, y = 1)

  temp_plot <- temp_od %>%
    ggplot(aes(x = measurement_time_hr, y = raw_absorbance)) +
    geom_line(aes(col = classification)) +
    geom_text(aes(x = 0, y = 1, label = id), temp_id, hjust = 0) +
    facet_grid(column ~ row, switch = 'y') +
    theme_bw() +
    labs(
      title = paste('Plate ID:', plate_ids[i]),
      x = 'Time (hr)',
      y = 'Absorbance (OD600)'
    ) +
    ylim(c(0, 1.25)) +
    scale_color_manual(values = cols)

  print(temp_plot)
}

dev.off()

# remove temporary items in plot
rm(list = c('temp_od', 'temp_plot'))

# Zoltan and I manually checked any values that were zoli_failed (where zoli said they were wrong and our methods said they were fine) and zoli_fine (where zoli said they were fine but our methods said they were wrong)
filter(d_metrics, classification %in% c('zoli_fine', 'zoli_failed')) |>
  select(file, well) %>%
  write.csv(
    'data/processed/05_wells_to_check.csv',
    row.names = FALSE
  )

filter(d_metrics, file == 'DR0407_16_1')

# we combine both of those classifications and combine them
d_dan <- read.csv('data/processed/05_wells_to_check_DP.csv') %>%
  rename(classification_dan = classification)
d_zoli <- read.csv(
  'data/processed/05_wells_to_check_ZE.csv'
) %>%
  rename(classification_zoli = classification)
d_check <- left_join(d_dan, d_zoli) %>%
  mutate(
    classification = case_when(
      classification_dan == classification_zoli ~ classification_dan,
      classification_dan == 'fine' | classification_zoli == 'fine' ~ 'fine',
      TRUE ~ 'failed'
    ),
    id = paste(file, well, sep = '_')
  )

for (i in 1:nrow(d_check)) {
  d_metrics[
    d_metrics$id2 == d_check$id[i],
    'classification'
  ] <- d_check$classification[i]
  d_processed[
    d_processed$id2 == d_check$id[i],
    'classification'
  ] <- d_check$classification[i]
}

# check numbers of each classification again
d_metrics %>%
  group_by(classification) %>%
  tally() %>%
  arrange(desc(n))

filter(d_metrics, classification == 'zoli_failed')

# save this out
write.csv(
  d_metrics,
  file.path('data/processed', '05_Kanamycin_metrics_cleaned.csv'),
  row.names = FALSE
)

write.csv(
  d_processed,
  file.path(
    'data/processed',
    '05_kanamycin_processed_cleaned.csv'
  ),
  row.names = FALSE
)

# convert well number to antibiotic concentration
# filter for only wells that are fine
d_metrics <-
  mutate(
    d_metrics,
    id2 = gsub('.c', '', id),
    well_num = parse_number(well),
    conc_µl = case_when(
      well_num == 1 ~ 64,
      well_num == 2 ~ 42.7,
      well_num == 3 ~ 28.4,
      well_num == 4 ~ 19,
      well_num == 5 ~ 12.4,
      well_num == 6 ~ 8.4,
      well_num == 7 ~ 5.6,
      well_num == 8 ~ 3.8,
      well_num == 9 ~ 2.5,
      well_num == 10 ~ 1.7,
      well_num == 11 ~ 0
    )
  ) %>%
  filter(., classification == 'fine')

# save this out
write.csv(
  d_metrics,
  'data/processed/05_kanamycin_metrics_analysis.csv',
  row.names = FALSE
)

# check the growth rates. If any crazy outliers still exist, work out a way to remove them.
ggplot(aes(x = temp, y = gr), data = d_metrics) +
  geom_point() +
  theme_bw() +
  facet_wrap(~ id2 + conc_µl)
# seems to be fine.
