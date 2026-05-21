#Code to define dataframe for cal, test 1, test 2

# Libraries ----
setwd("/home/dared/GitHub/albopictus_lagged_analysis/")
library(tidyverse)
library(lubridate)
library(skimr)
library(Epi) 
library(dlnm) 
library(tsModel)
library(splines) 
library(mgcv)
library(gnm) 
library(pbs) 
library(plyr)
library(dplyr) 
library(ggplot2)
library(pracma)
library(ggpubr)
library(scales)

# Load data frame ----
bio.matrix <- readRDS("data/localities_clim_weekly.rds") %>% 
  select(ID, Country, Region, year, lat, lon, week, date, value, 
         tmean, precip) %>% 
  dplyr::rename(eggs=value, medianTweek= tmean, cumPrecweek=precip) %>% 
  dplyr::mutate(month=month(date), 
                eggs=floor(eggs), 
                year  = as.integer(year)) %>% 
  as_tibble()

# Summary stats for calibration window decision  
# Years spanning and trap count per region
region_summary <- bio.matrix %>%
  dplyr::group_by(Region, Country) %>%
  dplyr::summarise(
    n_traps  = dplyr::n_distinct(ID),
    n_years  = dplyr::n_distinct(year),
    year_min = min(year),
    year_max = max(year),
    .groups  = "drop"
  ) %>%
  dplyr::arrange(year_min)

region_summary

# Data density per region per year
region_density <- bio.matrix %>%
  dplyr::group_by(Region, year) %>%
  dplyr::summarise(
    n_traps      = dplyr::n_distinct(ID),
    n_trap_weeks = dplyr::n(),
    pct_non_zero = round(mean(eggs > 0, na.rm = TRUE) * 100, 1),
    .groups      = "drop"
  ) %>%
  dplyr::arrange(Region, year)

as.data.frame(region_summary)
as.data.frame(region_density)

# Keep only regions with >= 3 years and >= 3 traps within your window
# Calibration: 9 regions, 2018-2023
calibration_regions_list <- c(
  "Autonomous Province of Bolzano",
  "Autonomous Province of Trento",
  "Emilia-Romagna",
  "Lazio",
  "Ostschweiz",
  "Puglia",
  "Ticino",
  "Tuscany",
  "Veneto"
)

bio.matrix_cal <- bio.matrix %>%
  dplyr::filter(
    Region %in% calibration_regions_list,
    year   %in% 2018:2023
  )

# Test set 1: held-out years from the same 9 regions (2010-2017 + 2024)
bio.matrix_test_years <- bio.matrix  %>%
  dplyr::mutate(year = as.integer(year)) %>%
  dplyr::filter(
    Region %in% calibration_regions_list,
    !year  %in% 2018:2023
  )

# Test set 2: held-out regions entirely (all years)
bio.matrix_test_regions <- bio.matrix  %>%
  dplyr::mutate(year = as.integer(year)) %>%
  dplyr::filter(!Region %in% calibration_regions_list)


# Corrected sanity check
cat("Calibration  — rows:", nrow(bio.matrix_cal),
    "| traps:", dplyr::n_distinct(bio.matrix_cal$ID),
    "| regions:", dplyr::n_distinct(bio.matrix_cal$Region), "\n")

cat("Test years   — rows:", nrow(bio.matrix_test_years),
    "| traps:", dplyr::n_distinct(bio.matrix_test_years$ID),
    "| regions:", dplyr::n_distinct(bio.matrix_test_years$Region), "\n")  # <-- fixed

cat("Test regions — rows:", nrow(bio.matrix_test_regions),
    "| traps:", dplyr::n_distinct(bio.matrix_test_regions$ID),
    "| regions:", dplyr::n_distinct(bio.matrix_test_regions$Region), "\n")

# Overlap checks
stopifnot(nrow(dplyr::intersect(
  bio.matrix_cal %>% dplyr::select(ID, year, week),
  bio.matrix_test_years %>% dplyr::select(ID, year, week))) == 0)

stopifnot(nrow(dplyr::intersect(
  bio.matrix_cal %>% dplyr::select(ID, year, week),
  bio.matrix_test_regions %>% dplyr::select(ID, year, week))) == 0)

# Find the actual overlapping rows
overlap <- dplyr::intersect(
  bio.matrix_cal %>% dplyr::select(ID, year, week),
  bio.matrix_test_regions %>% dplyr::select(ID, year, week)
)

cat("Overlapping rows:", nrow(overlap), "\n")

# Which IDs are causing it?
overlap_ids <- unique(overlap$ID)
cat("Overlapping IDs:", overlap_ids, "\n")

# Check what regions these IDs appear in
bio.matrix %>%
  dplyr::filter(ID %in% overlap_ids) %>%
  dplyr::distinct(ID, Region, Country) %>%
  as.data.frame()

bio.matrix  %>%
  dplyr::filter(ID == 981) %>%
  dplyr::distinct(ID, Region, Country, year) %>%
  dplyr::arrange(year) %>%
  as.data.frame()

# DOUBLE CHECK THIS WITH MARGO, NOW WORKAROUND
bio.matrix  <- bio.matrix  %>%
  dplyr::mutate(Region = dplyr::if_else(ID == 981, "Nordwestschweiz", Region))

#rebuild everything 
bio.matrix_cal <- bio.matrix  %>%
  dplyr::filter(Region %in% calibration_regions_list, year %in% 2018:2023)

bio.matrix_test_years <- bio.matrix  %>%
  dplyr::mutate(year = as.integer(year)) %>%
  dplyr::filter(Region %in% calibration_regions_list, !year %in% 2018:2023)

bio.matrix_test_regions <- bio.matrix  %>%
  dplyr::mutate(year = as.integer(year)) %>%
  dplyr::filter(!Region %in% calibration_regions_list)

cat("Calibration  — rows:", nrow(bio.matrix_cal),
    "| traps:", dplyr::n_distinct(bio.matrix_cal$ID),
    "| regions:", dplyr::n_distinct(bio.matrix_cal$Region), "\n")

cat("Test years   — rows:", nrow(bio.matrix_test_years),
    "| traps:", dplyr::n_distinct(bio.matrix_test_years$ID),
    "| regions:", dplyr::n_distinct(bio.matrix_test_years$Region), "\n")

cat("Test regions — rows:", nrow(bio.matrix_test_regions),
    "| traps:", dplyr::n_distinct(bio.matrix_test_regions$ID),
    "| regions:", dplyr::n_distinct(bio.matrix_test_regions$Region), "\n")

stopifnot(nrow(dplyr::intersect(
  bio.matrix_cal %>% dplyr::select(ID, year, week),
  bio.matrix_test_years %>% dplyr::select(ID, year, week))) == 0)

stopifnot(nrow(dplyr::intersect(
  bio.matrix_cal %>% dplyr::select(ID, year, week),
  bio.matrix_test_regions %>% dplyr::select(ID, year, week))) == 0)

cat("No overlap confirmed\n")


#### Table 1a ----

table1a1 <- bio.matrix_cal %>%
  dplyr::group_by(Region, year) %>%
  dplyr::summarise(
    n_traps      = dplyr::n_distinct(ID),
  ) %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    n_traps_year = round(mean(n_traps),1),
  ) %>% ungroup()

table1a2 <- bio.matrix_cal %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    Country = unique(Country),
    Years = paste(sort(unique(year)), collapse = ", "),
    n_traps      = dplyr::n_distinct(ID),
    mean_pos_obs = round(100*sum(eggs > 1, na.rm = T)/sum(!is.na(eggs)), 1)
  ) %>% ungroup() 

table1a <- left_join(table1a1, table1a2) %>%
  mutate(Set= "Calibration") %>%
  select(c("Set", "Country", "Region", "Years", "n_traps", "n_traps_year", "mean_pos_obs")) %>%
  arrange(-n_traps_year)

#### Table 1b ----

table1b1 <- bio.matrix_test_years %>%
  dplyr::group_by(Region, year) %>%
  dplyr::summarise(
    n_traps      = dplyr::n_distinct(ID),
  ) %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    n_traps_year = round(mean(n_traps),1),
  ) %>% ungroup()

table1b2 <- bio.matrix_test_years %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    Country = unique(Country),
    Years = paste(sort(unique(year)), collapse = ", "),
    n_traps      = dplyr::n_distinct(ID),
    mean_pos_obs = round(100*sum(eggs > 1, na.rm = T)/sum(!is.na(eggs)), 1)
  ) %>% ungroup() 

table1b <- left_join(table1b1, table1b2) %>%
  mutate(Set= "Test 1 (temporal)") %>%
  select(c("Set", "Country", "Region", "Years", "n_traps", "n_traps_year", "mean_pos_obs")) %>%
  arrange(-n_traps_year)

#### Table 1c ----

table1c1 <- bio.matrix_test_regions %>%
  dplyr::group_by(Region, year) %>%
  dplyr::summarise(
    n_traps      = dplyr::n_distinct(ID),
  ) %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    n_traps_year = round(mean(n_traps),1),
  ) %>% ungroup()

table1c2 <- bio.matrix_test_regions %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    Country = unique(Country),
    Years = paste(sort(unique(year)), collapse = ", "),
    n_traps      = dplyr::n_distinct(ID),
    mean_pos_obs = round(100*sum(eggs > 1, na.rm = T)/sum(!is.na(eggs)), 1)
  ) %>% ungroup() 

table1c <- left_join(table1c1, table1c2) %>%
  mutate(Set= "Test 2 (spatial)") %>%
  select(c("Set", "Country", "Region", "Years", "n_traps", "n_traps_year", "mean_pos_obs")) %>%
  arrange(-n_traps_year)

# write table

table1 = rbind(table1a, table1b, table1c)

write.csv(table1, "outputs/Table1.csv", row.names = FALSE)
print(table1)

# EXPORT DATASETS
saveRDS(bio.matrix_sel,          "data/bio.matrix_cal.rds")
saveRDS(bio.matrix_test_years,   "data/bio.matrix_test_years.rds")
saveRDS(bio.matrix_test_regions, "data/bio.matrix_test_regions.rds")

cat("Exported:\n")
cat("  Calibration : ", nrow(bio.matrix_sel),          "rows\n")
cat("  Test years  : ", nrow(bio.matrix_test_years),   "rows\n")
cat("  Test regions: ", nrow(bio.matrix_test_regions), "rows\n")