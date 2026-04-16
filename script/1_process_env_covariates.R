setwd("/home/dared/GitHub/albopictus_lagged_analysis/")
library(terra)
library(tidyverse)
library(lubridate)
library(sf)
library(mapview) 

# 0. HELPER FUNCTIONS ----
 parse_layer_times <- function(rast_obj) {
  as_date(
    as_datetime(
      as.numeric(sub(".*=", "", names(rast_obj))),
      origin = "1970-01-01", tz = "UTC"
    )
  )
}

build_layer_ref <- function(rast_obj, vars,
                            years = 2010:2024, weeks = 1:53) {  # <-- 1:53
  all_dates <- parse_layer_times(rast_obj)
  nms       <- names(rast_obj)
  
  dates_ref  <- all_dates[which(grepl(vars[1], nms))]
  
  layer_cols <- map(vars, ~ which(grepl(.x, nms))) %>%
    set_names(paste0("layer_", names(vars))) %>%
    as_tibble()
  
  bind_cols(tibble(date = dates_ref), layer_cols) %>%
    mutate(
      year = isoyear(date),
      week = isoweek(date)
    ) %>%
    filter(year %in% years, week %in% weeks)
}

extract_climate <- function(rast_obj, layer_ref, varnames, localities_df, loc_vect,
                            id_col = "ID") {
  map(varnames, function(v) {
    lyr_col <- paste0("layer_", v)
    layers  <- layer_ref[[lyr_col]]
    
    terra::extract(rast_obj[[layers]], loc_vect) %>%
      select(-ID) %>%
      mutate(loc_id = localities_df[[id_col]]) %>%
      pivot_longer(-loc_id,
                   names_to  = "col_idx",
                   values_to = v) %>%
      group_by(loc_id) %>%
      mutate(row = row_number()) %>%
      ungroup() %>%
      left_join(
        layer_ref %>%
          mutate(row = row_number()) %>%
          select(row, year, week, date),
        by = "row"
      ) %>%
      select(loc_id, year, week, date, all_of(v)) %>%
      rename(!!id_col := loc_id)
    
  }) %>%
    reduce(left_join, by = c(id_col, "year", "week", "date"))
}

 
# 1. LOAD OBSERVATIONS ----
localities <- readRDS("data/aggregated_albopictus_observations_2026-03-23.RDS")

unique_locs <- localities %>%
  ungroup() %>%
  distinct(ID, lon, lat)

cat("Unique locations:", nrow(unique_locs), "\n")
 
obs_years <- localities %>%
  mutate(year_int = as.integer(year)) %>%
  pull(year_int) %>%
  unique() %>%
  sort()

cat("Years in observations:", paste(obs_years, collapse = ", "), "\n")

localities %>%
  ungroup() %>%
  mutate(year_int = as.integer(year)) %>%
  group_by(ID) %>%
  summarise(
    n_years   = n_distinct(year_int),
    year_min  = min(year_int),
    year_max  = max(year_int),
    span      = year_max - year_min + 1,      # total window length
    gaps      = span - n_years                 # years missing within window
  ) %>%
  arrange(desc(n_years), year_min) %>%
  print(n = Inf)

localities %>%
  ungroup() %>%
  mutate(year_int = as.integer(year)) %>%
  filter(year_int %in% 2018:2023) %>%
  group_by(ID) %>%
  summarise(
    n_years  = n_distinct(year_int),
    year_min = min(year_int),
    year_max = max(year_int),
    span     = year_max - year_min + 1,
    gaps     = span - n_years
  ) %>%
  arrange(desc(n_years)) %>%
  count(n_years, name = "n_traps") %>%
  print(n = Inf)


localities %>%
  ungroup() %>%
  distinct(ID, year) %>%
  count(ID, name = "n_years") %>%
  ggplot(aes(x = n_years)) +
  geom_histogram(binwidth = 1, fill = "steelblue", colour = "white") +
  scale_x_continuous(breaks = 1:15) +
  labs(title = "Distribution of sampling duration per trap",
       x = "Number of years sampled", y = "Number of traps") +
  theme_minimal()


loc_vect <- vect(unique_locs,
                 geom = c("lon", "lat"),
                 crs  = "EPSG:4326")
mapview(st_as_sf(loc_vect))

# 2. LOAD RASTERS ----
indir <- "/home/dared/Documents/PoD/marieCurie2022/ZanZemap/ML/Download_Process_ClimData/2025/data/weekly/"

tas <- rast(paste0(indir, "t2m_weekly_stats_201001-202512_C.nc"))
tp  <- rast(paste0(indir, "tp_weekly_cumulative_201001-202512.nc"))

# ── Coastal gap filling ───────────────────────────────────────────────────────
# 5×5 focal mean applied only to NA cells (coastal/border locations)
# pad = TRUE ensures edge cells are handled correctly
w   <- matrix(1, 5, 5)
tas <- focal(tas, w, fun = mean, na.rm = TRUE, na.policy = "only", pad = TRUE)
tp  <- focal(tp,  w, fun = mean, na.rm = TRUE, na.policy = "only", pad = TRUE)

cat("tas layers:", nlyr(tas), "\n")
cat("tp  layers:", nlyr(tp),  "\n")

# 3. BUILD LAYER REFERENCES ----
tas_vars <- c(tmean = "t2m_mean", tmin = "t2m_min", tmax = "t2m_max")
tp_vars  <- c(precip = "tp")

layer_ref_tas <- build_layer_ref(tas, tas_vars, years = obs_years)
layer_ref_tp  <- build_layer_ref(tp,  tp_vars,  years = obs_years)

cat("tas layer ref rows:", nrow(layer_ref_tas), "\n")
cat("tp  layer ref rows:", nrow(layer_ref_tp),  "\n")

# 4. EXTRACT ----
clim_temp   <- extract_climate(tas, layer_ref_tas, names(tas_vars),
                               unique_locs, loc_vect, id_col = "ID")

clim_precip <- extract_climate(tp,  layer_ref_tp,  names(tp_vars),
                               unique_locs, loc_vect, id_col = "ID")

clim <- left_join(clim_temp, clim_precip,
                  by = c("ID", "year", "week", "date"))

cat("Climate rows:", nrow(clim), "\n")
# Expect: nrow(unique_locs) * nrow(layer_ref_tas)

# 5. JOIN CLIMATE BACK TO FULL OBSERVATIONS -----
localities_clim <- localities %>%
  mutate(year_int = as.integer(year),
         week_int = as.integer(week)) %>%
  left_join(
    clim %>% rename(year_int = year, week_int = week),
    by = c("ID", "year_int", "week_int")
  ) %>%
  select(-year_int, -week_int)

cat("Final rows:", nrow(localities_clim), "\n")   # must be 59,734
cat("Final cols:", ncol(localities_clim), "\n")

# 6. SANITY CHECKS ----
localities_clim %>%
  summarise(across(c(tmean, tmin, tmax, precip), ~ sum(is.na(.))))

localities_clim %>%
  ggplot(aes(x = as.integer(week), y = tmean, col = year)) +
  geom_line(alpha = 0.3) +
  facet_wrap(~ Region) +
  labs(title = "Weekly mean temperature by region", y = "tmean (°C)") +
  theme_minimal()

localities_clim %>%
  ggplot(aes(x = as.integer(week), y = precip, col = year)) +
  geom_bar(alpha = 0.3, stat="identity") +
  facet_grid( year ~ Region) +
  labs(title = "Weekly cumulative precipitation by region", y = "Cumulative precipitation (mm)") +
  theme_minimal()

# 7. SAVE ----
str(localities_clim)
saveRDS(localities_clim, "data/localities_clim_weekly.rds")

cat("Saved:", nrow(localities_clim), "rows,", ncol(localities_clim), "columns\n")
cat("Unique IDs:", n_distinct(localities_clim$ID), "\n")
