#Code to run a model selection based by a cross-validation upon non linear lagged model to study the relation between weather and oviposition

# Libraries ----
setwd("/home/ddare/GitHub/albopictus_lagged_analysis/")
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

bio.matrix <- bio.matrix %>%
  dplyr::mutate(
    Country = dplyr::if_else(Region == "Lombardia", "CH", Country),
    Region  = dplyr::if_else(Region == "Lombardia", "Ticino", Region)
  )
# no trap-week should be duplicated by the merge
bio.matrix %>% dplyr::count(ID, year, week) %>% dplyr::filter(n > 1)   # expect 0 rows

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
  "Provincia Autonoma di Bolzano/Bozen",
  "Provincia Autonoma di Trento",
  "Emilia-Romagna",
  "Lazio",
  "Puglia",
  "Ticino",
  "Toscana",
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

#traps subset
bio.matrix_test_years %>%
  dplyr::filter(ID %in% unique(bio.matrix_cal$ID)) %>%
  dplyr::summarise(rows = dplyr::n(), traps = dplyr::n_distinct(ID))

bio.matrix_test_years %>%
  dplyr::filter(!ID %in% unique(bio.matrix_cal$ID)) %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(traps = dplyr::n_distinct(ID), first_year = min(year), last_year = max(year))

dplyr::n_distinct(bio.matrix_cal$ID)                  # [N] calibration traps

bio.matrix_test_years %>%
  dplyr::summarise(rows = dplyr::n(), traps = dplyr::n_distinct(ID))   # temporal test set

bio.matrix_test_years %>%
  dplyr::filter(ID %in% unique(bio.matrix_cal$ID)) %>%
  dplyr::summarise(rows = dplyr::n(), traps = dplyr::n_distinct(ID))   # known traps, [M]


### add nweek and weekn
bio.matrix_sel <- bio.matrix_cal%>%
  dplyr::arrange(ID, year, week) %>% 
  dplyr::group_by(ID) %>%
  dplyr::mutate(nweek = 1:n()) %>%
  dplyr::mutate(weekn = as.numeric(week)) %>%
  dplyr::ungroup()

# Confirm nweek is consecutive and starts at 1 per trap
bio.matrix_sel %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(
    nweek_min  = min(nweek),
    nweek_max  = max(nweek),
    n_rows     = dplyr::n(),
    is_consecutive = (nweek_max - nweek_min + 1) == n_rows,
    .groups = "drop"
  ) %>%
  dplyr::summarise(
    all_start_at_1     = all(nweek_min == 1),
    all_consecutive    = all(is_consecutive),
    min_nweek_max      = min(nweek_max),
    max_nweek_max      = max(nweek_max)
  ) %>%
  as.data.frame()

# print tab 3 
tab_data <- dplyr::bind_rows(
  bio.matrix_cal          %>% dplyr::mutate(Set = "Calibration"),
  bio.matrix_test_years   %>% dplyr::mutate(Set = "Temporal test"),
  bio.matrix_test_regions %>% dplyr::mutate(Set = "Spatial test")
) %>%
  dplyr::mutate(
    Country = dplyr::case_when(
      Region %in% c("Freiburg", "Liechtenstein") ~ "Switzerland",
      Country == "IT" ~ "Italy",
      Country == "CH" ~ "Switzerland",
      Country == "FR" ~ "France",
      Country == "AL" ~ "Albania"
    ),
    Region = dplyr::recode(Region,
                           "Provincia Autonoma di Bolzano/Bozen" = "Aut. P. Bolzano",
                           "Provincia Autonoma di Trento"        = "Aut. P. Trento",
                           "Toscana"                             = "Tuscany",
                           "Sicilia"                             = "Sicily",
                           "Provence-Alpes-Côte d’Azur"          = "Côte d’Azur",
                           "Région lémanique"                    = "Région Lémanique")
  )

# years: runs of 3+ consecutive years as a range, shorter runs listed individually
years_tab <- tab_data %>%
  dplyr::distinct(Set, Country, Region, year) %>%
  dplyr::arrange(Set, Region, year) %>%
  dplyr::group_by(Set, Country, Region) %>%
  dplyr::mutate(run = cumsum(c(TRUE, diff(year) > 1))) %>%
  dplyr::group_by(Set, Country, Region, run) %>%
  dplyr::summarise(y_min = min(year), y_max = max(year), n_y = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(label = dplyr::case_when(
    n_y >= 3 ~ paste0(y_min, "-", y_max),
    n_y == 2 ~ paste0(y_min, ", ", y_max),
    TRUE     ~ as.character(y_min))) %>%
  dplyr::group_by(Set, Country, Region) %>%
  dplyr::summarise(Years = paste(label, collapse = ", "), .groups = "drop")

# traps and detection rate
metrics_tab <- tab_data %>%
  dplyr::group_by(Set, Country, Region) %>%
  dplyr::summarise(
    N_traps      = dplyr::n_distinct(ID),
    pct_non_zero = round(mean(eggs > 0, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  )

table1 <- dplyr::left_join(years_tab, metrics_tab, by = c("Set", "Country", "Region")) %>%
  dplyr::mutate(Set = factor(Set, levels = c("Calibration", "Temporal test", "Spatial test"))) %>%
  dplyr::arrange(Set, Country, dplyr::desc(N_traps))

as.data.frame(table1)
write.csv(table1, "outputs/Table3_dataset_summary.csv", row.names = FALSE)

# Define cross-valitation settings ----
n_ids <- nrow(bio.matrix_sel)
reps <- 20
prop_fit <- 0.7

#settings lag
tlag <- 8

# add n row + and new column for stratified sampling
bio.matrix_sel <- bio.matrix_sel%>%
  dplyr::mutate(IDdata = 1:n_ids) %>%
  dplyr::mutate(IDy = paste0(ID, "_", year))

# Performances df ----
#we define a data.frame with 10 rows, one for random repetitions, and n+1 columns (n = models)

# RMSE: root mean square error
RMSE_val_df <- data.frame(rep = 1:reps,
                          seas = NA,
                          seas.trend = NA,
                          Temp = NA,
                          T.seas = NA,
                          T.P.seas = NA,
                          T.seas.trend = NA,
                          T.P.seas.trend = NA)

RMSE_fit_df <- RMSE_val_df

# MAE: mean absolute error
MAE_val_df <- RMSE_val_df
MAE_fit_df <- RMSE_val_df

# RMSLE: root mean square logarithmic error
RMSLE_val_df <- RMSE_val_df
RMSLE_fit_df <- RMSE_val_df

# Residual deviance
ResDev_df <- RMSE_val_df

# Overdispersion
Dispersion_df <-RMSE_val_df

# Cycle to fit and validate the model ----
set.seed(123)   # reproducible splits
tic()
for(i in 1:reps){
  
  # for stratified sampling:
  
  # I must leave the same IDs and the same years in both subsets
  # I should not remove weather from any subsets (but only observations)
  
  # extract random batch from fit and validation  ----
  # ids_i = bio.matrix_sel %>%
  #   dplyr::group_by(IDy)%>%
  #   sample_n(size = 52*prop_fit) %>% 
  #   pull(IDdata)
  ids_i <- bio.matrix_sel %>%
    dplyr::group_by(IDy) %>%
    dplyr::filter(n() >= 2) %>%          # safety: skip trivially small groups
    dplyr::sample_frac(prop_fit) %>%
    dplyr::pull(IDdata)
  
  bio.matrix_fit_i = bio.matrix_sel
  bio.matrix_fit_i$eggs[-ids_i] = NA
  
  bio.matrix_val_i = bio.matrix_sel
  bio.matrix_val_i$eggs[ids_i] = NA
  
  # Extract obs_val for validation
  obs_val = bio.matrix_val_i$eggs
  obs_val_na = which(is.na(obs_val))
  
  # 1. Seas and year ----
  
  # natural spline over seasonality (df = 4)
  bseas <- ns(bio.matrix_fit_i$weekn, df=4)
  
  # predictive formulation 
  bseas_df <- as.data.frame(bseas)
  colnames(bseas_df) <- paste0("bs", seq_len(ncol(bseas_df)))
  
  # and cbind
  bio.matrix_fit_i_ext_seas <- cbind(bio.matrix_fit_i, bseas_df) 
  
  # Build formula with those columns
  form <- as.formula(paste("eggs ~",
                           paste(c(colnames(bseas_df), "year"), collapse = " + ")))
  
  # model definition 
  # with only spline trend & seasonality
  mod_seas_i <- gnm(form,
                    eliminate=as.factor(ID), #grouping factor
                    data=bio.matrix_fit_i_ext_seas,
                    family=quasipoisson(link="log"),
                    na.action= "na.exclude")
  
  # assess fitting error
  
  pred_seas_i = mod_seas_i$fitted
  obs_fit_i = mod_seas_i$model$eggs   # observed counts (fitted + residuals is not y: residuals are working residuals)
  stopifnot(length(pred_seas_i) == length(obs_fit_i))
  
  RMSE_fit_df$seas[i] = sqrt(mean((pred_seas_i - obs_fit_i)^2, na.rm = T))
  RMSLE_fit_df$seas[i]  = sqrt(mean((log(pred_seas_i+1) - log(obs_fit_i+1))^2, na.rm = T))
  MAE_fit_df$seas[i]  = mean(abs(pred_seas_i - obs_fit_i), na.rm = T)
  
  ## assess model over validation sample ----
  
  # natural spline over trend 
  bseas_val <- ns(bio.matrix_val_i$weekn, df=4)
  
  bseas_val_df <- as.data.frame(bseas_val)
  colnames(bseas_val_df) <- paste0("bs", seq_len(ncol(bseas_val_df)))
  
  bio.matrix_val_i_seas <- cbind(bio.matrix_val_i, bseas_val_df)
  
  # predict
  pred_seas <- predict(mod_seas_i, newdata = bio.matrix_val_i_seas, type = "response")
  
  RMSE_val_df$seas[i] = sqrt(mean((pred_seas[-(obs_val_na)] - obs_val[-(obs_val_na)])^2, na.rm = T))
  RMSLE_val_df$seas[i]  = sqrt(mean((log(pred_seas[-(obs_val_na)]+1) - log(obs_val[-(obs_val_na)]+1))^2, na.rm = T))
  MAE_val_df$seas[i]  = mean(abs(pred_seas[-(obs_val_na)] - obs_val[-(obs_val_na)]), na.rm = T)
  
  # 2. Seas.trend ----
  
  # df trend
  dftrend <- 4 * length(unique(bio.matrix_fit_i$year))
  
  # natural spline over trend
  btrend <- ns(bio.matrix_fit_i$nweek, knots=equalknots(bio.matrix_fit_i$nweek, dftrend-1))
  
  # predictive formulation 
  btrend_df <- as.data.frame(btrend)
  colnames(btrend_df) <- paste0("bt", seq_len(ncol(btrend_df)))
  
  # and cbind
  bio.matrix_fit_i_ext_seas.trend <- cbind(bio.matrix_fit_i, btrend_df) 
  
  # Build formula with those columns
  form <- as.formula(paste("eggs ~",
                           paste(colnames(btrend_df), collapse = " + ")))
  
  # model definition 
  # with only spline trend & seasonality
  mod_seas.trend_i <- gnm(form,
                          eliminate=as.factor(ID), #grouping factor
                          data=bio.matrix_fit_i_ext_seas.trend,
                          family=quasipoisson(link="log"),
                          na.action= "na.exclude")
  
  # assess fitting error
  
  pred_seas.trend_i = mod_seas.trend_i$fitted
  obs_fit_i = mod_seas.trend_i$model$eggs   # observed counts (fitted + residuals is not y: residuals are working residuals)
  stopifnot(length(pred_seas.trend_i) == length(obs_fit_i))
  
  RMSE_fit_df$seas.trend[i] = sqrt(mean((pred_seas.trend_i - obs_fit_i)^2, na.rm = T))
  RMSLE_fit_df$seas.trend[i]  = sqrt(mean((log(pred_seas.trend_i+1) - log(obs_fit_i+1))^2, na.rm = T))
  MAE_fit_df$seas.trend[i]  = mean(abs(pred_seas.trend_i - obs_fit_i), na.rm = T)
  
  ## assess model over validation sample ----
  
  # natural spline over trend (same dftrend as before)
  btrend_val <- ns(bio.matrix_val_i$nweek, knots=equalknots(bio.matrix_val_i$nweek, dftrend-1))
  
  btrend_val_df <- as.data.frame(btrend_val)
  colnames(btrend_val_df) <- paste0("bt", seq_len(ncol(btrend_val_df)))
  
  bio.matrix_val_i_seas.trend <- cbind(bio.matrix_val_i, btrend_val_df)
  
  # predict
  pred_seas.trend <- predict(mod_seas.trend_i, newdata = bio.matrix_val_i_seas.trend, type = "response")
  
  RMSE_val_df$seas.trend[i] = sqrt(mean((pred_seas.trend[-(obs_val_na)] - obs_val[-(obs_val_na)])^2, na.rm = T))
  RMSLE_val_df$seas.trend[i]  = sqrt(mean((log(pred_seas.trend[-(obs_val_na)]+1) - log(obs_val[-(obs_val_na)]+1))^2, na.rm = T))
  MAE_val_df$seas.trend[i]  = mean(abs(pred_seas.trend[-(obs_val_na)] - obs_val[-(obs_val_na)]), na.rm = T)
  
  # 3 Temp  ----
  
  # knots for temperature
  tknots<-unname(quantile(bio.matrix_fit_i$medianTweek, probs = c(0.25, 0.75), na.rm = TRUE))
  
  # cross basis over temperature
  cbTemp <- crossbasis(bio.matrix_fit_i$medianTweek, lag=tlag,
                       argvar=list(fun="bs", degree=2, knots=tknots),
                       arglag=list(knots=c(1,4)),
                       group=bio.matrix_fit_i$Region)
  
  # define model (preditive formulation)
  cbTemp_df <- as.data.frame(cbTemp)
  colnames(cbTemp_df) <- paste0("cbt", seq_len(ncol(cbTemp_df)))
  
  # and cbind
  bio.matrix_fit_i_ext_T <- cbind(bio.matrix_fit_i, cbTemp_df) 
  
  # Build formula with those columns
  form <- as.formula(paste("eggs ~",
                           paste(colnames(cbTemp_df), collapse = " + ")))
  
  # model definition 
  # with temperature, spline trend & seasonality
  mod_T_i <- gnm(form,
                 eliminate=as.factor(ID), #grouping factor
                 data=bio.matrix_fit_i_ext_T,
                 family=quasipoisson(link="log"),
                 na.action= "na.exclude")
  
  # assess fitting error
  
  pred_T_i = mod_T_i$fitted
  obs_fit_i = mod_T_i$model$eggs   # observed counts (fitted + residuals is not y: residuals are working residuals)
  stopifnot(length(pred_T_i) == length(obs_fit_i))
  
  RMSE_fit_df$Temp[i] = sqrt(mean((pred_T_i - obs_fit_i)^2, na.rm = T))
  RMSLE_fit_df$Temp[i]  = sqrt(mean((log(pred_T_i+1) - log(obs_fit_i+1))^2, na.rm = T))
  MAE_fit_df$Temp[i]  = mean(abs(pred_T_i - obs_fit_i), na.rm = T)
  
  ## assess model over validation sample ----
  
  # Build new crossbasis
  cbTemp_val <- crossbasis(bio.matrix_val_i$medianTweek, lag=tlag,
                           argvar=list(fun="bs", degree=2, knots=tknots),
                           arglag=list(knots=c(1,4)),
                           group=bio.matrix_val_i$Region)
  
  cbTemp_val_df <- as.data.frame(cbTemp_val)
  colnames(cbTemp_val_df) <- paste0("cbt", seq_len(ncol(cbTemp_val_df)))
  
  bio.matrix_val_i_T <- cbind(bio.matrix_val_i, cbTemp_val_df)
  
  # predict
  pred_T<- predict(mod_T_i, newdata = bio.matrix_val_i_T, type = "response")
  
  RMSE_val_df$Temp[i] = sqrt(mean((pred_T[-(obs_val_na)] - obs_val[-(obs_val_na)])^2, na.rm = T))
  RMSLE_val_df$Temp[i]  = sqrt(mean((log(pred_T[-(obs_val_na)]+1) - log(obs_val[-(obs_val_na)]+1))^2, na.rm = T))
  MAE_val_df$Temp[i]  = mean(abs(pred_T[-(obs_val_na)] - obs_val[-(obs_val_na)]), na.rm = T)
  
  # 4 T.seas  and year ----
  
  # cbind
  bio.matrix_fit_i_ext_T.seas <- cbind(bio.matrix_fit_i, cbTemp_df, bseas_df) 
  
  # Build formula with those columns
  form <- as.formula(paste("eggs ~",
                           paste(c(colnames(cbind(cbTemp_df, bseas_df)), "year"), collapse = " + ")))
  
  # model definition 
  # with temperature, spline seasonality
  mod_T.seas_i <- gnm(form,
                      eliminate=as.factor(ID), #grouping factor
                      data=bio.matrix_fit_i_ext_T.seas,
                      family=quasipoisson(link="log"),
                      na.action= "na.exclude")
  
  # assess fitting error
  
  pred_T.seas_i = mod_T.seas_i$fitted
  obs_fit_i = mod_T.seas_i$model$eggs   # observed counts (fitted + residuals is not y: residuals are working residuals)
  stopifnot(length(pred_T.seas_i) == length(obs_fit_i))
  
  RMSE_fit_df$T.seas[i] = sqrt(mean((pred_T.seas_i - obs_fit_i)^2, na.rm = T))
  RMSLE_fit_df$T.seas[i]  = sqrt(mean((log(pred_T.seas_i+1) - log(obs_fit_i+1))^2, na.rm = T))
  MAE_fit_df$T.seas[i]  = mean(abs(pred_T.seas_i - obs_fit_i), na.rm = T)
  
  ## assess model over validation sample ----
  
  bio.matrix_val_i_T.seas <- cbind(bio.matrix_val_i, cbTemp_val_df, bseas_val_df)
  
  # predict
  pred_T.seas<- predict(mod_T.seas_i, newdata = bio.matrix_val_i_T.seas, type = "response")
  
  RMSE_val_df$T.seas[i] = sqrt(mean((pred_T.seas[-(obs_val_na)] - obs_val[-(obs_val_na)])^2, na.rm = T))
  RMSLE_val_df$T.seas[i]  = sqrt(mean((log(pred_T.seas[-(obs_val_na)]+1) - log(obs_val[-(obs_val_na)]+1))^2, na.rm = T))
  MAE_val_df$T.seas[i]  = mean(abs(pred_T.seas[-(obs_val_na)] - obs_val[-(obs_val_na)]), na.rm = T)
  
  # 5. T.P.seas and year ----
  
  # knots for precipitation
  pknots<-unname(quantile(bio.matrix_fit_i$cumPrecweek, probs = c(0.25, 0.75, 0.9), na.rm = TRUE))
  
  cbPrec <- crossbasis(bio.matrix_fit_i$cumPrecweek, lag=tlag,
                       argvar=list(fun="bs", degree=2, knots=pknots),
                       arglag=list(knots=c(1,4)),
                       group=bio.matrix_fit_i$Region)
  
  # define model (preditive formulation)
  cbPrec_df <- as.data.frame(cbPrec)
  colnames(cbPrec_df) <- paste0("cbp", seq_len(ncol(cbPrec_df)))
  
  # cbind
  bio.matrix_fit_i_ext_T.P.seas <- cbind(bio.matrix_fit_i, cbTemp_df, cbPrec_df, bseas_df) 
  
  # Build formula with those columns
  form <- as.formula(paste("eggs ~",
                           paste(c(colnames(cbind(cbTemp_df, cbPrec_df,bseas_df)), "year"), collapse = " + ")))
  
  # model definition 
  # with temperature, spline seasonality
  mod_T.P.seas_i <- gnm(form,
                        eliminate=as.factor(ID), #grouping factor
                        data=bio.matrix_fit_i_ext_T.P.seas,
                        family=quasipoisson(link="log"),
                        na.action= "na.exclude")
  
  # assess fitting error
  
  pred_T.P.seas_i = mod_T.P.seas_i$fitted
  obs_fit_i = mod_T.P.seas_i$model$eggs   # observed counts (fitted + residuals is not y: residuals are working residuals)
  stopifnot(length(pred_T.P.seas_i) == length(obs_fit_i))
  
  RMSE_fit_df$T.P.seas[i] = sqrt(mean((pred_T.P.seas_i - obs_fit_i)^2, na.rm = T))
  RMSLE_fit_df$T.P.seas[i]  = sqrt(mean((log(pred_T.P.seas_i+1) - log(obs_fit_i+1))^2, na.rm = T))
  MAE_fit_df$T.P.seas[i]  = mean(abs(pred_T.P.seas_i - obs_fit_i), na.rm = T)
  
  ## assess model over validation sample ----
  
  # Build new crossbasis for Prec
  cbPrec_val <- crossbasis(bio.matrix_val_i$cumPrecweek, lag=tlag,
                           argvar=list(fun="bs", degree=2, knots=pknots),
                           arglag=list(knots=c(1,4)),
                           group=bio.matrix_val_i$Region)
  
  cbPrec_val_df <- as.data.frame(cbPrec_val)
  colnames(cbPrec_val_df) <- paste0("cbp", seq_len(ncol(cbPrec_val_df)))
  
  bio.matrix_val_i_T.P.seas <- cbind(bio.matrix_val_i, cbTemp_val_df, cbPrec_val_df, bseas_val_df)
  
  # predict
  pred_T.P.seas<- predict(mod_T.P.seas_i, newdata = bio.matrix_val_i_T.P.seas, type = "response")
  
  RMSE_val_df$T.P.seas[i] = sqrt(mean((pred_T.P.seas[-(obs_val_na)] - obs_val[-(obs_val_na)])^2, na.rm = T))
  RMSLE_val_df$T.P.seas[i]  = sqrt(mean((log(pred_T.P.seas[-(obs_val_na)]+1) - log(obs_val[-(obs_val_na)]+1))^2, na.rm = T))
  MAE_val_df$T.P.seas[i]  = mean(abs(pred_T.P.seas[-(obs_val_na)] - obs_val[-(obs_val_na)]), na.rm = T)
  
  # 6. T.seas.trend ----
  
  # and cbind
  bio.matrix_fit_i_ext_T.seas.trend <- cbind(bio.matrix_fit_i, cbTemp_df, btrend_df) 
  
  # Build formula with those columns
  form <- as.formula(paste("eggs ~",
                           paste(colnames(cbind(cbTemp_df, btrend_df)), collapse = " + ")))
  
  
  # model definition 
  # with temperature, spline trend & seasonality
  mod_T.seas.trend_i <- gnm(form,
                            eliminate=as.factor(ID), #grouping factor
                            data=bio.matrix_fit_i_ext_T.seas.trend,
                            family=quasipoisson(link="log"),
                            na.action= "na.exclude")
  
  # assess fitting error
  
  pred_mod_T.seas.trend_i = mod_T.seas.trend_i$fitted
  obs_fit_i = mod_T.seas.trend_i$model$eggs   # observed counts (fitted + residuals is not y: residuals are working residuals)
  stopifnot(length(pred_mod_T.seas.trend_i) == length(obs_fit_i))
  
  RMSE_fit_df$T.seas.trend[i] = sqrt(mean((pred_mod_T.seas.trend_i - obs_fit_i)^2, na.rm = T))
  RMSLE_fit_df$T.seas.trend[i]  = sqrt(mean((log(pred_mod_T.seas.trend_i+1) - log(obs_fit_i+1))^2, na.rm = T))
  MAE_fit_df$T.seas.trend[i]  = mean(abs(pred_mod_T.seas.trend_i - obs_fit_i), na.rm = T)
  
  ## assess model over validation sample ----
  
  bio.matrix_val_i_T.seas.trend <- cbind(bio.matrix_val_i, cbTemp_val_df, btrend_val_df)
  
  # predict
  pred_T.seas.trend <- predict(mod_T.seas.trend_i, newdata = bio.matrix_val_i_T.seas.trend, type = "response")
  
  RMSE_val_df$T.seas.trend[i] = sqrt(mean((pred_T.seas.trend[-(obs_val_na)] - obs_val[-(obs_val_na)])^2, na.rm = T))
  RMSLE_val_df$T.seas.trend[i]  = sqrt(mean((log(pred_T.seas.trend[-(obs_val_na)]+1) - log(obs_val[-(obs_val_na)]+1))^2, na.rm = T))
  MAE_val_df$T.seas.trend[i]  = mean(abs(pred_T.seas.trend[-(obs_val_na)] - obs_val[-(obs_val_na)]), na.rm = T)
  
  # 7. T.P.seas.trend ----
  
  # define cbind
  bio.matrix_fit_i_ext_T.P.seas.trend <- cbind(bio.matrix_fit_i, cbTemp_df, cbPrec_df, btrend_df) 
  
  # Build formula with those columns
  form <- as.formula(paste("eggs ~",
                           paste(colnames(cbind(cbTemp_df, cbPrec_df, btrend_df)), collapse = " + ")))
  
  # model definition 
  # with temperature, precipitaiton, with spline trend & seasonality
  mod_T.P.seas.trend_i <- gnm(form,
                              eliminate=as.factor(ID), #grouping factor
                              data=bio.matrix_fit_i_ext_T.P.seas.trend,
                              family=quasipoisson(link="log"),
                              na.action= "na.exclude")
  
  # assess fitting error
  
  pred_mod_T.P.seas.trend_i = mod_T.P.seas.trend_i$fitted
  obs_fit_i = mod_T.P.seas.trend_i$model$eggs   # observed counts (fitted + residuals is not y: residuals are working residuals)
  stopifnot(length(pred_mod_T.P.seas.trend_i) == length(obs_fit_i))
  
  RMSE_fit_df$T.P.seas.trend[i] = sqrt(mean((pred_mod_T.P.seas.trend_i - obs_fit_i)^2, na.rm = T))
  RMSLE_fit_df$T.P.seas.trend[i]  = sqrt(mean((log(pred_mod_T.P.seas.trend_i+1) - log(obs_fit_i+1))^2, na.rm = T))
  MAE_fit_df$T.P.seas.trend[i]  = mean(abs(pred_mod_T.P.seas.trend_i - obs_fit_i), na.rm = T)
  
  ## assess model over validation sample ----
  
  bio.matrix_val_i_T.P.seas.trend <- cbind(bio.matrix_val_i, cbTemp_val_df, cbPrec_val_df, btrend_val_df)
  
  # predict
  pred_T.P.seas.trend <- predict(mod_T.P.seas.trend_i, newdata = bio.matrix_val_i_T.P.seas.trend, type = "response")
  
  RMSE_val_df$T.P.seas.trend[i] = sqrt(mean((pred_T.P.seas.trend[-(obs_val_na)] - obs_val[-(obs_val_na)])^2, na.rm = T))
  RMSLE_val_df$T.P.seas.trend[i]  = sqrt(mean((log(pred_T.P.seas.trend[-(obs_val_na)]+1) - log(obs_val[-(obs_val_na)]+1))^2, na.rm = T))
  MAE_val_df$T.P.seas.trend[i]  = mean(abs(pred_T.P.seas.trend[-(obs_val_na)] - obs_val[-(obs_val_na)]), na.rm = T)
  
  
  ## Other metrics ----
  # Overdispersion parameters and residual deviance
  
  Dispersion_df$seas[i] = summary(mod_seas_i)$dispersion
  Dispersion_df$seas.trend[i] = summary(mod_seas.trend_i)$dispersion
  Dispersion_df$Temp[i] = summary(mod_T_i)$dispersion
  Dispersion_df$T.seas[i] = summary(mod_T.seas_i)$dispersion
  Dispersion_df$T.P.seas[i] = summary(mod_T.P.seas_i)$dispersion
  Dispersion_df$T.seas.trend[i] = summary(mod_T.seas.trend_i)$dispersion
  Dispersion_df$T.P.seas.trend[i] = summary(mod_T.P.seas.trend_i)$dispersion
  
  ResDev_df$seas[i] = summary(mod_seas_i)$deviance
  ResDev_df$seas.trend[i] = summary(mod_seas.trend_i)$deviance
  ResDev_df$Temp[i] = summary(mod_T_i)$deviance
  ResDev_df$T.seas[i] = summary(mod_T.seas_i)$deviance
  ResDev_df$T.P.seas[i] = summary(mod_T.P.seas_i)$deviance
  ResDev_df$T.seas.trend[i] = summary(mod_T.seas.trend_i)$deviance
  ResDev_df$T.P.seas.trend[i] = summary(mod_T.P.seas.trend_i)$deviance
}
toc() # roughly 4 min to run 

# Check on the last replicate: fitted + residuals is not the observed count ----
all.equal(as.numeric(mod_T.P.seas_i$fitted.values + mod_T.P.seas_i$residuals),
          as.numeric(mod_T.P.seas_i$model$eggs))      # expect a difference

# Tidy results ----
# all metrics in one long table: rep, metric, dataset, model, value
model_levels <- c("seas", "seas.trend", "Temp", "T.seas", "T.P.seas", "T.seas.trend", "T.P.seas.trend")

cv_long <- dplyr::bind_rows(
  RMSE_fit_df    %>% dplyr::mutate(metric = "RMSE",       dataset = "Fitting"),
  RMSE_val_df    %>% dplyr::mutate(metric = "RMSE",       dataset = "Validation"),
  RMSLE_fit_df   %>% dplyr::mutate(metric = "RMSLE",      dataset = "Fitting"),
  RMSLE_val_df   %>% dplyr::mutate(metric = "RMSLE",      dataset = "Validation"),
  MAE_fit_df     %>% dplyr::mutate(metric = "MAE",        dataset = "Fitting"),
  MAE_val_df     %>% dplyr::mutate(metric = "MAE",        dataset = "Validation"),
  ResDev_df      %>% dplyr::mutate(metric = "ResDev",     dataset = "Fitting"),
  Dispersion_df  %>% dplyr::mutate(metric = "Dispersion", dataset = "Fitting")
) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(model_levels),
                      names_to = "model", values_to = "value") %>%
  dplyr::mutate(model = factor(model, levels = model_levels))

# Summary table ----
# median (IQR) over the repetitions
cv_table <- cv_long %>%
  dplyr::group_by(metric, dataset, model) %>%
  dplyr::summarise(
    median = round(median(value, na.rm = TRUE), 3),
    Q25    = round(quantile(value, 0.25, na.rm = TRUE), 3),
    Q75    = round(quantile(value, 0.75, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  dplyr::mutate(median_IQR = paste0(median, " (", Q25, "-", Q75, ")")) %>%
  dplyr::select(metric, dataset, model, median_IQR) %>%
  tidyr::pivot_wider(names_from = model, values_from = median_IQR) %>%
  dplyr::arrange(metric, dataset) %>%
  as.data.frame()

cv_table

# Main table: validation metrics
table2 <- cv_table %>%
  dplyr::filter(dataset == "Validation") %>%
  dplyr::select(-dataset) %>%
  dplyr::arrange(factor(metric, levels = c("RMSLE", "RMSE", "MAE")))

table2

# Supplementary table: fitting metrics and diagnostics
tableS1 <- cv_table %>%
  dplyr::filter(dataset == "Fitting") %>%
  dplyr::select(-dataset) %>%
  dplyr::arrange(factor(metric, levels = c("RMSLE", "RMSE", "MAE", "ResDev", "Dispersion")))

tableS1

# Paired Wilcoxon test against T.P.seas ----
# same random splits for every model, so the comparison is paired by repetition
# diff = other model - T.P.seas: positive means T.P.seas has the lower error
paired <- cv_long %>%
  dplyr::filter(dataset == "Validation", metric %in% c("RMSLE", "RMSE", "MAE")) %>%
  dplyr::select(metric, rep, model, value) %>%
  tidyr::pivot_wider(names_from = model, values_from = value)

paired_test <- paired %>%
  tidyr::pivot_longer(cols = -c(metric, rep, T.P.seas),
                      names_to = "model", values_to = "other") %>%
  dplyr::mutate(diff = other - T.P.seas) %>%
  dplyr::group_by(metric, model) %>%
  dplyr::summarise(
    median_diff        = round(median(diff), 4),
    reps_TPseas_better = sum(diff > 0),
    p_wilcoxon         = signif(wilcox.test(other, T.P.seas, paired = TRUE)$p.value, 2),
    .groups = "drop"
  ) %>%
  dplyr::arrange(factor(metric, levels = c("RMSLE", "RMSE", "MAE")), p_wilcoxon)

as.data.frame(paired_test)

# Plots ----
# validation: RMSLE, RMSE, MAE
p_val <- cv_long %>%
  dplyr::filter(dataset == "Validation", metric %in% c("RMSLE", "RMSE", "MAE")) %>%
  dplyr::mutate(metric = factor(metric, levels = c("RMSLE", "RMSE", "MAE"))) %>%
  ggplot(aes(x = model, y = value, fill = model)) +
  geom_boxplot(alpha = 0.5, width = 0.4) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = "Validation error") +
  theme_test() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none")

p_val

# fitting: all metrics, log scale
p_fit <- cv_long %>%
  dplyr::filter(dataset == "Fitting") %>%
  ggplot(aes(x = model, y = value, fill = model)) +
  geom_boxplot(alpha = 0.5, width = 0.4) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_y_continuous(trans = "log10") +
  labs(x = NULL, y = "Fitting metric (log scale)") +
  theme_test() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none")

p_fit

ggsave("outputs/FigS_cv_validation.png", p_val, width = 12, height = 4.5, dpi = 300)
ggsave("outputs/FigS_cv_fitting.png",    p_fit, width = 16, height = 4.5, dpi = 300)

# Export tables
write.csv(table2,      "outputs/Table2_cv_validation.csv", row.names = FALSE)
write.csv(tableS1,     "outputs/TableS1_cv_fitting.csv",   row.names = FALSE)
write.csv(paired_test, "outputs/TableS2_cv_wilcoxon.csv",  row.names = FALSE)

# EXPORT DATASETS
saveRDS(bio.matrix_sel,          "data/bio.matrix_cal.rds")
saveRDS(bio.matrix_test_years,   "data/bio.matrix_test_years.rds")
saveRDS(bio.matrix_test_regions, "data/bio.matrix_test_regions.rds")

cat("Exported:\n")
cat("  Calibration : ", nrow(bio.matrix_sel),          "rows\n")
cat("  Test years  : ", nrow(bio.matrix_test_years),   "rows\n")
cat("  Test regions: ", nrow(bio.matrix_test_regions), "rows\n")