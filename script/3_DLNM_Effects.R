# Best model: T.P.seas (temperature + precipitation + seasonality)

# Libraries 
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
library(sf)
library(terra)
library(mapview)
library(tidyterra)
library(patchwork)
library(RColorBrewer)

# 1. Load data ----
bio.matrix_sel <- readRDS("data/bio.matrix_cal.rds")
str(bio.matrix_sel)

unique_locs <- bio.matrix_sel %>%
  ungroup() %>%
  distinct(ID, lon, lat)

loc_vect <- vect(unique_locs,
                 geom = c("lon", "lat"),
                 crs  = "EPSG:4326")
# mapview(st_as_sf(loc_vect))

# 2. Crossbasis definition ----
tlag <- 8

bseas <- ns(bio.matrix_sel$weekn, df = 4)

tknots <- unname(quantile(bio.matrix_sel$medianTweek,
                          probs = c(0.25, 0.75), na.rm = TRUE))

cbTemp <- crossbasis(bio.matrix_sel$medianTweek, lag = tlag,
                     argvar = list(fun = "bs", degree = 2, knots = tknots),
                     arglag = list(knots = c(1, 4)),
                     group  = bio.matrix_sel$Region)

pknots <- unname(quantile(bio.matrix_sel$cumPrecweek,
                          probs = c(0.25, 0.75, 0.9), na.rm = TRUE))

cbPrec <- crossbasis(bio.matrix_sel$cumPrecweek, lag = tlag,
                     argvar = list(fun = "bs", degree = 2, knots = pknots),
                     arglag = list(knots = c(1, 4)),
                     group  = bio.matrix_sel$Region)

cat("tknots:", tknots, "\n")
cat("pknots:", pknots, "\n")

## CHECK 1a: lag structure in the calibration data ----
# crossbasis() builds lags from the previous rows within each group (here Region), in data order.
# A "break" is a row where the previous row of the same region is another trap or not the previous week:
# the first tlag rows after a break take their lagged weather from the wrong trap or the wrong period.
is.unsorted(order(bio.matrix_sel$ID, bio.matrix_sel$year, bio.matrix_sel$week))   # expect FALSE

bio.matrix_sel %>%
  dplyr::group_by(Region) %>%
  dplyr::mutate(brk = dplyr::coalesce(ID != dplyr::lag(ID) |
                                        as.numeric(date - dplyr::lag(date)) > 7, FALSE),
                run = cumsum(brk)) %>%
  dplyr::group_by(Region, run) %>%
  dplyr::mutate(pos = dplyr::row_number()) %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    rows              = dplyr::n(),
    breaks            = sum(brk),
    rows_lag_bridged  = sum(run > 0 & pos <= tlag),
    eggs_lag_bridged  = sum(run > 0 & pos <= tlag & !is.na(eggs)),
    eggs_sum_bridged  = sum(eggs[run > 0 & pos <= tlag], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  as.data.frame()
# rows_lag_bridged: rows whose lags cross a trap or time break; eggs_sum_bridged shows whether they matter

# 3. Fit T.P.seas ----
mod_T.P.seas <- gnm(eggs ~ cbTemp + cbPrec + bseas + year,
                    eliminate = as.factor(ID),
                    data      = bio.matrix_sel,
                    family    = quasipoisson(link = "log"),
                    na.action = "na.exclude")

cat("Converged:", mod_T.P.seas$converged, "\n")
cat("Dispersion:", summary(mod_T.P.seas)$dispersion, "\n")
cat("Residual deviance:", summary(mod_T.P.seas)$deviance, "\n")

# Null model for context, fitted on the same rows as the full model
# (the full model also drops rows with NA lags or NA weather)
used_rows <- setdiff(seq_len(nrow(bio.matrix_sel)), as.integer(mod_T.P.seas$na.action))
cat("Rows with eggs:      ", sum(!is.na(bio.matrix_sel$eggs)), "\n")
cat("Rows used full model:", length(used_rows), "\n")

mod_null <- gnm(eggs ~ 1,
                eliminate = as.factor(ID),
                data      = bio.matrix_sel[used_rows, ],
                family    = quasipoisson(link = "log"),
                na.action = "na.exclude")

stopifnot(length(mod_null$fitted.values) == length(mod_T.P.seas$fitted.values))

cat("Null deviance:     ", summary(mod_null)$deviance, "\n")
cat("Residual deviance: ", summary(mod_T.P.seas)$deviance, "\n")
cat("Deviance explained:", round(
  (1 - summary(mod_T.P.seas)$deviance / summary(mod_null)$deviance) * 100, 1),
  "%\n")

# Save for reporting
cat("Model summary for paper:\n")
cat("  Quasi-Poisson dispersion:     ", round(summary(mod_T.P.seas)$dispersion, 1), "\n")
cat("  Null deviance (trap-adjusted):", round(summary(mod_null)$deviance, 0), "\n")
cat("  Residual deviance:            ", round(summary(mod_T.P.seas)$deviance, 0), "\n")
cat("  Deviance explained:           ", round(
  (1 - summary(mod_T.P.seas)$deviance / summary(mod_null)$deviance) * 100, 1), "%\n")

# 4. Temperature effects ----
myCen_T <- 15
cptmean <- crosspred(cbTemp, mod_T.P.seas, cen = myCen_T, by = 1, at = 6:30)

# Overall cumulative effect
png("outputs/temp_overall.png", width = 800, height = 600, res = 120)
plot(cptmean, ptype = "overall",
     xlab = "Temperature (°C)", ylab = "IRR",
     col  = 2, main = "Temperature — overall cumulative effect")
dev.off()

# Lag profiles at fixed temperatures
png("outputs/temp_lag_profiles.png", width = 900, height = 800, res = 120)
layout(matrix(1:4, ncol = 2, nrow = 2))
for (tmp in c(6, 14, 22, 30)) {
  plot(cptmean, var = tmp,
       xlab = "Lag (weeks)", ylab = "IRR", col = 2,
       main = paste("T =", tmp, "°C"), ylim = c(0, 3))
}
layout(1)
dev.off()

# Exposure-response at fixed lags
png("outputs/temp_exposure_response.png", width = 900, height = 800, res = 120)
layout(matrix(1:4, ncol = 2, nrow = 2))
for (lg in c(0, 2, 4, 8)) {
  plot(cptmean, lag = lg,
       xlab = "Temperature (°C)", ylab = "IRR", col = 2,
       main = paste("Lag =", lg, "weeks"), ylim = c(0, 3))
}
layout(1)
dev.off()

plot(cptmean, ptype = "3d",
     xlab = "Temperature (°C)", ylab = "Lag (weeks)", zlab = "IRR",
     main = "Temperature lag-response surface")


# 5. Precipitation effects ----
myCen_P <- 0
P_max   <- unname(quantile(bio.matrix_sel$cumPrecweek, 0.95, na.rm = TRUE))
cat("Precipitation Q95:", P_max, "mm\n")

cppmean <- crosspred(cbPrec, mod_T.P.seas,
                     cen = myCen_P, by = 1, at = 0:ceiling(P_max))

# Overall cumulative effect
png("outputs/precip_overall.png", width = 800, height = 600, res = 120)
plot(cppmean, ptype = "overall",
     xlab = "Precipitation (mm)", ylab = "IRR",
     col  = 2, main = "Precipitation — overall cumulative effect",
     ylim = c(0.5, 3))
dev.off()

# Lag profiles at fixed precipitation values
png("outputs/precip_lag_profiles.png", width = 900, height = 800, res = 120)
layout(matrix(1:4, ncol = 2, nrow = 2))
for (prc in c(5, 10, 20, 50)) {
  plot(cppmean, var = prc,
       xlab = "Lag (weeks)", ylab = "IRR", col = 2,
       main = paste("Precip =", prc, "mm"), ylim = c(0.75, 1.25))
}
layout(1)
dev.off()

# Exposure-response at fixed lags
png("outputs/precip_exposure_response.png", width = 900, height = 800, res = 120)
layout(matrix(1:4, ncol = 2, nrow = 2))
for (lg in c(0, 2, 4, 8)) {
  plot(cppmean, lag = lg,
       xlab = "Precipitation (mm)", ylab = "IRR", col = 2,
       main = paste("Lag =", lg, "weeks"), ylim = c(0, 2))
}
layout(1)
dev.off()

plot(cppmean, ptype = "3d",
     xlab = "Precipitation (mm)", ylab = "Lag (weeks)", zlab = "IRR",
     main = "Precipitation lag-response surface")

# 6. Validation on test sets ----

# Load test datasets
bio.matrix_test_years   <- readRDS("data/bio.matrix_test_years.rds") %>%
  dplyr::mutate(weekn = as.numeric(week))

unique_locs <- bio.matrix_test_years %>%
  ungroup() %>%
  distinct(ID, lon, lat)

loc_vect <- vect(unique_locs,
                 geom = c("lon", "lat"),
                 crs  = "EPSG:4326")
# mapview(st_as_sf(loc_vect))

bio.matrix_test_regions <- readRDS("data/bio.matrix_test_regions.rds") %>%
  dplyr::mutate(weekn = as.numeric(week))

unique_locs <- bio.matrix_test_regions %>%
  ungroup() %>%
  distinct(ID, lon, lat)

loc_vect <- vect(unique_locs,
                 geom = c("lon", "lat"),
                 crs  = "EPSG:4326")
# mapview(st_as_sf(loc_vect))

cat("Test years   — rows:", nrow(bio.matrix_test_years),
    "| traps:", dplyr::n_distinct(bio.matrix_test_years$ID),
    "| regions:", dplyr::n_distinct(bio.matrix_test_years$Region), "\n")

cat("Test regions — rows:", nrow(bio.matrix_test_regions),
    "| traps:", dplyr::n_distinct(bio.matrix_test_regions$ID),
    "| regions:", dplyr::n_distinct(bio.matrix_test_regions$Region), "\n")

## CHECK 1b: row order and lag breaks in the test sets ----
# the test sets are never sorted: if FALSE is not returned, the lags are scrambled
is.unsorted(order(bio.matrix_test_years$ID, bio.matrix_test_years$year, bio.matrix_test_years$week))     # expect FALSE
is.unsorted(order(bio.matrix_test_regions$ID, bio.matrix_test_regions$year, bio.matrix_test_regions$week)) # expect FALSE

# temporal test set: breaks include the 2017 -> 2024 jump and missing years within traps
bio.matrix_test_years %>%
  dplyr::group_by(Region) %>%
  dplyr::mutate(brk = dplyr::coalesce(ID != dplyr::lag(ID) |
                                        as.numeric(date - dplyr::lag(date)) > 7, FALSE),
                run = cumsum(brk)) %>%
  dplyr::group_by(Region, run) %>%
  dplyr::mutate(pos = dplyr::row_number()) %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    rows             = dplyr::n(),
    breaks           = sum(brk),
    rows_lag_bridged = sum(run > 0 & pos <= tlag),
    eggs_sum_bridged = sum(eggs[run > 0 & pos <= tlag], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  as.data.frame()

# spatial test set
bio.matrix_test_regions %>%
  dplyr::group_by(Region) %>%
  dplyr::mutate(brk = dplyr::coalesce(ID != dplyr::lag(ID) |
                                        as.numeric(date - dplyr::lag(date)) > 7, FALSE),
                run = cumsum(brk)) %>%
  dplyr::group_by(Region, run) %>%
  dplyr::mutate(pos = dplyr::row_number()) %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    rows             = dplyr::n(),
    breaks           = sum(brk),
    rows_lag_bridged = sum(run > 0 & pos <= tlag),
    eggs_sum_bridged = sum(eggs[run > 0 & pos <= tlag], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  as.data.frame()

# Fixed-effect linear predictor on calibration data  
coef_fixed <- coef(mod_T.P.seas)[!is.na(coef(mod_T.P.seas))]

X_cal <- cbind(
  as.matrix(as.data.frame(cbTemp)),
  as.matrix(as.data.frame(cbPrec)),
  as.matrix(as.data.frame(bseas)),
  bio.matrix_sel$year
)

lp_fixed_cal <- as.numeric(X_cal %*% coef_fixed)

# Expand fitted values to full data length (na.exclude compresses them)
na_idx      <- as.integer(mod_T.P.seas$na.action)
used_idx    <- setdiff(seq_len(nrow(bio.matrix_sel)), na_idx)
lp_full_cal <- rep(NA_real_, nrow(bio.matrix_sel))
lp_full_cal[used_idx] <- log(mod_T.P.seas$fitted.values)

cat("Rows used in model:", length(used_idx), "\n")
cat("Non-NA in lp_full: ", sum(!is.na(lp_full_cal)), "\n")

## 6.1 Trap-specific intercepts ----
# taken directly from gnm (eliminated coefficients); reconstructing them from
# log(fitted) is biased where fitted values hit the floor of the log link
elim <- attr(coef(mod_T.P.seas), "eliminated")
trap_intercepts_df <- data.frame(ID = as.numeric(names(elim)), intercept = as.numeric(elim))
stopifnot(nrow(trap_intercepts_df) == dplyr::n_distinct(bio.matrix_sel$ID))

# check: fixed part + trap intercept reproduces the fitted values wherever they are not floored
lp_rec <- lp_fixed_cal + trap_intercepts_df$intercept[match(bio.matrix_sel$ID, trap_intercepts_df$ID)]
ok     <- !is.na(lp_full_cal) & lp_full_cal > log(.Machine$double.eps) + 1e-6
max(abs(lp_rec[ok] - lp_full_cal[ok]))                       # expect ~0
sum(!is.na(lp_full_cal) & !ok)                               # rows at the floor

# trap baselines on the count scale (not expected to be 0)
summary(exp(trap_intercepts_df$intercept))

cat("Trap intercepts computed:", nrow(trap_intercepts_df), "\n")

# Population-average intercept for unknown trap IDs
# row-weighted mean of the trap intercepts
n_used_trap   <- table(bio.matrix_sel$ID[!is.na(lp_full_cal)])
avg_intercept <- weighted.mean(trap_intercepts_df$intercept,
                               as.numeric(n_used_trap[as.character(trap_intercepts_df$ID)]))
cat("Population-average intercept:", round(avg_intercept, 4), "\n")

##6.2 Filter test_years to known trap IDs only ----
known   <- intersect(unique(bio.matrix_test_years$ID), unique(bio.matrix_sel$ID))
unknown <- setdiff(unique(bio.matrix_test_years$ID),   unique(bio.matrix_sel$ID))
cat("Test years — known intercept:", length(known),
    "| excluded (not in calibration):", length(unknown),
    "| unknown IDs:", paste(unknown, collapse = ", "), "\n")

bio.matrix_test_years_known <- bio.matrix_test_years %>%
  dplyr::filter(ID %in% known)

# Prediction function  
predict_with_intercept <- function(newdata, mod, bseas_cal,
                                   tknots, pknots,
                                   intercept_lookup,
                                   avg_intercept,
                                   tlag = 8) {
  cbT_new <- crossbasis(newdata$medianTweek, lag = tlag,
                        argvar = list(fun = "bs", degree = 2, knots = tknots),
                        arglag = list(knots = c(1, 4)),
                        group  = newdata$Region)
  
  cbP_new <- crossbasis(newdata$cumPrecweek, lag = tlag,
                        argvar = list(fun = "bs", degree = 2, knots = pknots),
                        arglag = list(knots = c(1, 4)),
                        group  = newdata$Region)
  
  bseas_new <- predict(bseas_cal, newx = newdata$weekn)
  
  X <- cbind(
    as.matrix(as.data.frame(cbT_new)),
    as.matrix(as.data.frame(cbP_new)),
    as.matrix(bseas_new),
    newdata$year
  )
  
  coef_fixed <- coef(mod)[!is.na(coef(mod))]
  
  if (ncol(X) != length(coef_fixed))
    stop(paste("Dimension mismatch: X has", ncol(X),
               "cols but model has", length(coef_fixed), "fixed coefs"))
  
  lp_fixed <- as.numeric(X %*% coef_fixed)
  
  alpha <- dplyr::left_join(
    data.frame(ID = newdata$ID), intercept_lookup, by = "ID"
  ) %>% dplyr::pull(intercept)
  
  alpha[is.na(alpha)] <- avg_intercept
  
  as.numeric(exp(lp_fixed + alpha))
}

## 6.3 Evaluation metrics ----
eval_metrics <- function(obs, pred, label) {
  obs  <- as.numeric(obs)
  pred <- as.numeric(pred)
  idx  <- which(!is.na(obs) & !is.na(pred) & is.finite(pred))
  obs  <- obs[idx]
  pred <- pred[idx]
  cat(label, "\n")
  cat("  N:        ", length(idx), "\n")
  cat("  RMSE:     ", round(sqrt(mean((obs - pred)^2)),                   3), "\n")
  cat("  MAE:      ", round(mean(abs(obs - pred)),                        3), "\n")
  cat("  RMSLE:    ", round(sqrt(mean((log(pred + 1) - log(obs + 1))^2)), 4), "\n")
  cat("  Pearson r:", round(cor(obs, pred, method = "pearson"),           4), "\n")
  cat("  Spearman rho:", round(cor(obs, pred, method = "spearman"),      4), "\n\n")
}

## 6.4 Run predictions ----
pred_years_known <- predict_with_intercept(
  bio.matrix_test_years_known, mod_T.P.seas, bseas,
  tknots, pknots, trap_intercepts_df, avg_intercept
)

pred_regions <- predict_with_intercept(
  bio.matrix_test_regions, mod_T.P.seas, bseas,
  tknots, pknots, trap_intercepts_df, avg_intercept
)

## CHECK 2: extreme predictions and extrapolation beyond the calibration range ----
summary(pred_years_known)
summary(pred_regions)
cat("Max fitted value in calibration:", round(max(mod_T.P.seas$fitted.values), 1), "\n")

# calibration ranges of the exposures (bs() extrapolates quadratically beyond these)
rng_T <- range(bio.matrix_sel$medianTweek, na.rm = TRUE)
rng_P <- range(bio.matrix_sel$cumPrecweek, na.rm = TRUE)
cat("Calibration temperature range (°C):", round(rng_T, 2), "\n")
cat("Calibration precipitation range (mm):", round(rng_P, 2), "\n")

# share of test weeks outside the calibration range
cat("Test years, temperature outside range:  ",
    round(mean(bio.matrix_test_years_known$medianTweek < rng_T[1] |
                 bio.matrix_test_years_known$medianTweek > rng_T[2], na.rm = TRUE), 4), "\n")
cat("Test years, precipitation outside range:",
    round(mean(bio.matrix_test_years_known$cumPrecweek < rng_P[1] |
                 bio.matrix_test_years_known$cumPrecweek > rng_P[2], na.rm = TRUE), 4), "\n")
cat("Test regions, temperature outside range:  ",
    round(mean(bio.matrix_test_regions$medianTweek < rng_T[1] |
                 bio.matrix_test_regions$medianTweek > rng_T[2], na.rm = TRUE), 4), "\n")
cat("Test regions, precipitation outside range:",
    round(mean(bio.matrix_test_regions$cumPrecweek < rng_P[1] |
                 bio.matrix_test_regions$cumPrecweek > rng_P[2], na.rm = TRUE), 4), "\n")

# the ten largest predictions in the temporal test set and their exposures
bio.matrix_test_years_known %>%
  dplyr::mutate(pred = pred_years_known) %>%
  dplyr::arrange(dplyr::desc(pred)) %>%
  dplyr::select(ID, Region, year, week, eggs, pred, medianTweek, cumPrecweek) %>%
  head(10) %>%
  as.data.frame()

# the ten largest predictions in the spatial test set
bio.matrix_test_regions %>%
  dplyr::mutate(pred = pred_regions) %>%
  dplyr::arrange(dplyr::desc(pred)) %>%
  dplyr::select(ID, Region, year, week, eggs, pred, medianTweek, cumPrecweek) %>%
  head(10) %>%
  as.data.frame()

# linear year term: multiplicative shift when predicting outside 2018-2023
b_year <- coef(mod_T.P.seas)["year"]
cat("Year coefficient:", round(b_year, 4), "\n")
cat("Multiplier 2010 vs 2018:", round(exp(b_year * (2010 - 2018)), 3),
    "| 2024 vs 2023:", round(exp(b_year * (2024 - 2023)), 3), "\n")

## 6.5 Pooled evaluation ----
eval_metrics(
  bio.matrix_test_years_known$eggs[bio.matrix_test_years_known$year %in% 2010:2017],
  pred_years_known[bio.matrix_test_years_known$year %in% 2010:2017],
  "Test set 1a — retrospective 2010-2017"
)

eval_metrics(
  bio.matrix_test_years_known$eggs[bio.matrix_test_years_known$year == 2024],
  pred_years_known[bio.matrix_test_years_known$year == 2024],
  "Test set 1b — forward validation 2024"
)

eval_metrics(bio.matrix_test_regions$eggs, pred_regions,
             "Test set 2 — held-out regions")

## 6.6 Per-trap evaluation — test_years ----
per_trap_r <- bio.matrix_test_years_known %>%
  dplyr::mutate(pred = pred_years_known) %>%
  dplyr::filter(!is.na(eggs), is.finite(pred)) %>%
  dplyr::group_by(ID, Region) %>%
  dplyr::summarise(
    n         = dplyr::n(),
    Pearson_r = cor(eggs, pred, method = "pearson", use = "complete.obs"),
    RMSLE     = sqrt(mean((log(pred + 1) - log(eggs + 1))^2)),
    .groups   = "drop"
  ) %>%
  dplyr::filter(n >= 10) %>%
  dplyr::arrange(desc(Pearson_r))

cat("Per-trap validation summary:\n")
cat("  Median Pearson r:", round(median(per_trap_r$Pearson_r, na.rm = TRUE), 3), "\n")
cat("  IQR:             ", round(quantile(per_trap_r$Pearson_r, 0.25, na.rm = TRUE), 3),
    "—", round(quantile(per_trap_r$Pearson_r, 0.75, na.rm = TRUE), 3), "\n")
cat("  Median RMSLE:    ", round(median(per_trap_r$RMSLE, na.rm = TRUE), 3), "\n")
cat("  Traps r > 0.5:   ", sum(per_trap_r$Pearson_r > 0.5, na.rm = TRUE),
    "/", sum(!is.na(per_trap_r$Pearson_r)), "\n\n")

## 6.6b Per-trap Spearman — test_years ----
per_trap_rho <- bio.matrix_test_years_known %>%
  dplyr::mutate(pred = pred_years_known) %>%
  dplyr::filter(!is.na(eggs), is.finite(pred)) %>%
  dplyr::group_by(ID, Region) %>%
  dplyr::filter(dplyr::n() >= 10, sd(eggs) > 0, sd(pred) > 0) %>%
  dplyr::summarise(
    n     = dplyr::n(),
    rho   = cor(eggs, pred, method = "spearman"),
    p_rho = cor.test(eggs, pred, method = "spearman", exact = FALSE)$p.value,
    RMSLE = sqrt(mean((log(pred + 1) - log(eggs + 1))^2)),
    .groups = "drop"
  )

cat("Traps evaluated:", nrow(per_trap_rho), "\n")
cat("Median rho:", round(median(per_trap_rho$rho), 2),
    "| IQR:", round(quantile(per_trap_rho$rho, 0.25), 2), "-", round(quantile(per_trap_rho$rho, 0.75), 2), "\n")
cat("Traps with p < 0.001:", sum(per_trap_rho$p_rho < 0.001), "/", nrow(per_trap_rho), "\n")

# per region
per_trap_rho %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(n_traps      = dplyr::n(),
                   median_rho   = round(median(rho), 2),
                   median_RMSLE = round(median(RMSLE), 2),
                   .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(median_rho)) %>%
  as.data.frame()

# pooled across all known traps and test years
with(bio.matrix_test_years_known %>% dplyr::mutate(pred = pred_years_known) %>%
       dplyr::filter(!is.na(eggs), is.finite(pred)),
     round(cor(eggs, pred, method = "spearman"), 2))

## 6.7 Per-region summary — test_years ----
per_trap_r %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    n_traps      = dplyr::n(),
    median_r     = round(median(Pearson_r, na.rm = TRUE), 3),
    median_RMSLE = round(median(RMSLE,     na.rm = TRUE), 3),
    .groups      = "drop"
  ) %>%
  dplyr::arrange(desc(median_r)) %>%
  as.data.frame()

##6.8 Per-region breakdown — test_regions ----
bio.matrix_test_regions %>%
  dplyr::mutate(pred = pred_regions) %>%
  dplyr::filter(!is.na(eggs), is.finite(pred)) %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    n         = dplyr::n(),
    RMSLE     = round(sqrt(mean((log(pred + 1) - log(eggs + 1))^2)), 4),
    Pearson_r = round(cor(eggs, pred, method = "pearson"), 4),
    .groups   = "drop"
  ) %>%
  dplyr::arrange(RMSLE) %>%
  as.data.frame()

# Table 3: Per-region — temporal test set  
table3 <- bio.matrix_test_years_known %>%
  dplyr::mutate(pred = pred_years_known) %>%
  dplyr::filter(!is.na(eggs), is.finite(pred)) %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    n_traps      = dplyr::n_distinct(ID),
    n_trap_weeks = dplyr::n(),
    Pearson_r    = round(cor(eggs, pred, method = "pearson"), 3),
    RMSLE        = round(sqrt(mean((log(pred + 1) - log(eggs + 1))^2)), 3),
    .groups      = "drop"
  ) %>%
  dplyr::arrange(RMSLE)

# add per-trap medians (the values quoted in the text): pooled RMSLE and median per-trap RMSLE differ
table3 <- table3 %>%
  dplyr::left_join(per_trap_rho %>%
                     dplyr::group_by(Region) %>%
                     dplyr::summarise(median_rho   = round(median(rho), 2),
                                      median_RMSLE = round(median(RMSLE), 2),
                                      .groups = "drop"),
                   by = "Region")

write.csv(table3, "outputs/Table3_temporal_validation.csv", row.names = FALSE)
print(table3)

# Table 4: Per-region — spatial test set  
table4 <- bio.matrix_test_regions %>%
  dplyr::mutate(pred = pred_regions) %>%
  dplyr::filter(!is.na(eggs), is.finite(pred)) %>%
  dplyr::group_by(Region, Country) %>%
  dplyr::summarise(
    n_traps      = dplyr::n_distinct(ID),
    n_trap_weeks = dplyr::n(),
    Pearson_r    = round(cor(eggs, pred, method = "pearson"), 3),
    RMSLE        = round(sqrt(mean((log(pred + 1) - log(eggs + 1))^2)), 3),
    .groups      = "drop"
  ) %>%
  dplyr::arrange(RMSLE)

write.csv(table4, "outputs/Table4_spatial_validation.csv", row.names = FALSE)
print(table4)

## 6.9 Visual check — example trap ----
# representative trap: median Pearson r among evaluated traps
ID_check <- per_trap_r %>%
  dplyr::filter(!is.na(Pearson_r)) %>%
  dplyr::arrange(Pearson_r) %>%
  dplyr::slice(ceiling(dplyr::n() / 2)) %>%
  dplyr::pull(ID)

stopifnot(ID_check %in% bio.matrix_test_years_known$ID)
cat("ID_check:", ID_check, "\n")

check_df <- bio.matrix_test_years_known %>%
  dplyr::filter(ID == ID_check) %>%
  dplyr::mutate(pred = pred_years_known[bio.matrix_test_years_known$ID == ID_check])

plot(check_df$date, check_df$eggs,
     pch = 16, xlab = "Date", ylab = "Eggs",
     main = paste("ID", ID_check, "— test years"))
lines(check_df$date, check_df$pred, col = "darkgreen", lwd = 2)
legend("topright", legend = c("Observed", "Predicted"),
       col = c("black", "darkgreen"), pch = c(16, NA), lty = c(NA, 1))

cat("Pearson r ID", ID_check, ":",
    round(cor(check_df$eggs, check_df$pred,
              use = "complete.obs", method = "pearson"), 4), "\n")


## 6.10 Visual check — spatial test set, best trap per region ----

# Add predictions to test_regions
bio.matrix_test_regions_pred <- bio.matrix_test_regions %>%
  dplyr::mutate(pred = pred_regions)

# Find best trap per region by Pearson r
best_traps <- bio.matrix_test_regions_pred %>%
  dplyr::filter(!is.na(eggs), is.finite(pred)) %>%
  dplyr::group_by(Region, ID) %>%
  dplyr::summarise(
    n         = dplyr::n(),
    Pearson_r = cor(eggs, pred, method = "pearson", use = "complete.obs"),
    .groups   = "drop"
  ) %>%
  dplyr::filter(n >= 10) %>%
  dplyr::group_by(Region) %>%
  dplyr::slice_max(Pearson_r, n = 1) %>%
  dplyr::ungroup()

cat("Best trap per region:\n")
print(as.data.frame(best_traps))

## 6.10b Spearman — spatial test set ----
per_trap_rho_reg <- bio.matrix_test_regions_pred %>%
  dplyr::filter(!is.na(eggs), is.finite(pred)) %>%
  dplyr::group_by(ID, Region) %>%
  dplyr::filter(dplyr::n() >= 10, sd(eggs) > 0, sd(pred) > 0) %>%
  dplyr::summarise(
    n     = dplyr::n(),
    rho   = cor(eggs, pred, method = "spearman"),
    p_rho = cor.test(eggs, pred, method = "spearman", exact = FALSE)$p.value,
    RMSLE = sqrt(mean((log(pred + 1) - log(eggs + 1))^2)),
    .groups = "drop"
  )

cat("Traps evaluated:", nrow(per_trap_rho_reg), "of", dplyr::n_distinct(bio.matrix_test_regions$ID), "\n")
cat("Median rho:", round(median(per_trap_rho_reg$rho), 2),
    "| IQR:", round(quantile(per_trap_rho_reg$rho, 0.25), 2), "-", round(quantile(per_trap_rho_reg$rho, 0.75), 2), "\n")
cat("Traps with p < 0.001:", sum(per_trap_rho_reg$p_rho < 0.001), "/", nrow(per_trap_rho_reg), "\n")

# per region: pooled within region, as in Table 4
rho_region <- bio.matrix_test_regions_pred %>%
  dplyr::filter(!is.na(eggs), is.finite(pred)) %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(Spearman_rho = round(cor(eggs, pred, method = "spearman"), 2),
                   .groups = "drop")

table4 <- table4 %>% dplyr::left_join(rho_region, by = "Region")
write.csv(table4, "outputs/Table4_spatial_validation.csv", row.names = FALSE)
as.data.frame(table4)

# representative traps from spatial test set — low to high Pearson r
# region labels and r are taken from the current predictions
selected_ids <- c(22983, 14499, 17040, 16578, 9987, 4208)

selected_traps <- bio.matrix_test_regions_pred %>%
  dplyr::filter(ID %in% selected_ids, !is.na(eggs), is.finite(pred)) %>%
  dplyr::group_by(Region, ID) %>%
  dplyr::summarise(Pearson_r = round(cor(eggs, pred, method = "pearson"), 2),
                   .groups = "drop") %>%
  dplyr::mutate(Region = dplyr::recode(Region,
                                       "Sicilia"                    = "Sicily",
                                       "Provence-Alpes-Côte d’Azur" = "Côte d’Azur",
                                       "Région lémanique"           = "Région Lémanique")) %>%
  dplyr::arrange(Pearson_r) %>%
  as.data.frame()

stopifnot(nrow(selected_traps) == length(selected_ids))

cat("Selected traps:\n")
print(selected_traps)

# Plot
png("outputs/Fig_spatial_validation_6traps.png",
    width = 4800, height = 3200, res = 200)

par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

for (i in seq_len(nrow(selected_traps))) {
  
  id_check  <- selected_traps$ID[i]
  reg       <- selected_traps$Region[i]
  r_val     <- selected_traps$Pearson_r[i]
  
  df <- bio.matrix_test_regions_pred %>%
    dplyr::filter(ID == id_check, !is.na(eggs))
  
  # y axis upper limit — leave some headroom above max
  y_max <- max(c(df$eggs, df$pred), na.rm = TRUE) * 1.15
  
  plot(df$date, df$eggs,
       pch  = 16, cex = 0.8,
       xlab = "Date", ylab = "Eggs per trap per week",
       ylim = c(0, y_max),
       main = paste0(reg, "  (r = ", r_val, ")"),
       cex.main = 0.9, cex.axis = 0.85, cex.lab = 0.85)
  
  lines(df$date, df$pred, col = "darkgreen", lwd = 2)
  
  if (i == 1) {
    legend("topright",
           legend = c("Observed", "Predicted"),
           col    = c("black", "darkgreen"),
           pch    = c(16, NA), lty = c(NA, 1),
           bty    = "n", cex = 0.8)
  }
}

par(mfrow = c(1, 1))
dev.off()

# 7. Spatial extrapolation — oviposition seasonality ----
# Load spatial environmental matrix
env.matrix <- readRDS("explorative//envMatrix_2018_2023.RDS") %>%
  dplyr::rename(medianTweek = tas, cumPrecweek = tp) %>%
  dplyr::filter(year > 2018) %>%
  dplyr::mutate(weekn = as.numeric(week)) %>%
  as_tibble()

cat("Env matrix rows:", nrow(env.matrix), "\n")
cat("Years:          ", paste(unique(env.matrix$year), collapse = ", "), "\n")
cat("Columns:        ", paste(names(env.matrix), collapse = ", "), "\n")

## CHECK 1c: row order of the grid ----
# with a single group, lag 1 must be the same cell one week earlier:
# the matrix must be sorted cell by cell (x, y), then by year and week
head(env.matrix %>% dplyr::select(x, y, year, week), 10)
is.unsorted(order(env.matrix$x, env.matrix$y, env.matrix$year, env.matrix$week))   # expect FALSE
cat("Grid cells:", dplyr::n_distinct(paste(env.matrix$x, env.matrix$y)), "\n")
# note: filter(year > 2018) drops 2018, so early-2019 lags cross into the previous cell's series

## CHECK 2b: grid exposures outside the calibration range ----
cat("Grid, temperature outside range:  ",
    round(mean(env.matrix$medianTweek < rng_T[1] | env.matrix$medianTweek > rng_T[2], na.rm = TRUE), 4), "\n")
cat("Grid, precipitation outside range:",
    round(mean(env.matrix$cumPrecweek < rng_P[1] | env.matrix$cumPrecweek > rng_P[2], na.rm = TRUE), 4), "\n")

## 7.1 Predict on spatial grid ----
# env.matrix has no trap IDs — use avg_intercept for all locations
# Region column needed for crossbasis group argument — check if present
cat("Region column present:", "Region" %in% names(env.matrix), "\n")

# If no Region column, create a dummy — crossbasis group just breaks lag continuity
# at region boundaries; a single group is valid for a spatial grid
if (!"Region" %in% names(env.matrix)) {
  env.matrix <- env.matrix %>% dplyr::mutate(Region = "grid")
  cat("Region column added as dummy\n")
}

# Build crossbases — same knots as calibration
cbT_grid <- crossbasis(env.matrix$medianTweek, lag = tlag,
                       argvar = list(fun = "bs", degree = 2, knots = tknots),
                       arglag = list(knots = c(1, 4)),
                       group  = env.matrix$Region)

cbP_grid <- crossbasis(env.matrix$cumPrecweek, lag = tlag,
                       argvar = list(fun = "bs", degree = 2, knots = pknots),
                       arglag = list(knots = c(1, 4)),
                       group  = env.matrix$Region)

# Project seasonality spline onto calibration basis
bseas_grid <- predict(bseas, newx = env.matrix$weekn)

# Build model matrix
X_grid <- cbind(
  as.matrix(as.data.frame(cbT_grid)),
  as.matrix(as.data.frame(cbP_grid)),
  as.matrix(bseas_grid),
  env.matrix$year
)

if (ncol(X_grid) != length(coef_fixed))
  stop(paste("Dimension mismatch: X_grid has", ncol(X_grid),
             "cols but model has", length(coef_fixed), "fixed coefs"))

# Predict using fixed effects + population-average intercept
lp_grid          <- as.numeric(X_grid %*% coef_fixed)
env.matrix$eggsPred <- as.numeric(exp(lp_grid + avg_intercept))
env.matrix$eggsPred <- round(pmax(env.matrix$eggsPred, 0))

cat("Prediction range:", round(range(env.matrix$eggsPred, na.rm = TRUE), 1), "\n")
print(quantile(env.matrix$eggsPred, probs = c(0.5, 0.9, 0.99, 0.999, 1), na.rm = TRUE))

## 7.2 Season onset and offset per grid cell per year ----
# Threshold: first/last week with predicted eggs > 10
eggThresh <- 10
seas.df <- env.matrix %>%
  dplyr::group_by(x, y, year) %>%
  dplyr::filter(eggsPred > eggThresh) %>%
  dplyr::summarise(
    start = min(week),
    end   = max(week),
    .groups = "drop"
  )

cat("Grid cells with detectable season:", dplyr::n_distinct(seas.df %>% dplyr::select(x, y)), "\n")

seas.df %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    n_cells       = dplyr::n(),
    mean_start    = round(mean(start), 1),
    mean_end      = round(mean(end),   1),
    mean_length   = round(mean(end - start), 1),
    .groups       = "drop"
  ) %>%
  as.data.frame()

## 7.3 Average season metrics across years ----
seas.av.df <- seas.df %>%
  dplyr::group_by(x, y) %>%
  dplyr::summarise(
    m.start  = mean(start),
    m.end    = mean(end),
    s.length = mean(end - start),
    .groups  = "drop"
  )

## 7.4 Rasterise and plot ----
slength <- seas.av.df %>%
  dplyr::select(x, y, s.length) %>%
  terra::rast(type = "xyz")

onset <- seas.av.df %>%
  dplyr::select(x, y, m.start) %>%
  terra::rast(type = "xyz")

offset <- seas.av.df %>%
  dplyr::select(x, y, m.end) %>%
  terra::rast(type = "xyz")

par(mfrow = c(1, 3))
plot(slength / 4,  main = "Average season length (months)")
plot(onset,        main = "Average season onset (week)")
plot(offset,       main = "Average season offset (week)")
par(mfrow = c(1, 1))

## 7.5 Maps plot ----
## 7.5 Maps plot ----
cat("Season length range (weeks):", range(seas.av.df$s.length, na.rm = TRUE), "\n")
cat("Onset range (week):         ", range(seas.av.df$m.start,  na.rm = TRUE), "\n")
cat("Offset range (week):        ", range(seas.av.df$m.end,    na.rm = TRUE), "\n")

# the grid is in longitude/latitude (EPSG:4326)
crs(slength) <- "EPSG:4326"
crs(onset)   <- "EPSG:4326"
crs(offset)  <- "EPSG:4326"

slength_m <- slength / 4                       # season length in months

# boundaries clipped to the extent of the cells with data (avoids artefacts at the map edges)
e <- as.vector(terra::ext(terra::trim(slength)))
box_sfc <- st_as_sfc(st_bbox(c(xmin = e[["xmin"]], xmax = e[["xmax"]],
                               ymin = e[["ymin"]], ymax = e[["ymax"]]),
                             crs = st_crs(4326)))

countries_all <- eurostat::get_eurostat_geospatial(output_class = "sf", resolution = "01",
                                                   nuts_level = 0, year = 2024)
nuts2_all     <- eurostat::get_eurostat_geospatial(output_class = "sf", resolution = "01",
                                                   nuts_level = 2, year = 2024)

countries_lines <- st_intersection(st_boundary(countries_all), box_sfc)
nuts2_lines     <- st_intersection(st_boundary(nuts2_all),     box_sfc)

# Common theme
map_theme <- theme_classic() +
  theme(
    legend.position       = "bottom",
    legend.title.position = "top",
    text                  = element_text(size = 12),
    axis.text             = element_text(size = 9),
    legend.text           = element_text(size = 10),
    legend.key.width      = unit(1.5, "cm"),
    legend.key.height     = unit(0.35, "cm")
  )

# one palette per panel; darker colours = more favourable in all panels:
# longer season, earlier onset, later offset

# Panel A: season length (months)
pA <- ggplot() +
  geom_spatraster(data = slength_m) +
  geom_sf(data = countries_lines, colour = "black", linewidth = 0.25, alpha = 0.7) +
  geom_sf(data = nuts2_lines,     colour = "black", linewidth = 0.18, alpha = 0.3) +
  scale_fill_viridis_b(option = "rocket", direction = -1,
                       limits = c(0, 8), breaks = 0:8,
                       oob = scales::squish, na.value = "transparent",
                       name = "Season length (months)") +
  coord_sf(xlim = c(e[["xmin"]], e[["xmax"]]), ylim = c(e[["ymin"]], e[["ymax"]]), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  map_theme

# Panel B: season onset (ISO week)
pB <- ggplot() +
  geom_spatraster(data = onset) +
  geom_sf(data = countries_lines, colour = "black", linewidth = 0.25, alpha = 0.7) +
  geom_sf(data = nuts2_lines,     colour = "black", linewidth = 0.18, alpha = 0.3) +
  scale_fill_viridis_b(option = "mako", direction = 1,
                       limits = c(16, 32), breaks = seq(16, 32, 2),
                       oob = scales::squish, na.value = "transparent",
                       name = "Season onset (ISO week)") +
  coord_sf(xlim = c(e[["xmin"]], e[["xmax"]]), ylim = c(e[["ymin"]], e[["ymax"]]), expand = FALSE) +
  labs(x = "Longitude", y = NULL) +
  map_theme

# Panel C: season offset (ISO week)
pC <- ggplot() +
  geom_spatraster(data = offset) +
  geom_sf(data = countries_lines, colour = "black", linewidth = 0.25, alpha = 0.7) +
  geom_sf(data = nuts2_lines,     colour = "black", linewidth = 0.18, alpha = 0.3) +
  scale_fill_viridis_b(option = "viridis", direction = -1,
                       limits = c(30, 46), breaks = seq(30, 46, 2),
                       oob = scales::squish, na.value = "transparent",
                       name = "Season offset (ISO week)") +
  ggspatial::annotation_scale(location = "br", width_hint = 0.25) +
  coord_sf(xlim = c(e[["xmin"]], e[["xmax"]]), ylim = c(e[["ymin"]], e[["ymax"]]), expand = FALSE) +
  labs(x = "Longitude", y = NULL) +
  map_theme

# Combine and save (caption goes in the manuscript, not in the image)
fig_seas <- (pA | pB | pC) + plot_annotation(tag_levels = "A")
fig_seas
ggsave("outputs/Figure4_seasonality.png",
       fig_seas,
       width  = 4800,
       height = 1600,
       units  = "px",
       dpi    = 300)

# 8. Save objects for additional analyses ----
saveRDS(bio.matrix_test_regions_pred, "data/bio.matrix_test_regions_pred.rds")
saveRDS(selected_traps,               "data/selected_traps.rds")
saveRDS(env.matrix,                   "data/env.matrix_pred.rds")
saveRDS(seas.av.df,                   "data/seas.av.df.rds")
saveRDS(seas.df,                      "data/seas.df.rds")
saveRDS(trap_intercepts_df,           "data/trap_intercepts_df.rds")
saveRDS(lp_fixed_cal,                 "data/lp_fixed_cal.rds")

cat("All objects saved\n")
