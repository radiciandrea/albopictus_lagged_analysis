# Best model: T.P.seas (temperature + precipitation + seasonality)

# Libraries 
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
library(sf)
library(terra)
library(mapview)
library(tidyterra)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(gcplyr)

# 1. Load data ----
bio.matrix_sel <- readRDS("data/bio.matrix_cal.rds")
str(bio.matrix_sel)

unique_locs <- bio.matrix_sel %>%
  ungroup() %>%
  distinct(ID, lon, lat)

loc_vect <- vect(unique_locs,
                 geom = c("lon", "lat"),
                 crs  = "EPSG:4326")
mapview(st_as_sf(loc_vect))

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

# 3. Fit T.P.seas ----
mod_T.P.seas <- gnm(eggs ~ cbTemp + cbPrec + bseas + year,
                    eliminate = as.factor(ID),
                    data      = bio.matrix_sel,
                    family    = quasipoisson(link = "log"),
                    na.action = "na.exclude")

cat("Converged:", mod_T.P.seas$converged, "\n")
cat("Dispersion:", summary(mod_T.P.seas)$dispersion, "\n")
cat("Residual deviance:", summary(mod_T.P.seas)$deviance, "\n")

# Null model for context
mod_null <- gnm(eggs ~ 1,
                eliminate = as.factor(ID),
                data      = bio.matrix_sel,
                family    = quasipoisson(link = "log"),
                na.action = "na.exclude")

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

step = 0.1
temp_vec = seq(6, 30, step)
cptmean <- crosspred(cbTemp, mod_T.P.seas, cen = myCen_T, by = step, at = temp_vec)

# max of the curve
max(cptmean$allRRfit)
temp_vec[which(cptmean$allRRfit == max(cptmean$allRRfit))]

# doubling temperature - tolerance of 0.02
temp_vec[which(abs(cptmean$allRRfit - 2) < 0.02)][1]

# derivative of the curve
allRRfit_deriv = calc_deriv(y = cptmean$allRRfit, x = temp_vec)

# and find the maximum (which is the temperature at the infelction point)
temp_vec[which(allRRfit_deriv == max(allRRfit_deriv, na.rm = T))]

# Overall cumulative effect
# png("outputs/temp_overall.png", width = 800, height = 600, res = 120)
plot(cptmean, ptype = "overall",
     xlab = "Temperature (°C)", ylab = "IRR",
     col  = 2, main = "a ) Temperature — overall cumulative effect")
# dev.off()

# Lag profiles at fixed temperatures
# png("outputs/temp_lag_profiles.png", width = 900, height = 800, res = 120)
# layout(matrix(1:4, ncol = 2, nrow = 2))
# for (tmp in c(6, 14, 22, 30)) {
#   plot(cptmean, var = tmp,
#        xlab = "Lag (weeks)", ylab = "IRR", col = 2,
#        main = paste("T =", tmp, "°C"), ylim = c(0, 3))
# }
# layout(1)
# dev.off()

# Exposure-response at fixed lags
# png("outputs/temp_exposure_response.png", width = 900, height = 800, res = 120)
# layout(matrix(1:4, ncol = 2, nrow = 2))
# for (lg in c(0, 2, 4, 8)) {
#   plot(cptmean, lag = lg,
#        xlab = "Temperature (°C)", ylab = "IRR", col = 2,
#        main = paste("Lag =", lg, "weeks"), ylim = c(0, 3))
# }
# layout(1)
# dev.off()

plot(cptmean, ptype = "3d",
     xlab = "Temperature (°C)", ylab = "Lag (weeks)", zlab = "IRR",
     main = "Temperature lag-response surface")


# 5. Precipitation effects ----
myCen_P <- 0
P_max   <- unname(quantile(bio.matrix_sel$cumPrecweek, 0.95, na.rm = TRUE))
cat("Precipitation Q95:", P_max, "mm\n")

step = 0.1
prec_vec = seq(0, ceiling(P_max), step)

cppmean <- crosspred(cbPrec, mod_T.P.seas,
                     cen = myCen_P, by = step, at = prec_vec)
# intersection 0
id_prec = which(abs(cppmean$allRRfit-1)< 0.001)[2] # the first os trivially 0
prec_vec[id_prec]

# Overall cumulative effect
# png("outputs/precip_overall.png", width = 800, height = 600, res = 120)
plot(cppmean, ptype = "overall",
     xlab = "Precipitation (mm)", ylab = "IRR",
     col  = 2, main = "Precipitation — overall cumulative effect",
     ylim = c(0.5, 3))
# dev.off()

# Lag profiles at fixed precipitation values
# png("outputs/precip_lag_profiles.png", width = 900, height = 800, res = 120)
# layout(matrix(1:4, ncol = 2, nrow = 2))
# for (prc in c(5, 10, 20, 50)) {
#   plot(cppmean, var = prc,
#        xlab = "Lag (weeks)", ylab = "IRR", col = 2,
#        main = paste("Precip =", prc, "mm"), ylim = c(0.75, 1.25))
# }
# layout(1)
# dev.off()

# Exposure-response at fixed lags
# png("outputs/precip_exposure_response.png", width = 900, height = 800, res = 120)
# layout(matrix(1:4, ncol = 2, nrow = 2))
# for (lg in c(0, 2, 4, 8)) {
#   plot(cppmean, lag = lg,
#        xlab = "Precipitation (mm)", ylab = "IRR", col = 2,
#        main = paste("Lag =", lg, "weeks"), ylim = c(0, 2))
# }
# layout(1)
# dev.off()

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
mapview(st_as_sf(loc_vect))

bio.matrix_test_regions <- readRDS("data/bio.matrix_test_regions.rds") %>%
  dplyr::mutate(weekn = as.numeric(week))

unique_locs <- bio.matrix_test_regions %>%
  ungroup() %>%
  distinct(ID, lon, lat)

loc_vect <- vect(unique_locs,
                 geom = c("lon", "lat"),
                 crs  = "EPSG:4326")
mapview(st_as_sf(loc_vect))

cat("Test years   — rows:", nrow(bio.matrix_test_years),
    "| traps:", dplyr::n_distinct(bio.matrix_test_years$ID),
    "| regions:", dplyr::n_distinct(bio.matrix_test_years$Region), "\n")

cat("Test regions — rows:", nrow(bio.matrix_test_regions),
    "| traps:", dplyr::n_distinct(bio.matrix_test_regions$ID),
    "| regions:", dplyr::n_distinct(bio.matrix_test_regions$Region), "\n")

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
trap_intercepts_df <- bio.matrix_sel %>%
  dplyr::mutate(lp_fixed = lp_fixed_cal, lp_full = lp_full_cal) %>%
  dplyr::filter(!is.na(lp_full)) %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(
    intercept = mean(lp_full - lp_fixed, na.rm = TRUE),
    .groups   = "drop"
  )

# summary of the exponential of the intercepts (this is to test that it is 0)
summary(exp(trap_intercepts_df$intercept))

cat("Trap intercepts computed:", nrow(trap_intercepts_df), "\n")

# Population-average intercept for unknown trap IDs  
valid_idx     <- which(!is.na(lp_full_cal))
avg_intercept <- mean(lp_full_cal[valid_idx] - lp_fixed_cal[valid_idx], na.rm = TRUE)
cat("Population-average intercept:", round(avg_intercept, 4), "\n")

##6.2 Filter test_years to known trap IDs only ----
known   <- intersect(unique(bio.matrix_test_years$ID), unique(bio.matrix_sel$ID))
unknown <- setdiff(unique(bio.matrix_test_years$ID),   unique(bio.matrix_sel$ID))
cat("Test years — known intercept:", length(known),
    "| using avg:", length(unknown),
    "| unknown IDs:", paste(unknown, collapse = ", "), "\n")

bio.matrix_test_years_known <- bio.matrix_test_years %>%
  dplyr::filter(ID %in% known)

# Prediction function  
predict_with_intercept <- function(newdata, mod,
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
  
  bseas_new <- ns(newdata$weekn, df = 4)
  
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
  cat("  RMSLE:    ", round(sqrt(mean((log10(pred + 1) - log10(obs + 1))^2)), 4), "\n")
  cat("  Spearman's rank r:", round(cor(obs, pred, method = "spearman"),           4), "\n\n")
}

## 6.4 Run predictions ----
pred_years_known <- predict_with_intercept(
  bio.matrix_test_years_known, mod_T.P.seas,
  tknots, pknots, trap_intercepts_df, avg_intercept
)

pred_regions <- predict_with_intercept(
  bio.matrix_test_regions, mod_T.P.seas,
  tknots, pknots, trap_intercepts_df, avg_intercept
)

pred_cal <- predict_with_intercept(
  bio.matrix_sel, mod_T.P.seas,
  tknots, pknots, trap_intercepts_df, avg_intercept
)

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

eval_metrics(
  bio.matrix_test_years_known$eggs,
  pred_years_known,
  "Test set 1 "
)

eval_metrics(bio.matrix_test_regions$eggs,
             pred_regions,
             "Test set 2 — held-out regions")

## 6.6 Per-trap evaluation  ----
per_trap_test_years <- bio.matrix_test_years_known %>%
  dplyr::mutate(pred = pred_years_known) %>%
  dplyr::filter(!is.na(eggs), is.finite(pred)) %>%
  dplyr::group_by(ID, Region) %>%
  dplyr::summarise(
    n         = dplyr::n(),
    Spearman_r = cor(eggs, pred, method = "spearman", use = "complete.obs"),
    p_value_Sr    = round((cor.test(eggs, pred, method = "spearman"))$p.value,4),
    RMSLE     = sqrt(mean((log10(pred + 1) - log(eggs + 1))^2)),
    .groups   = "drop"
  ) %>%
  # dplyr::filter(n >= 10) %>%
  dplyr::mutate(test = "years")

per_trap_test_regions <- bio.matrix_test_regions %>%
  dplyr::mutate(pred = pred_regions) %>%
  dplyr::filter(!is.na(eggs), is.finite(pred)) %>%
  dplyr::group_by(ID, Region) %>%
  dplyr::summarise(
    n         = dplyr::n(),
    Spearman_r = cor(eggs, pred, method = "spearman", use = "complete.obs"),
    p_value_Sr    = round((cor.test(eggs, pred, method = "spearman"))$p.value,4),
    RMSLE     = sqrt(mean((log10(pred + 1) - log(eggs + 1))^2)),
    .groups   = "drop"
  ) %>%
  # dplyr::filter(n >= 10) %>%
  dplyr::mutate(test = "regions")

per_trap_r <- per_trap_test_years %>% #rbind(per_trap_test_year, per_trap_test_regions) %>%
  dplyr::arrange(desc(Spearman_r))


cat("Per-trap validation summary:\n")
cat("  Median Spearman's rank r:", round(median(per_trap_r$Spearman_r, na.rm = TRUE), 3), "\n")
cat("  IQR:             ", round(quantile(per_trap_r$Spearman_r, 0.25, na.rm = TRUE), 3),
    "—", round(quantile(per_trap_r$Spearman_r, 0.75, na.rm = TRUE), 3), "\n")
cat("  Median RMSLE:    ", round(median(per_trap_r$RMSLE, na.rm = TRUE), 3), "\n")
cat("  Traps r > 0.5:   ", sum(per_trap_r$Spearman_r > 0.5, na.rm = TRUE),
    "/", sum(!is.na(per_trap_r$Spearman_r)), "\n\n")
cat("  Significative traps (0.001):   ", sum(per_trap_r$p_value_Sr < 0.001, na.rm = TRUE),
    "/", sum(!is.na(per_trap_r$Spearman_r)), "\n\n")

## 6.7 Per-region summary — test_years ----
bio.matrix_test_years_known %>%
  dplyr::mutate(pred = pred_years_known) %>%
  dplyr::filter(!is.na(eggs), is.finite(pred)) %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    n         = dplyr::n(),
    Spearman_r = round(cor(eggs, pred, method = "spearman"), 4),
    p_value_Sr    = round((cor.test(eggs, pred, method = "spearman"))$p.value,4),
    RMSLE     = round(sqrt(mean((log10(pred + 1) - log(eggs + 1))^2)), 4),
    .groups   = "drop"
  ) %>%
  dplyr::arrange(-Spearman_r) %>%
  as.data.frame()

##6.8 Per-region breakdown — test_regions ----
bio.matrix_test_regions %>%
  dplyr::mutate(pred = pred_regions) %>%
  dplyr::filter(!is.na(eggs), is.finite(pred)) %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    n         = dplyr::n(),
    Spearman_r = round(cor(eggs, pred, method = "spearman"), 4),
    p_value_Sr    = round((cor.test(eggs, pred, method = "spearman"))$p.value, 4),
    RMSLE     = round(sqrt(mean((log10(pred + 1) - log(eggs + 1))^2)), 4),
    .groups   = "drop"
  ) %>%
  dplyr::arrange(-Spearman_r)

# Table 4: Per-region — temporal test set  
table4 <- bio.matrix_test_years_known %>%
  dplyr::mutate(pred = pred_years_known) %>%
  dplyr::filter(!is.na(eggs), is.finite(pred)) %>%
  dplyr::group_by(Region, Country) %>%
  dplyr::summarise(
    n_traps      = dplyr::n_distinct(ID),
    n_trap_weeks = dplyr::n(),
    Spearman_r    = round(cor(eggs, pred, method = "spearman"), 3),
    RMSLE        = round(sqrt(mean((log10(pred + 1) - log(eggs + 1))^2)), 3)
  ) %>%
  dplyr::arrange(-Spearman_r)

write.csv(table4, "outputs/table4_temporal_validation.csv", row.names = FALSE)
print(table4)

# Table 5: Per-region — spatial test set  
table5 <- bio.matrix_test_regions %>%
  dplyr::mutate(pred = pred_regions) %>%
  dplyr::filter(!is.na(eggs), is.finite(pred)) %>%
  dplyr::group_by(Region, Country) %>%
  dplyr::summarise(
    n_traps      = dplyr::n_distinct(ID),
    n_trap_weeks = dplyr::n(),
    Spearman_r    = round(cor(eggs, pred, method = "spearman"), 3),
    RMSLE        = round(sqrt(mean((log10(pred + 1) - log(eggs + 1))^2)), 3),
  ) %>%
  dplyr::arrange(-Spearman_r)

write.csv(table5, "outputs/table5_spatial_validation.csv", row.names = FALSE)
print(table5)

## 6.9 Visual check — example trap ----
ID_check <- 3769 #4228

check_df <- bio.matrix_test_years_known %>%
  dplyr::filter(ID == ID_check) %>%
  dplyr::mutate(pred = pred_years_known[bio.matrix_test_years_known$ID == ID_check])

plot(check_df$date, check_df$eggs,
     pch = 16, xlab = "Date", ylab = "Eggs",
     main = paste("ID", ID_check, "— test years"))
lines(check_df$date, check_df$pred, col = "darkgreen", lwd = 2)
legend("topright", legend = c("Observed", "Predicted"),
       col = c("black", "darkgreen"), pch = c(16, NA), lty = c(NA, 1))

cat("Spearman r ID", ID_check, ":",
    round(cor(check_df$eggs, check_df$pred,
              use = "complete.obs", method = "spearman"), 4), "\n")


## 6.10 Visual check — spatial test set, best trap per region ----

# Add predictions to test_regions
bio.matrix_test_regions_pred <- bio.matrix_test_regions %>%
  dplyr::mutate(pred = pred_regions)

# Find best trap per region by Spearman r
best_traps <- bio.matrix_test_regions_pred %>%
  dplyr::filter(!is.na(eggs), is.finite(pred)) %>%
  dplyr::group_by(Region, ID) %>%
  dplyr::summarise(
    n         = dplyr::n(),
    Spearman_r = cor(eggs, pred, method = "spearman", use = "complete.obs"),
    .groups   = "drop"
  ) %>%
  dplyr::filter(n >= 10) %>%
  dplyr::group_by(Region) %>%
  dplyr::slice_max(Spearman_r, n = 1) %>%
  dplyr::ungroup()

cat("Best trap per region:\n")
print(as.data.frame(best_traps))

# representative traps from spatial test set — low to high Spearman r
selected_traps <- data.frame(
  Region    = c("Sicily", "Lezhe", "Fier", "Lushnje", "Cote Azur", "Région lémanique"),
  ID        = c(22983,    14499,   17040,  16578,     9987,        4208),
  Spearman_r = c(0.36,     0.48,    0.58,   0.69,      0.78,        0.96)
)

cat("Selected traps:\n")
print(selected_traps)

# Plot
png("outputs/Fig_spatial_validation_6traps.png",
    width = 4800, height = 3200, res = 200)

par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

for (i in seq_len(nrow(selected_traps))) {
  
  id_check  <- selected_traps$ID[i]
  reg       <- selected_traps$Region[i]
  r_val     <- selected_traps$Spearman_r[i]
  
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
cat("Season length range (weeks):", range(seas.av.df$s.length, na.rm = TRUE), "\n")
cat("Onset range (week):         ", range(seas.av.df$m.start,  na.rm = TRUE), "\n")
cat("Offset range (week):        ", range(seas.av.df$m.end,    na.rm = TRUE), "\n")
 
# Common theme  
map_theme <- theme_minimal(base_size = 11) +
  theme(
    axis.title        = element_blank(),
    axis.text         = element_text(size = 8, colour = "grey40"),
    panel.grid        = element_line(colour = "grey92", linewidth = 0.3),
    legend.position   = "bottom",
    legend.key.width  = unit(1.2, "cm"),
    legend.key.height = unit(0.4, "cm"),
    legend.title      = element_text(size = 9),
    legend.text       = element_text(size = 8),
    plot.title        = element_text(size = 11, face = "bold", hjust = 0),
    plot.caption      = element_text(size = 7, colour = "grey50", hjust = 0)
  )

threshold_note <- "Season defined as weeks with predicted egg count > 10 eggs per trap per week. Average over calibration period 2018–2023."

# Panel A: season length (months)  
pA <- ggplot() +
  geom_spatraster(data = slength / 4) +
  scale_fill_stepsn(
    name     = "Season length (months)",
    colours  = brewer.pal(9, "YlOrRd"),
    breaks   = seq(0, 9, 1),
    limits   = c(0, 9),
    na.value = "white",
    guide    = guide_colorsteps(
      show.limits = TRUE,
      barwidth    = unit(12, "cm"),
      barheight   = unit(0.4, "cm")
    )
  ) +
  coord_equal() +
  map_theme+
  labs(title = "A — Average season length")

# Panel B: season onset (week)  
pB <- ggplot() +
  geom_spatraster(data = onset) +
  scale_fill_stepsn(
    name     = "Season onset (ISO week)",
    colours  = brewer.pal(9, "Blues"),
    breaks   = seq(13, 32, 3),
    limits   = c(13, 32),
    na.value = "white",
    guide    = guide_colorsteps(
      show.limits = TRUE,
      barwidth    = unit(12, "cm"),
      barheight   = unit(0.4, "cm")
    )
  ) +
  coord_equal() +
  labs(title = "B — Average season onset") +
  map_theme

# Panel C: season offset (week)  
pC <- ggplot() +
  geom_spatraster(data = offset) +
  scale_fill_stepsn(
    name     = "Season offset (ISO week)",
    colours  = brewer.pal(9, "Greens"),
    breaks   = seq(33, 49, 2),
    limits   = c(33, 49),
    na.value = "white",
    guide    = guide_colorsteps(
      show.limits = TRUE,
      barwidth    = unit(12, "cm"),
      barheight   = unit(0.4, "cm")
    )
  ) +
  coord_equal() +
  labs(title = "C — Average season offset",  caption = threshold_note) +
  map_theme

# Combine and save  
fig_seas <- pA / pB / pC
fig_seas
ggsave("outputs/Figure4_seasonality.png",
       fig_seas,
       width  = 8,
       height = 14,
       dpi    = 300)


# 8. Save objects for additional analyses ----

# Add predictions to test_years
bio.matrix_test_years_pred <- bio.matrix_test_years_known %>%
  dplyr::mutate(pred = pred_years_known)

# Add predictions to test_regions
bio.matrix_sel_pred <- bio.matrix_sel %>%
  dplyr::mutate(pred = pred_cal)

saveRDS(bio.matrix_test_regions_pred, "data/bio.matrix_test_regions_pred.rds")
saveRDS(bio.matrix_test_years_pred, "data/bio.matrix_test_years_pred.rds")
saveRDS(bio.matrix_sel_pred, "data/bio.matrix_sel_pred.rds")
saveRDS(selected_traps,               "data/selected_traps.rds")
saveRDS(env.matrix,                   "data/env.matrix_pred.rds")
saveRDS(seas.av.df,                   "data/seas.av.df.rds")
saveRDS(seas.df,                      "data/seas.df.rds")
saveRDS(trap_intercepts_df,           "data/trap_intercepts_df.rds")
saveRDS(lp_fixed_cal,                 "data/lp_fixed_cal.rds")

cat("All objects saved\n")
 