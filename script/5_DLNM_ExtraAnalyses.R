# Libraries ----
library(tidyverse)
library(lubridate)
library(ggplot2)
library(patchwork)

# Load objects ----
bio.matrix_sel               <- readRDS("data/bio.matrix_cal.rds")
bio.matrix_test_regions_pred <- readRDS("data/bio.matrix_test_regions_pred.rds")
selected_traps               <- readRDS("data/selected_traps.rds")
env.matrix                   <- readRDS("data/env.matrix_pred.rds")
seas.av.df                   <- readRDS("data/seas.av.df.rds")
seas.df                      <- readRDS("data/seas.df.rds")
trap_intercepts_df           <- readRDS("data/trap_intercepts_df.rds")
lp_fixed_cal                 <- readRDS("data/lp_fixed_cal.rds")

tknots    <- unname(quantile(bio.matrix_sel$medianTweek,
                             probs = c(0.25, 0.75), na.rm = TRUE))
pknots    <- unname(quantile(bio.matrix_sel$cumPrecweek,
                             probs = c(0.25, 0.75, 0.9), na.rm = TRUE))
eggThresh <- 10

# 1. Seasonal anticipation bias ----
anticipation_df <- bio.matrix_test_regions_pred %>%
  dplyr::filter(ID %in% selected_traps$ID, !is.na(eggs), is.finite(pred)) %>%
  dplyr::group_by(Region, ID, year) %>%
  dplyr::summarise(
    obs_peak_week  = week[which.max(eggs)],
    pred_peak_week = week[which.max(pred)],
    bias           = pred_peak_week - obs_peak_week,
    .groups        = "drop"
  ) %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    n_years   = dplyr::n(),
    mean_bias = round(mean(bias, na.rm = TRUE), 1),
    sd_bias   = round(sd(bias,   na.rm = TRUE), 1),
    .groups   = "drop"
  )

cat("Seasonal anticipation bias (positive = model predicts peak before observed):\n")
as.data.frame(anticipation_df)

# Overall summary
cat("Overall mean bias (weeks):", round(mean(anticipation_df$mean_bias, na.rm = TRUE), 1), "\n")
cat("Overall sd   bias (weeks):", round(sd(anticipation_df$mean_bias,   na.rm = TRUE), 1), "\n")

# 2. Thermal threshold crossing vs season onset ----
threshold_cross <- env.matrix %>%
  dplyr::group_by(x, y, year) %>%
  dplyr::arrange(week) %>%
  dplyr::summarise(
    threshold_week = week[which(medianTweek >= 17.5)[1]],
    .groups        = "drop"
  )

threshold_av <- threshold_cross %>%
  dplyr::group_by(x, y) %>%
  dplyr::summarise(
    mean_threshold_week = mean(threshold_week, na.rm = TRUE),
    .groups             = "drop"
  )

onset_vs_threshold <- dplyr::left_join(
  seas.av.df %>% dplyr::select(x, y, m.start),
  threshold_av,
  by = c("x", "y")
) %>% dplyr::filter(!is.na(m.start), !is.na(mean_threshold_week))

cor_result <- cor.test(onset_vs_threshold$m.start,
                       onset_vs_threshold$mean_threshold_week,
                       method = "pearson")

cat("Pearson r (onset vs threshold crossing):", round(cor_result$estimate, 3), "\n")
cat("p-value:                                ", format(cor_result$p.value, scientific = TRUE), "\n")

pt <- ggplot(onset_vs_threshold,
       aes(x = mean_threshold_week, y = m.start)) +
  geom_point(alpha = 0.1, size = 0.5, col = "steelblue") +
  geom_smooth(method = "lm", col = "firebrick", linewidth = 1) +
  labs(x     = "Week of 17.5°C threshold crossing",
       y     = "Predicted season onset (week)",
       title = "Thermal threshold vs predicted season onset") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave( "outputs/Fig_threshold_vs_onset.png", pt,  width = 7, height = 6, dpi = 300)

# 3. Interannual variability in season length ----
seas_annual_summary <- seas.df %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    n_cells       = dplyr::n(),
    mean_length   = round(mean(end - start, na.rm = TRUE) / 4, 2),
    median_length = round(median(end - start, na.rm = TRUE) / 4, 2),
    sd_length     = round(sd(end - start, na.rm = TRUE) / 4, 2),
    mean_onset    = round(mean(start, na.rm = TRUE), 1),
    mean_offset   = round(mean(end,   na.rm = TRUE), 1),
    .groups       = "drop"
  )

cat("Interannual variability in season length:\n")
as.data.frame(seas_annual_summary)

seas.df %>%
  dplyr::mutate(length_months = (end - start) / 4) %>%
  ggplot(aes(x = factor(year), y = length_months)) +
  geom_boxplot(fill = "steelblue", alpha = 0.5, outlier.size = 0.5) +
  labs(x     = "Year",
       y     = "Season length (months)",
       title = "Interannual variability in predicted season length") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave("outputs/Fig_interannual_season_length.png", width = 7, height = 5, dpi = 300)

# 4. Precipitation threshold — fraction of active weeks above 12 mm by region ----
# Compute calibration predictions using trap-specific intercepts
cal_pred_df <- bio.matrix_sel %>%
  dplyr::left_join(trap_intercepts_df, by = "ID") %>%
  dplyr::mutate(
    pred_cal = exp(lp_fixed_cal + intercept)
  ) %>%
  dplyr::filter(!is.na(pred_cal), pred_cal > eggThresh)

precip_threshold <- cal_pred_df %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    n_active_weeks = dplyr::n(),
    pct_above_12mm = round(mean(cumPrecweek >= 12, na.rm = TRUE) * 100, 1),
    median_precip  = round(median(cumPrecweek, na.rm = TRUE), 1),
    .groups        = "drop"
  ) %>%
  dplyr::arrange(pct_above_12mm)

cat("Precipitation during active season by region:\n")
as.data.frame(precip_threshold)

precip_threshold <- precip_threshold %>%
  ggplot(aes(x = reorder(Region, pct_above_12mm), y = pct_above_12mm)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_hline(yintercept = 50, linetype = "dashed", col = "grey40") +
  coord_flip() +
  labs(x     = NULL,
       y     = "% of active-season weeks with precip > 12 mm",
       title = "Precipitation threshold during active oviposition season") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave("outputs/Fig_precip_threshold_by_region.png",precip_threshold,  width = 8, height = 5, dpi = 300)


# Combine and save  
fig_seas <- pt + precip_threshold
fig_seas
ggsave("outputs/Figure5.png",
       fig_seas,
       width  = 8,
       height = 14,
       dpi    = 300)
