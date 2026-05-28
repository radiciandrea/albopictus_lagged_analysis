
# define a starting season indicator also for observations
library(dplyr)
library(ggplot2)

# Load objects ----
bio.matrix_sel               <- readRDS("data/bio.matrix_sel_pred.rds") %>%
  select(c("ID", "year", "week", "medianTweek", "cumPrecweek", "eggs", "pred"))
bio.matrix_test_regions <- readRDS("data/bio.matrix_test_regions_pred.rds") %>%
  select(c("ID", "year", "week", "medianTweek", "cumPrecweek", "eggs", "pred"))
bio.matrix_test_years        <- readRDS("data/bio.matrix_test_years_pred.rds") %>%
  select(c("ID", "year", "week", "medianTweek", "cumPrecweek", "eggs", "pred"))

# there are duplicates in bio.matrix_test_regions_pred (ID 981, y 2019)
bio.matrix_test_regions_pred <- bio.matrix_test_regions_pred %>%
  filter(ID != 981)

bio.matrix <- rbind(bio.matrix_sel, bio.matrix_test_regions, bio.matrix_test_years)

# with temperature treshold only ----

predicted_cross_sel <- bio.matrix  %>%
  group_by(ID) %>% #group_by(ID, year) %>%
  summarise(threshold_thermal = week[which(medianTweek >= 17.5)[1]]) %>%
  filter(!is.na(threshold_thermal)) %>%
  ungroup()

observed_cross_sel <- bio.matrix %>%
  filter(!is.na(eggs)) %>%
  filter(week < 32) %>% # to exclude autumn monitoring noise
  group_by(ID, year) %>%
  mutate(consec = 1*(week == lag(week, default = 1)+1)) %>%
  mutate(egg_prec = 1*(lag(eggs, default =0))) %>%
  mutate(threshold_eggs = week*(consec == 1)*(egg_prec < 10)*(eggs > 10)) %>%
  filter(threshold_eggs >0) %>%
  summarise(threshold_eggs = min(threshold_eggs)) %>% # if more than one in the same year
  group_by(ID) %>% # to summarize per ID
  summarise(threshold_eggs = mean(threshold_eggs)) %>% # to summarize per ID
  ungroup()

# observed_predicted_cross = left_join(observed_cross_sel, predicted_cross_sel, by = c("ID", "year")) %>%
#   filter(!is.na(threshold_thermal))

observed_predicted_cross = left_join(observed_cross_sel, predicted_cross_sel, by = c("ID")) %>%
  filter(!is.na(threshold_thermal))

#plot e corr

cor_result <- cor.test(observed_predicted_cross$threshold_eggs,
                       observed_predicted_cross$threshold_thermal,
                       method = "pearson")

rmse = sqrt(mean((observed_predicted_cross$threshold_thermal - observed_predicted_cross$threshold_eggs )^2))
mae = mean(abs(observed_predicted_cross$threshold_thermal - observed_predicted_cross$threshold_eggs ))

cat("Pearson r (onset vs threshold crossing):", round(cor_result$estimate, 3), "\n")
cat("p-value:                                ", format(cor_result$p.value, scientific = TRUE), "\n")
cat("rmse:                                ", format(rmse, 3), "\n")
cat("mae:                                ", format(mae, 3), "\n")

ggplot(observed_predicted_cross,
             aes(x = threshold_thermal, y = threshold_eggs)) +
  geom_point(alpha = 0.2, size = 1, col = "orange") + #
  geom_smooth(method = "lm", col = "firebrick", linewidth = 1) +
  labs(x     = "Week of 17.5°C threshold crossing",
       y     = "Observed season onset (week)",
       title = "Thermal threshold vs observed season onset") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave( "outputs/Fig_observed_onset_vs_thermal_threshold_onset.png", width = 7, height = 6, dpi = 300)

# ### test esempio
# r = 1
# observed_predicted_cross[r,]
# plot(bio.matrix %>% filter(year == observed_predicted_cross$year[r],
#                                ID == observed_predicted_cross$ID[r]) %>% pull(eggs))


# only predicted (with the whole model) vs observed ----

predicted_cross_sel <- bio.matrix  %>%
  group_by(ID) %>% #group_by(ID, year) %>%
  summarise(model_predicted = week[which(pred > 10)[1]]) %>%
  filter(!is.na(model_predicted)) %>%
  ungroup()

observed_predicted_cross = left_join(observed_cross_sel, predicted_cross_sel, by = c("ID")) %>%
  filter(!is.na(model_predicted))

#plot e corr

cor_result <- cor.test(observed_predicted_cross$model_predicted,
                       observed_predicted_cross$threshold_eggs,
                       method = "pearson")

rmse = sqrt(mean((observed_predicted_cross$model_predicted - observed_predicted_cross$threshold_eggs )^2))
mae = mean(abs(observed_predicted_cross$model_predicted - observed_predicted_cross$threshold_eggs ))

cat("Pearson r (onset vs threshold crossing):", round(cor_result$estimate, 3), "\n")
cat("p-value:                                ", format(cor_result$p.value, scientific = TRUE), "\n")
cat("rmse:                                ", format(rmse, 3), "\n")
cat("mae:                                ", format(mae, 3), "\n")

ggplot(observed_predicted_cross,
       aes(x = model_predicted, y = threshold_eggs)) +
  geom_point(alpha = 0.2, size = 1, col = "orange") + #
  geom_smooth(method = "lm", col = "firebrick", linewidth = 1) +
  labs(x     = "Predicted season onset (week)",
       y     = "Observed season onset (week)",
       title = "predicted vs observed season onset") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())


ggsave( "outputs/Fig_observed_onset_vs_modelled_onset.png", width = 7, height = 6, dpi = 300)

# # with precipitation only ----
# 
# predicted_cross_sel <- bio.matrix  %>%
#   group_by(ID) %>% #group_by(ID, year) %>%
#   summarise(threshold_rain = week[which(cumPrecweek > 12)[1]]) %>%
#   filter(!is.na(threshold_rain)) %>%
#   ungroup()
# 
# observed_predicted_cross = left_join(observed_cross_sel, predicted_cross_sel, by = c("ID")) %>%
#   filter(!is.na(threshold_rain))
# 
# #plot e corr
# 
# cor_result <- cor.test(observed_predicted_cross$threshold_eggs,
#                        observed_predicted_cross$threshold_rain,
#                        method = "pearson")
# 
# # significativa e negativa, pensa te
# 
# cat("Pearson r (onset vs threshold crossing):", round(cor_result$estimate, 3), "\n")
# cat("p-value:                                ", format(cor_result$p.value, scientific = TRUE), "\n")
# 
# ggplot(observed_predicted_cross,
#        aes(x = threshold_rain, y = threshold_eggs)) +
#   geom_point(alpha = 0.2, size = 1, col = "orange") + #
#   geom_smooth(method = "lm", col = "firebrick", linewidth = 1) +
#   labs(x     = "Week of 17.5°C threshold crossing",
#        y     = "Observed season onset (week)",
#        title = "rain threshold vs observed season onset") +
#   theme_minimal(base_size = 11) +
#   theme(panel.grid.minor = element_blank())
