# Prediction of the start and end of the season based on the model and a threshold of 17.5 °C
# And comparison with observations

# library ----
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
bio.matrix_test_regions <- bio.matrix_test_regions %>%
  filter(ID != 981)

bio.matrix <- rbind(bio.matrix_sel, bio.matrix_test_regions, bio.matrix_test_years)

# ONSET ----

## temperature treshold only ----

thermal_onset_df <- bio.matrix  %>%
  group_by(ID) %>% #group_by(ID, year) %>%
  summarise(thermal_onset = week[which(medianTweek >= 17.5)[1]]) %>%
  filter(!is.na(thermal_onset)) %>%
  ungroup()

observed_onset_df <- bio.matrix %>%
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

observed_thermal_onset_df = left_join(observed_onset_df, thermal_onset_df, by = c("ID")) %>%
  filter(!is.na(thermal_onset))

#plot e corr

cor_result_thermal_onset <- cor.test(observed_thermal_onset_df$threshold_eggs,
                       observed_thermal_onset_df$thermal_onset,
                       method = "pearson")

rmse_thermal_onset = sqrt(mean((observed_thermal_onset_df$thermal_onset - observed_thermal_onset_df$threshold_eggs )^2))
mae_thermal_onset = mean(abs(observed_thermal_onset_df$thermal_onset - observed_thermal_onset_df$threshold_eggs ))

cat("Pearson r (onset vs threshold crossing):", round(cor_result_thermal_onset$estimate, 3), "\n")
cat("p-value:                                ", format(cor_result_thermal_onset$p.value, scientific = TRUE), "\n")
cat("rmse:                                ", format(rmse_thermal_onset, digits = 3), "\n")
cat("mae:                                ", format(mae_thermal_onset, digits =3), "\n")

p1 <- ggplot(observed_thermal_onset_df,
             aes(x = thermal_onset, y = threshold_eggs)) +
  geom_point(alpha = 0.2, size = 1, col = "orange") + #
  geom_smooth(method = "lm", col = "firebrick", linewidth = 1) +
  labs(x     = "Week of 17.5°C threshold crossing",
       y     = "Observed season onset (week)",
       title = "a) Thermal threshold vs observed onset") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank()) + 
  annotate("text", x = 33, y = 18, hjust = 1, label = paste0("Pearson's r: ", 
                                                  round(cor_result_thermal_onset$estimate, 3),
                                                  "\nMAE: ", format(mae_thermal_onset, digits = 3),
                                                  "\nRMSE: ", format(rmse_thermal_onset, digits = 3)))+ 
  theme(text=element_text(family="sans"))

p1

## predicted with the whole model ----

modelled_onset_df <- bio.matrix  %>%
  group_by(ID) %>% #group_by(ID, year) %>%
  summarise(model_onset = week[which(pred > 10)[1]]) %>%
  filter(!is.na(model_onset)) %>%
  ungroup()

observed_modelled_onset_df = left_join(observed_onset_df, modelled_onset_df, by = c("ID")) %>%
  filter(!is.na(model_onset))

#plot e corr

cor_result_modelled_onset <- cor.test(observed_modelled_onset_df$threshold_eggs,
                                     observed_modelled_onset_df$model_onset,
                                     method = "pearson")

rmse_modelled_onset = sqrt(mean((observed_modelled_onset_df$model_onset - observed_modelled_onset_df$threshold_eggs )^2))
mae_modelled_onset = mean(abs(observed_modelled_onset_df$model_onset - observed_modelled_onset_df$threshold_eggs ))

cat("Pearson r (onset vs threshold crossing):", round(cor_result_modelled_onset$estimate, 3), "\n")
cat("p-value:                                ", format(cor_result_modelled_onset$p.value, scientific = TRUE), "\n")
cat("rmse:                                ", format(rmse_modelled_onset, digits = 3), "\n")
cat("mae:                                ", format(mae_modelled_onset, digits =3), "\n")

p2 <- ggplot(observed_modelled_onset_df,
             aes(x = model_onset, y = threshold_eggs)) +
  geom_point(alpha = 0.2, size = 1, col = "orange") + #
  geom_smooth(method = "lm", col = "firebrick", linewidth = 1) +
  labs(x     = "T.P.seas.year-predicted onset (week)",
       y     = "Observed season onset (week)",
       title = "b) Full-model predicted vs observed onset") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank()) + 
  annotate("text", x = 33, y = 18, hjust = 1, label = paste0("Pearson's r: ", 
                                                             round(cor_result_modelled_onset$estimate, 3),
                                                             "\nMAE: ", format(mae_modelled_onset, digits = 3),
                                                             "\nRMSE: ", format(rmse_modelled_onset, digits = 3)))+ 
  theme(text=element_text(family="sans"))

p1+p2

# OFFSET ----

## temperature treshold only ----

thermal_offset_df <- bio.matrix  %>%
  filter(week > 32) %>%
  group_by(ID) %>% #group_by(ID, year) %>%
  summarise(thermal_offset = week[which(medianTweek <= 17.5)[1]]) %>%
  filter(!is.na(thermal_offset)) %>%
  ungroup()

observed_offset_df <- bio.matrix %>%
  filter(!is.na(eggs)) %>%
  filter(week > 32) %>% # to exclude autumn monitoring noise
  group_by(ID, year) %>%
  mutate(consec = 1*(week == lag(week, default = 1)+1)) %>%
  mutate(egg_prec = 1*(lag(eggs, default =0))) %>%
  mutate(threshold_eggs = week*(consec == 1)*(egg_prec > 10)*(eggs < 10)) %>%
  filter(threshold_eggs >0) %>%
  summarise(threshold_eggs = min(threshold_eggs)) %>% # if more than one in the same year
  group_by(ID) %>% # to summarize per ID
  summarise(threshold_eggs = mean(threshold_eggs)) %>% # to summarize per ID
  ungroup()

observed_thermal_offset_df = left_join(observed_offset_df, thermal_offset_df, by = c("ID")) %>%
  filter(!is.na(thermal_offset))

#plot e corr

cor_result_thermal_offset <- cor.test(observed_thermal_offset_df$threshold_eggs,
                                      observed_thermal_offset_df$thermal_offset,
                                      method = "pearson")

rmse_thermal_offset = sqrt(mean((observed_thermal_offset_df$thermal_offset - observed_thermal_offset_df$threshold_eggs )^2))
mae_thermal_offset = mean(abs(observed_thermal_offset_df$thermal_offset - observed_thermal_offset_df$threshold_eggs ))

cat("Pearson r (offset vs threshold crossing):", round(cor_result_thermal_offset$estimate, 3), "\n")
cat("p-value:                                ", format(cor_result_thermal_offset$p.value, scientific = TRUE), "\n")
cat("rmse:                                ", format(rmse_thermal_offset, digits = 3), "\n")
cat("mae:                                ", format(mae_thermal_offset, digits =3), "\n")

p3 <- ggplot(observed_thermal_offset_df,
             aes(x = thermal_offset, y = threshold_eggs)) +
  geom_point(alpha = 0.2, size = 1, col = "steelblue") + #
  geom_smooth(method = "lm", col = "firebrick", linewidth = 1) +
  labs(x     = "Week of 17.5°C threshold crossing",
       y     = "Observed season offset (week)",
       title = "c) Thermal threshold vs observed offset") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank()) + 
  annotate("text", x = 33, y = 45, hjust = 0, label = paste0("Pearson's r: ", 
                                                             round(cor_result_thermal_offset$estimate, 3),
                                                             "\nMAE: ", format(mae_thermal_offset, digits = 3),
                                                             "\nRMSE: ", format(rmse_thermal_offset, digits = 3)))+ 
  theme(text=element_text(family="sans"))

p3

## predicted with the whole model ----

modelled_offset_df <- bio.matrix  %>%
  filter(week > 32) %>%
  group_by(ID) %>% #group_by(ID, year) %>%
  summarise(model_offset = week[which(pred < 10)[1]]) %>%
  filter(!is.na(model_offset)) %>%
  ungroup()

observed_modelled_offset_df = left_join(observed_offset_df, modelled_offset_df, by = c("ID")) %>%
  filter(!is.na(model_offset))

#plot e corr

cor_result_modelled_offset <- cor.test(observed_modelled_offset_df$threshold_eggs,
                                       observed_modelled_offset_df$model_offset,
                                       method = "pearson")

rmse_modelled_offset = sqrt(mean((observed_modelled_offset_df$model_offset - observed_modelled_offset_df$threshold_eggs )^2))
mae_modelled_offset = mean(abs(observed_modelled_offset_df$model_offset - observed_modelled_offset_df$threshold_eggs ))

cat("Pearson r (offset vs threshold crossing):", round(cor_result_modelled_offset$estimate, 3), "\n")
cat("p-value:                                ", format(cor_result_modelled_offset$p.value, scientific = TRUE), "\n")
cat("rmse:                                ", format(rmse_modelled_offset, digits = 3), "\n")
cat("mae:                                ", format(mae_modelled_offset, digits =3), "\n")

p4 <- ggplot(observed_modelled_offset_df,
             aes(x = model_offset, y = threshold_eggs)) +
  geom_point(alpha = 0.2, size = 1, col = "steelblue") + #
  geom_smooth(method = "lm", col = "firebrick", linewidth = 1) +
  labs(x     = "T.P.seas.year-predicted offeset (week)",
       y     = "Observed season offset (week)",
       title = "d) Full-model predicted vs observed offset") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank()) + 
  annotate("text", x = 33, y = 45, hjust = 0, label = paste0("Pearson's r: ", 
                                                             round(cor_result_modelled_offset$estimate, 3),
                                                             "\nMAE: ", format(mae_modelled_offset, digits = 3),
                                                             "\nRMSE: ", format(rmse_modelled_offset, digits = 3)))+ 
  theme(text=element_text(family="sans"))

fig_4 = p1 + p2 + p3 + p4 


ggsave("outputs/Figure3.png",
       fig_4,
       width  = 9.5,
       height = 8,
       dpi    = 300)








# with precipitation threshold only ----
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
