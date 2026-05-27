
# define a starting season indicator also for observations
library(dplyr)
library(ggplot2)

# Load objects ----
bio.matrix_sel               <- readRDS("data/bio.matrix_cal.rds") %>%
  select(c("ID", "year", "week", "medianTweek", "eggs"))
bio.matrix_test_regions_pred <- readRDS("data/bio.matrix_test_regions_pred.rds") %>%
  select(c("ID", "year", "week", "medianTweek", "eggs"))
bio.matrix_test_years        <- readRDS("data/bio.matrix_test_years.rds") %>%
  select(c("ID", "year", "week", "medianTweek", "eggs"))

# there are duplicates in bio.matrix_test_regions_pred (ID 981, y 2019)
bio.matrix_test_regions_pred <- bio.matrix_test_regions_pred %>%
  filter(ID != 981)

bio.matrix <- rbind(bio.matrix_sel, bio.matrix_test_regions_pred, bio.matrix_test_years)

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

cat("Pearson r (onset vs threshold crossing):", round(cor_result$estimate, 3), "\n")
cat("p-value:                                ", format(cor_result$p.value, scientific = TRUE), "\n")

ggplot(observed_predicted_cross,
             aes(x = threshold_thermal, y = threshold_eggs)) +
  geom_point(alpha = 0.2, size = 1, col = "orange") + #
  geom_smooth(method = "lm", col = "firebrick", linewidth = 1) +
  labs(x     = "Week of 17.5°C threshold crossing",
       y     = "Observed season onset (week)",
       title = "Thermal threshold vs observed season onset") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave( "outputs/Fig_threshold_vs_observed_onset.png", width = 7, height = 6, dpi = 300)

### test esempio
r = 1
observed_predicted_cross[r,]
plot(bio.matrix %>% filter(year == observed_predicted_cross$year[r],
                               ID == observed_predicted_cross$ID[r]) %>% pull(eggs))

