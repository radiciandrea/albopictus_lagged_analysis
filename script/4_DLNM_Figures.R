while (!is.null(dev.list())) dev.off()

width_set = 1600
height_set = width_set*8/9
res_set = 300

# Temperature effects

png("outputs/Fig2A_temp_overall.png", width = width_set, height = height_set, res = res_set)
# par(mar = c(5, 5, 4, 2))
plot(cptmean, ptype = "overall",
     xlab = "Temperature (°C)", ylab = "IRR",
     col = "firebrick", lwd = 2,
     main = "A — Overall cumulative effect")
abline(h = 1, col = "grey40", lwd = 1, lty = 2)
dev.off()

codeLetters = c("B", "C", "D", "E")

png("outputs/Fig2BCDE_temp_exposure_lags.png", width = width_set, height = height_set, res = res_set)
# par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
par(mfrow = c(2, 2))
cols <- c("steelblue", "darkgreen", "firebrick", "darkorange")
lags <- c(0, 2, 4, 8)
for (j in seq_along(lags)) {
  plot(cptmean, lag = lags[j],
       xlab = "Temperature (°C)", ylab = "IRR",
       col = cols[j], lwd = 2,
       main = paste0(codeLetters[j], " - Lag = ", lags[j], " weeks"),
       ylim = c(0, 2))
  abline(h = 1, col = "grey40", lwd = 1, lty = 2)
  abline(v = 15, col = "grey40", lwd = 1, lty = 2)
}
par(mfrow = c(1, 1))
dev.off()

png("outputs/Fig2F_temp_contour.png", width = width_set, height = height_set, res = res_set)
# par(mar = c(5, 5, 4, 2))
plot(cptmean, "contour",
     xlab = "Temperature (°C)", ylab = "Lag (weeks)",
     xlim = c(6, 30),                              # explicit xlim so abline lands correctly
     main = "F — Lag-response surface")
abline(v = 15 - 4.2, col = "black", lwd = 2, lty = 2) # corrected by and
dev.off()

# Precipitation effects

png("outputs/Fig2G_precip_overall.png", width = width_set, height = height_set, res = res_set)
# par(mar = c(5, 5, 4, 2))
plot(cppmean, ptype = "overall",
     xlab = "Precipitation (mm)", ylab = "IRR",
     col = "steelblue", lwd = 2,
     main = "G — Overall cumulative effect",
     ylim = c(0.5, 3))
abline(h = 1, col = "grey40", lwd = 1, lty = 2)   # dashed at IRR = 1
dev.off()

codeLetters = c("H", "I", "J", "K")

png("outputs/Fig2HIJK_precip_exposure_lags.png", width = width_set, height = height_set, res = res_set)
# par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
par(mfrow = c(2, 2))
cols <- c("steelblue", "darkgreen", "firebrick", "darkorange")
lags <- c(0, 2, 4, 8)
for (j in seq_along(lags)) {
  plot(cppmean, lag = lags[j],
       xlab = "Precipitation (mm)", ylab = "IRR",
       col = cols[j], lwd = 2,
       main = paste0(codeLetters[j], " - Lag = ", lags[j], " weeks"),
       ylim = c(0.5, 1.5))
  abline(h = 1, col = "grey40", lwd = 1, lty = 2)
}
par(mfrow = c(1, 1))
dev.off()

png("outputs/Fig2L_precip_contour.png", width = width_set, height = height_set, res = res_set)
# par(mar = c(5, 5, 4, 2))
plot(cppmean, "contour",
     xlab = "Precipitation (mm)", ylab = "Lag (weeks)",
     xlim = c(0, ceiling(P_max)),                  # explicit xlim so abline lands correctly
     main = "L — Lag-response surface")
abline(v = 0-5.2, col = "black", lwd = 2, lty = 2) # aggiustato a mano
dev.off()
