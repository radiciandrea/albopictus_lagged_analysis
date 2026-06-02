while (!is.null(dev.list())) dev.off()

width_set = 1600
height_set = width_set*8/9
res_set = 300

# Temperature effects

png("outputs/Fig2A_temp_overall.png", width = width_set, height = height_set, res = res_set)
# par(mar = c(5, 5, 4, 2))
plot(cptmean, ptype = "overall",
     xlab = "Average temperature (°C)", ylab = "IRR",
     col = "firebrick", lwd = 2,
     main = "A — Overall cumulative effect of temperature")
abline(h = 1, col = "grey40", lwd = 1, lty = 2)
dev.off()

png("outputs/Fig2B_temp_exposure_lags.png", width = width_set, height = height_set, res = res_set)
par(mfrow = c(2, 2), mar = c(4.5, 4, 2.5, 1))

# layout(mat = matrix(c(1, 2, 3, 4), ncol = 2, byrow = T))
# par(mfrow = c(2, 2))
cols <- c("steelblue", "darkgreen", "firebrick", "darkorange")
lags <- c(0, 2, 4, 8)
for (j in seq_along(lags)) {
  plot(cptmean, lag = lags[j],
       xlab = "Average temperature (°C)", ylab = "IRR",
       col = cols[j], lwd = 2,
       ylim = c(0, 2))
  text(x = 6, y = 1.8, label = paste0("Lag = ", lags[j], "\nweeks"), cex = 1, adj = 0)
  
  abline(h = 1, col = "grey40", lwd = 1, lty = 2)
  abline(v = 15, col = "grey40", lwd = 1, lty = 2)
}
title("B - Effect of average temperature at different lags", outer = TRUE, line = -1.5)
# par(mfrow = c(1, 1))
dev.off()

png("outputs/Fig2C_temp_contour.png", width = width_set, height = height_set, res = res_set)
# par(mar = c(5, 5, 4, 2))
plot(cptmean, "contour",
     xlab = "Average temperature (°C)", ylab = "Lag (weeks)",
     xlim = c(6, 30),                              # explicit xlim so abline lands correctly
     main = "C — Lag-response surface")
abline(v = 15 - 4.2, col = "black", lwd = 2, lty = 2) # corrected by and
dev.off()

# Precipitation effects

png("outputs/Fig2D_precip_overall.png", width = width_set, height = height_set, res = res_set)
# par(mar = c(5, 5, 4, 2))
plot(cppmean, ptype = "overall",
     xlab = "Cumulative precipitation (mm)", ylab = "IRR",
     col = "steelblue", lwd = 2,
     main = "D — Overall cumulative effect of precipitation",
     ylim = c(0.5, 3))
abline(h = 1, col = "grey40", lwd = 1, lty = 2)   # dashed at IRR = 1
dev.off()

png("outputs/Fig2E_precip_exposure_lags.png", width = width_set, height = height_set, res = res_set)
par(mfrow = c(2, 2), mar = c(4.5, 4, 2.5, 1))
# par(mfrow = c(2, 2))
cols <- c("steelblue", "darkgreen", "firebrick", "darkorange")
lags <- c(0, 2, 4, 8)
for (j in seq_along(lags)) {
  plot(cppmean, lag = lags[j],
       xlab = "Cumulative precipitation (mm)", ylab = "IRR",
       col = cols[j], lwd = 2,
       ylim = c(0.5, 1.5))
  text(x = 2, y = 1.4, label = paste0("Lag = ", lags[j], " weeks"), cex = 1, adj = 0)
  
  abline(h = 1, col = "grey40", lwd = 1, lty = 2)
}
title("E - Effect of cumulative precipitation at different lags", outer = TRUE, line = -1.5)
# par(mfrow = c(1, 1))
dev.off()

png("outputs/Fig2F_precip_contour.png", width = width_set, height = height_set, res = res_set)
# par(mar = c(5, 5, 4, 2))
plot(cppmean, "contour",
     xlab = "Cumulative precipitation (mm)", ylab = "Lag (weeks)",
     xlim = c(0, ceiling(P_max)),                  # explicit xlim so abline lands correctly
     main = "F — Lag-response surface")
abline(v = 0-5.2, col = "black", lwd = 2, lty = 2) # aggiustato a mano
dev.off()
