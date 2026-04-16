while (!is.null(dev.list())) dev.off()

# Temperature effects

png("outputs/Fig2A_temp_overall.png", width = 1800, height = 1600, res = 200)
par(mar = c(5, 5, 4, 2))
plot(cptmean, ptype = "overall",
     xlab = "Temperature (°C)", ylab = "IRR",
     col = "firebrick", lwd = 2,
     main = "A — Overall cumulative effect")
abline(h = 1, col = "grey40", lwd = 1, lty = 2)
dev.off()

png("outputs/Fig2B_temp_exposure_lags.png", width = 1800, height = 1600, res = 200)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
cols <- c("steelblue", "darkgreen", "firebrick", "darkorange")
lags <- c(0, 2, 4, 8)
for (j in seq_along(lags)) {
  plot(cptmean, lag = lags[j],
       xlab = "Temperature (°C)", ylab = "IRR",
       col = cols[j], lwd = 2,
       main = paste0("Lag = ", lags[j], " weeks"),
       ylim = c(0, 3))
  abline(h = 1, col = "grey40", lwd = 1, lty = 2)
  abline(v = 15, col = "grey40", lwd = 1, lty = 2)
}
par(mfrow = c(1, 1))
dev.off()

png("outputs/Fig2C_temp_contour.png", width = 1800, height = 1600, res = 200)
par(mar = c(5, 5, 4, 2))
plot(cptmean, "contour",
     xlab = "Temperature (°C)", ylab = "Lag (weeks)",
     xlim = c(6, 30),                              # explicit xlim so abline lands correctly
     main = "C — Lag-response surface")
abline(v = 15, col = "black", lwd = 2, lty = 2)
dev.off()

# Precipitation effects

png("outputs/Fig3A_precip_overall.png", width = 1800, height = 1600, res = 200)
par(mar = c(5, 5, 4, 2))
plot(cppmean, ptype = "overall",
     xlab = "Precipitation (mm)", ylab = "IRR",
     col = "steelblue", lwd = 2,
     main = "A — Overall cumulative effect",
     ylim = c(0.5, 3))
abline(h = 1, col = "grey40", lwd = 1, lty = 2)   # dashed at IRR = 1
dev.off()

png("outputs/Fig3B_precip_exposure_lags.png", width = 1800, height = 1600, res = 200)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
cols <- c("steelblue", "darkgreen", "firebrick", "darkorange")
lags <- c(0, 2, 4, 8)
for (j in seq_along(lags)) {
  plot(cppmean, lag = lags[j],
       xlab = "Precipitation (mm)", ylab = "IRR",
       col = cols[j], lwd = 2,
       main = paste0("Lag = ", lags[j], " weeks"),
       ylim = c(0, 2))
  abline(h = 1, col = "grey40", lwd = 1, lty = 2)
}
par(mfrow = c(1, 1))
dev.off()

png("outputs/Fig3C_precip_contour.png", width = 1800, height = 1600, res = 200)
par(mar = c(5, 5, 4, 2))
plot(cppmean, "contour",
     xlab = "Precipitation (mm)", ylab = "Lag (weeks)",
     xlim = c(0, ceiling(P_max)),                  # explicit xlim so abline lands correctly
     main = "C — Lag-response surface")
abline(v = 0, col = "black", lwd = 2, lty = 2)
dev.off()