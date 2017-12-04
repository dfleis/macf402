setwd("~/drive/concordia/2015-2016/1_fall2015/macf402/assignment5/")
dat <- read.csv("./data/q6_dat.csv")
dat_clean <- dat

z_alpha <- 1.645
dat$short_ci_lo <- dat$short - z_alpha * dat$short_se
dat$short_ci_hi <- dat$short + z_alpha * dat$short_se
dat$long_ci_lo <- dat$long - z_alpha * dat$long_se
dat$long_ci_hi <- dat$long + z_alpha * dat$long_se

ymin <- min(dat[,-1])
ymax <- max(dat[,-1])

pdf("./plots/q6_varbreaks_all.pdf", height = 5, width = 5)
par(mfrow = c(1,1), oma=c(0,0,2,0))
plot(dat$short~dat$N, log = 'x', type = 'b', ylim = c(ymin, ymax), lwd = 2, 
     xlab = "N", ylab = "VaR Break Proportion", xaxt = "n")
axis(side = 1, at = dat$N, labels = paste0("1e",log10(dat$N)))
lines(dat$long~dat$N, type = 'b', pch = 2, lwd = 2)
legend('topright', c("Short", "Long"), pch = c(1, 2), lwd = 2)
mtext("95% VaR Breaks\nShort & Long European Call Positions\n", cex = 1.25)
dev.off()

pdf("./plots/q6_varbreaks_ci.pdf", height = 4, width = 8.5)
par(mfrow = c(1,2), oma=c(0,0,2,0))
plot(dat$short~dat$N, log = 'x', type = 'b', ylim = c(ymin, ymax), lwd = 2,
     xlab = "N", ylab = "VaR Break Proportion", main = "Short Position", xaxt = "n")
axis(side = 1, at = dat$N, labels = paste0("1e",log10(dat$N)))
segments(x0 = dat$N, y0 = dat$short_ci_lo, x1 = dat$N, y1 = dat$short_ci_hi, lwd = 2)
plot(dat$long~dat$N, log = 'x', type = 'b', ylim = c(ymin, ymax), lwd = 2,
     xlab = "N", ylab = "VaR Break Rate", main = "Long Position", xaxt = "n")
axis(side = 1, at = dat$N, labels = paste0("1e",log10(dat$N)))
segments(x0 = dat$N, y0 = dat$long_ci_lo, x1 = dat$N, y1 = dat$long_ci_hi, lwd = 2)
mtext("95% VaR Breaks: 95% Confidence Intervals", outer = T, cex = 1.25)
dev.off()

pdf("./plots/q6_varbreaks_diff.pdf", height = 5, width = 5)
par(mfrow = c(1,1), oma=c(0,0,2,0))
plot(dat$short - dat$long ~ dat$N, log = 'x', type = 'b', lwd = 2,
     xlab = "N", ylab = "Difference", xaxt = "n")
axis(side = 1, at = dat$N, labels = paste0("1e",log10(dat$N)))
mtext("Short - Long Difference\n95% VaR Breaks\n", cex = 1.25)
dev.off()

names(dat_clean) <- c("N", "Short Position", "Short Std Error", "Long Position", "Long Std Error")
dat_clean$N <- paste0("1e", log10(dat_clean$N))
write.csv(x = dat_clean, file = "./data/q6_dat_clean.csv", quote = F, row.names = F)




