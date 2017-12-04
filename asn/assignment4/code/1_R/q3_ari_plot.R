## Set working directory
setwd("~/drive/concordia/2015-2016/1_fall2015/macf402/assignment4/")

## Load data
dat_crude <- read.csv("./data/q3_crude_est.csv")
dat_anti <- read.csv("./data/q3_anti_est.csv")
dat_cont <- read.csv("./data/q3_cont_est.csv")

## Get overall extrema for the y coordinates
y_min <- min(dat_crude$CI_lo, dat_anti$CI_lo, dat_cont$CI_lo)
y_max <- max(dat_crude$CI_hi, dat_anti$CI_hi, dat_cont$CI_hi)

## Plot means w/ upper and lower confidence intervals as a function of N sims
pdf("./plots/q3/ari_asian_call_sim_est.pdf", height = 7.5, width = 6)
par(mfrow = c(3,2),  oma = c(0, 0, 2, 0))
plot(dat_crude$mean ~ dat_crude$N, log = 'x', type = 'b', ylim = c(y_min, y_max),
     main = 'Crude Estimator: 95% CI', xlab = 'N', ylab = 'Asian Call Value' )
segments(x0 = dat_crude$N, y0 = dat_crude$CI_lo, x1 = dat_crude$N, y1 = dat_crude$CI_hi, lwd = 2)
axis(1, at = dat_crude$N)
plot(dat_crude$mean ~ dat_crude$N, log = 'x', type = 'b', ylim = c(min(dat_crude$CI_lo), max(dat_crude$CI_hi)),
     main = 'Crude Estimator: 95% CI', xlab = 'N', ylab = 'Asian Call Value' )
segments(x0 = dat_crude$N, y0 = dat_crude$CI_lo, x1 = dat_crude$N, y1 = dat_crude$CI_hi, lwd = 2)
axis(1, at = dat_crude$N)

plot(dat_anti$mean ~ dat_anti$N, log = 'x', type = 'b', ylim = c(y_min, y_max),
     main = 'Antithetic Estimator: 95% CI', xlab = 'N', ylab = 'Asian Call Value')
segments(x0 = dat_anti$N, y0 = dat_anti$CI_lo, x1 = dat_anti$N, y1 = dat_anti$CI_hi, lwd = 2)
axis(1, at = dat_anti$N)
plot(dat_anti$mean ~ dat_anti$N, log = 'xy', type = 'b', ylim = c(min(dat_anti$CI_lo), max(dat_anti$CI_hi)),
     main = 'Antithetic Estimator: 95% CI', xlab = 'N', ylab = 'Asian Call Value')
segments(x0 = dat_anti$N, y0 = dat_anti$CI_lo, x1 = dat_anti$N, y1 = dat_anti$CI_hi, lwd = 2)
axis(1, at = dat_anti$N)

plot(dat_cont$mean ~ dat_cont$N, log = 'x', type = 'b', ylim = c(y_min, y_max), 
     main = 'Control Variate Estimator: 95% CI', xlab = 'N', ylab = 'Asian Call Value',)
segments(x0 = dat_cont$N, y0 = dat_cont$CI_lo, x1 = dat_cont$N, y1 = dat_cont$CI_hi, lwd = 2)
axis(1, at = dat_cont$N)
plot(dat_cont$mean ~ dat_cont$N, log = 'x', type = 'b', ylim = c(min(dat_cont$CI_lo), max(dat_cont$CI_hi)),
     main = 'Control Variate Estimator: 95% CI', xlab = 'N', ylab = 'Asian Call Value',)
segments(x0 = dat_cont$N, y0 = dat_cont$CI_lo, x1 = dat_cont$N, y1 = dat_cont$CI_hi, lwd = 2)
axis(1, at = dat_cont$N)
mtext("Arithmetic Asian Call Option Price", outer = T, cex = 1.25)
dev.off()

## Compute confidence interval width/difference of max - min
dat_crude$diff <- dat_crude$CI_hi - dat_crude$CI_lo
dat_anti$diff <- dat_anti$CI_hi - dat_anti$CI_lo
dat_cont$diff <- dat_cont$CI_hi - dat_cont$CI_lo

## Find extrema for plotting
ymin_diff <- min(dat_crude$diff, dat_anti$diff, dat_cont$diff)
ymax_diff <- max(dat_crude$diff, dat_anti$diff, dat_cont$diff)

## Plot confidence interval widths as a function of N sims
pdf("./plots/q3/ari_asian_call_CI_widths.pdf", height = 4, width = 4)
par(mfrow = c(1,1))
plot(dat_crude$diff ~ dat_crude$N, log = 'xy', type = 'b', ylim = c(ymin_diff, ymax_diff), 
     pch = 1, lwd = 2, main = '95% Confidence Interval Width', xlab = 'N', ylab = 'Width')
lines(dat_anti$diff ~ dat_anti$N, type = 'b', pch = 2, lwd = 2)
lines(dat_cont$diff ~ dat_cont$N, type = 'b', pch = 0, lwd = 2)
legend('topright', legend = c('Crude', 'Antithetic', 'Control'), pch = c(1,2,0))
axis(1, at = dat_cont$N)
dev.off()

## Find time extrema for plotting
ymin_time <- min(dat_crude$time, dat_anti$time, dat_cont$time)
ymax_time <- max(dat_crude$time, dat_anti$time, dat_cont$time)

## Plot confidence interval widths as a function of N sims
pdf("./plots/q3/ari_asian_call_time.pdf", height = 4, width = 4)
par(mfrow = c(1,1))
plot(dat_crude$time ~ dat_crude$N, log = 'xy', type = 'b', ylim = c(ymin_time, ymax_time), 
     pch = 1, lwd = 2, main = 'Computation Time', xlab = 'N', ylab = 'Time')
lines(dat_anti$time ~ dat_anti$N, type = 'b', pch = 2, lwd = 2)
lines(dat_cont$time ~ dat_cont$N, type = 'b', pch = 0, lwd = 2)
legend('topleft', legend = c('Crude', 'Antithetic', 'Control'), pch = c(1,2,0))
axis(1, at = dat_crude$N)
dev.off()

## Compute efficiency w/o computation time
dat_crude$eff_wo_time <- 1/(dat_crude$var)
dat_anti$eff_wo_time <- 1/(dat_anti$var)
dat_cont$eff_wo_time <- 1/(dat_cont$var)

ymin_eff_wo <- min(dat_crude$eff_wo_time, dat_anti$eff_wo_time, dat_cont$eff_wo_time)
ymax_eff_wo <- max(dat_crude$eff_wo_time, dat_anti$eff_wo_time, dat_cont$eff_wo_time)

pdf("./plots/q3/ari_asian_call_eff_wo_time.pdf", height = 4, width = 4)
par(mfrow = c(1,1))
plot(dat_crude$eff_wo_time ~ dat_crude$N, log = 'xy', type = 'b', lwd = 2, ylim = c(ymin_eff_wo, 100*10^ceiling(log10(ymax_eff_wo))),
     main = 'Efficiency without\nComputation Time', xlab = 'N', ylab = 'Efficiency')
lines(dat_anti$eff_wo_time ~ dat_anti$N, type = 'b', pch = 2, lwd = 2)
lines(dat_cont$eff_wo_time ~ dat_cont$N, type = 'b', pch = 2, lwd = 2)
legend('topright', legend = c('Crude', 'Antithetic', 'Control'), pch = c(1,2,0))
axis(1, at = dat_anti$N)
dev.off()










