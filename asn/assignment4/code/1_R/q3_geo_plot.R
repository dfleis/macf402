## Set working directory
setwd("~/drive/concordia/2015-2016/1_fall2015/macf402/assignment4/")

## Load data
dat_crude <- read.csv("./data/q3_geo_crude_est.csv")
dat_anti <- read.csv("./data/q3_geo_anti_est.csv")

## true geometric asian call price
true_call <- dat_crude$geo_price[1]

## Get overall extrema for the y coordinates
y_min <- min(dat_crude$CI_lo, dat_anti$CI_lo)
y_max <- max(dat_crude$CI_hi, dat_anti$CI_hi)
err_min <- min( abs(c(dat_crude$mean, dat_anti$mean) - true_call) )
err_max <- max( abs(c(dat_crude$mean, dat_anti$mean) - true_call) )

## Plot means w/ upper and lower confidence intervals as a function of N sims
pdf("./plots/q3/geo_asian_call_sim_est.pdf", height = 6, width = 6)
par(mfrow = c(2,2),  oma = c(0, 0, 2, 0))
plot(dat_crude$mean ~ dat_crude$N, log = 'x', type = 'b', ylim = c(y_min, y_max),
     main = 'Crude Estimator: 95% CI', xlab = 'N', ylab = 'Asian Call Value')
abline(h = dat_crude$geo_price, lwd = 1, col = 'gray70')
lines(dat_crude$mean ~ dat_crude$N, log = 'x', type = 'b')
segments(x0 = dat_crude$N, y0 = dat_crude$CI_lo, x1 = dat_crude$N, y1 = dat_crude$CI_hi, lwd = 2)
axis(1, at = dat_crude$N)

plot(abs(dat_crude$mean - true_call) ~ dat_crude$N, log = 'xy', type = 'b', ylim = c(err_min, err_max),
     main = 'Crude Estimator\nError to Analytic Price', xlab = 'N', ylab = 'Price Error' )
axis(1, at = dat_crude$N)

plot(dat_anti$mean ~ dat_anti$N, log = 'x', type = 'b', ylim = c(y_min, y_max),
     main = 'Antithetic Estimator: 95% CI', xlab = 'N', ylab = 'Asian Call Value')
abline(h = dat_crude$geo_price, lwd = 1, col = "gray70")
lines(dat_anti$mean ~ dat_anti$N, log = 'x', type = 'b')
segments(x0 = dat_anti$N, y0 = dat_anti$CI_lo, x1 = dat_anti$N, y1 = dat_anti$CI_hi, lwd = 2)
axis(1, at = dat_anti$N)

plot(abs(dat_anti$mean - true_call) ~ dat_anti$N, log = 'xy', type = 'b', ylim = c(err_min, err_max),
     main = 'Antithetic Estimator:\nError to Analytic Price', xlab = 'N', ylab = 'Price Error' )
axis(1, at = dat_anti$N)
mtext("Geometric Asian Call Option Price", outer = T, cex = 1.25)
dev.off()

## Compute confidence interval width/difference of max - min
dat_crude$diff <- dat_crude$CI_hi - dat_crude$CI_lo
dat_anti$diff <- dat_anti$CI_hi - dat_anti$CI_lo

## Find extrema for plotting
ymin_diff <- min(dat_crude$diff, dat_anti$diff)
ymax_diff <- max(dat_crude$diff, dat_anti$diff)

## Plot confidence interval widths as a function of N sims
pdf("./plots/q3/geo_asian_call_CI_widths.pdf", height = 4, width = 4)
par(mfrow = c(1,1))
plot(dat_crude$diff ~ dat_crude$N, log = 'xy', type = 'b', ylim = c(ymin_diff, ymax_diff), 
     pch = 1, lwd = 2, main = '95% Confidence Interval Width', xlab = 'N', ylab = 'Width')
lines(dat_anti$diff ~ dat_anti$N, type = 'b', pch = 2, lwd = 2)
legend('topright', legend = c('Crude', 'Antithetic'), pch = c(1,2))
axis(1, at = dat_crude$N)
dev.off()

## Find time extrema for plotting
ymin_time <- min(dat_crude$time, dat_anti$time)
ymax_time <- max(dat_crude$time, dat_anti$time)

## Plot confidence interval widths as a function of N sims
pdf("./plots/q3/geo_asian_call_time.pdf", height = 4, width = 4)
par(mfrow = c(1,1))
plot(dat_crude$time ~ dat_crude$N, log = 'xy', type = 'b', ylim = c(ymin_time, ymax_time), 
     pch = 1, lwd = 2, main = 'Computation Time', xlab = 'N', ylab = 'Time')
lines(dat_anti$time ~ dat_anti$N, type = 'b', pch = 2, lwd = 2)
legend('topleft', legend = c('Crude', 'Antithetic'), pch = c(1,2))
axis(1, at = dat_crude$N)
dev.off()

## Compute efficiency w/o computation time
dat_crude$eff_wo_time <- 1/(dat_crude$var)
dat_anti$eff_wo_time <- 1/(dat_anti$var)

ymin_eff_wo <- min(dat_crude$eff_wo_time, dat_anti$eff_wo_time)
ymax_eff_wo <- max(dat_crude$eff_wo_time, dat_anti$eff_wo_time)

pdf("./plots/q3/geo_asian_call_eff_wo_time.pdf", height = 4, width = 4)
par(mfrow = c(1,1))
plot(dat_crude$eff_wo_time ~ dat_crude$N, log = 'xy', type = 'b', lwd = 2, ylim = c(ymin_eff_wo, 10^ceiling(log10(ymax_eff_wo))),
     main = 'Efficiency without\nComputation Time', xlab = 'N', ylab = 'Efficiency')
lines(dat_anti$eff_wo_time ~ dat_anti$N, type = 'b', pch = 2, lwd = 2)
legend('topright', legend = c('Crude', 'Antithetic'), pch = c(1,2))
axis(1, at = dat_anti$N)
dev.off()










