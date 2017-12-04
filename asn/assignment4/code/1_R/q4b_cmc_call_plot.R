## Set working directory
setwd("~/drive/concordia/2015-2016/1_fall2015/macf402/assignment4/")

## Load data
dat_crude_n7 <- read.csv("./data/q4_heston_cmc_est_-7.csv")
dat_crude_0 <- read.csv("./data/q4_heston_cmc_est_0.csv")
dat_crude_7 <- read.csv("./data/q4_heston_cmc_est_7.csv")

## Get overall extrema for the y coordinates
y_min <- min(dat_crude_n7$CI_lo, dat_crude_0$CI_lo, dat_crude_7$CI_lo)
y_max <- max(dat_crude_n7$CI_hi, dat_crude_0$CI_hi, dat_crude_7$CI_hi)

## Plot means w/ upper and lower confidence intervals as a function of N sims
pdf("./plots/q4/heston_cmc_call_est.pdf", height = 7.5, width = 5)
par(mfrow = c(3,1),  oma = c(0, 0, 2, 0))
plot(dat_crude_n7$mean ~ dat_crude_n7$N, log = 'x', type = 'b', ylim = c(min(dat_crude_n7$CI_lo), max(dat_crude_n7$CI_hi)),
     main = substitute(paste(rho, " = ", rh), list(rh = -0.7)), xlab = 'N', ylab = 'Value' )
segments(x0 = dat_crude_n7$N, y0 = dat_crude_n7$CI_lo, x1 = dat_crude_n7$N, y1 = dat_crude_n7$CI_hi, lwd = 2)
axis(1, at = dat_crude_n7$N)

plot(dat_crude_0$mean ~ dat_crude_0$N, log = 'x', type = 'b', ylim = c(min(dat_crude_0$CI_lo), max(dat_crude_0$CI_hi)),
     main = substitute(paste(rho, " = ", rh), list(rh = 0)), xlab = 'N', ylab = 'Value' )
segments(x0 = dat_crude_0$N, y0 = dat_crude_0$CI_lo, x1 = dat_crude_0$N, y1 = dat_crude_0$CI_hi, lwd = 2)
axis(1, at = dat_crude_0$N)

plot(dat_crude_7$mean ~ dat_crude_7$N, log = 'x', type = 'b', ylim = c(min(dat_crude_7$CI_lo), max(dat_crude_7$CI_hi)),
     main = substitute(paste(rho, " = ", rh), list(rh = 0.7)), xlab = 'N', ylab = 'Value' )
segments(x0 = dat_crude_7$N, y0 = dat_crude_7$CI_lo, x1 = dat_crude_7$N, y1 = dat_crude_7$CI_hi, lwd = 2)
axis(1, at = dat_crude_7$N)
mtext("Heston Model: European Call Option Price (CMC)", outer = T, cex = 1.25)
dev.off()
# 
# ## Compute confidence interval width/difference of max - min
# dat_crude$diff <- dat_crude$CI_hi - dat_crude$CI_lo
# dat_anti$diff <- dat_anti$CI_hi - dat_anti$CI_lo
# dat_cont$diff <- dat_cont$CI_hi - dat_cont$CI_lo
# 
# ## Find extrema for plotting
# ymin_diff <- min(dat_crude$diff, dat_anti$diff, dat_cont$diff)
# ymax_diff <- max(dat_crude$diff, dat_anti$diff, dat_cont$diff)
# 
# ## Plot confidence interval widths as a function of N sims
# pdf("./plots/q3/ari_asian_call_CI_widths.pdf", height = 4, width = 4)
# par(mfrow = c(1,1))
# plot(dat_crude$diff ~ dat_crude$N, log = 'xy', type = 'b', ylim = c(ymin_diff, ymax_diff), 
#      pch = 1, lwd = 2, main = '95% Confidence Interval Width', xlab = 'N', ylab = 'Width')
# lines(dat_anti$diff ~ dat_anti$N, type = 'b', pch = 2, lwd = 2)
# lines(dat_cont$diff ~ dat_cont$N, type = 'b', pch = 0, lwd = 2)
# legend('topright', legend = c('Crude', 'Antithetic', 'Control'), pch = c(1,2,0))
# axis(1, at = dat_cont$N)
# dev.off()
# 
# ## Find time extrema for plotting
# ymin_time <- min(dat_crude$time, dat_anti$time, dat_cont$time)
# ymax_time <- max(dat_crude$time, dat_anti$time, dat_cont$time)
# 
# ## Plot confidence interval widths as a function of N sims
# pdf("./plots/q3/ari_asian_call_time.pdf", height = 4, width = 4)
# par(mfrow = c(1,1))
# plot(dat_crude$time ~ dat_crude$N, log = 'xy', type = 'b', ylim = c(ymin_time, ymax_time), 
#      pch = 1, lwd = 2, main = 'Computation Time', xlab = 'N', ylab = 'Time')
# lines(dat_anti$time ~ dat_anti$N, type = 'b', pch = 2, lwd = 2)
# lines(dat_cont$time ~ dat_cont$N, type = 'b', pch = 0, lwd = 2)
# legend('topleft', legend = c('Crude', 'Antithetic', 'Control'), pch = c(1,2,0))
# axis(1, at = dat_crude$N)
# dev.off()
# 
# ## Compute efficiency w/o computation time
# dat_crude$eff_wo_time <- 1/(dat_crude$var)
# dat_anti$eff_wo_time <- 1/(dat_anti$var)
# dat_cont$eff_wo_time <- 1/(dat_cont$var)
# 
# ymin_eff_wo <- min(dat_crude$eff_wo_time, dat_anti$eff_wo_time, dat_cont$eff_wo_time)
# ymax_eff_wo <- max(dat_crude$eff_wo_time, dat_anti$eff_wo_time, dat_cont$eff_wo_time)
# 
# pdf("./plots/q3/ari_asian_call_eff_wo_time.pdf", height = 4, width = 4)
# par(mfrow = c(1,1))
# plot(dat_crude$eff_wo_time ~ dat_crude$N, log = 'xy', type = 'b', lwd = 2, ylim = c(ymin_eff_wo, 100*10^ceiling(log10(ymax_eff_wo))),
#      main = 'Efficiency without\nComputation Time', xlab = 'N', ylab = 'Efficiency')
# lines(dat_anti$eff_wo_time ~ dat_anti$N, type = 'b', pch = 2, lwd = 2)
# lines(dat_cont$eff_wo_time ~ dat_cont$N, type = 'b', pch = 2, lwd = 2)
# legend('topright', legend = c('Crude', 'Antithetic', 'Control'), pch = c(1,2,0))
# axis(1, at = dat_anti$N)
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
