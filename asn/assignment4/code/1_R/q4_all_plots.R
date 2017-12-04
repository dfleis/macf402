## Set working directory
setwd("~/drive/concordia/2015-2016/1_fall2015/macf402/assignment4/")

## Load data
dat_crude <- read.csv("./data/q4a_clean_heston_crude_est.csv")
dat_cmc <- read.csv("./data/q4b_clean_heston_cmc_est.csv")

## add ID variable
dat_crude$est <- "Crude"
dat_cmc$est <- "CMC"

## combine data
all_dat <- rbind(dat_crude, dat_cmc)

## Compute confidence interval width/difference of max - min
all_dat$diff <- all_dat$CI.Hi - all_dat$CI.Lo
## Find extrema for plotting
ymin_price <- min(all_dat$CI.Low)
ymax_price <- max(all_dat$CI.High)
ymin_diff <- min(all_dat$diff)
ymax_diff <- max(all_dat$diff)
ymin_time <- min(all_dat$Time)
ymax_time <- max(all_dat$Time)
ymin_eff <- min(all_dat$Eff.)
ymax_eff <- max(all_dat$Eff.)

## Plot means w/ upper and lower confidence intervals as a function of N sims
pdf("./plots/q4/heston_call_est.pdf", height = 7.5, width = 6)
par(mfrow = c(3,2),  oma = c(0, 0, 2, 0))
for (j in 1:length(unique(all_dat$rho))) {
  for (i in 1:length(unique(all_dat$est))) {
    sub_dat <- subset(all_dat, rho == unique(all_dat$rho)[j] & est == unique(all_dat$est)[i])
    plot(sub_dat$Mean ~ sub_dat$N, log = 'x', type = 'b', ylim = c(ymin_price, ymax_price), lwd = 1.5,  
         main = substitute(paste(est, " Estimator: ", rho, " = ", rh), 
                                 list(est = unique(all_dat$est)[i] ,rh = unique(all_dat$rho)[j])), 
         xlab = 'N', ylab = 'Value')
    segments(x0 = sub_dat$N, y0 = sub_dat$CI.Low, x1 = sub_dat$N, y1 = sub_dat$CI.High, lwd = 2)
  }
}
mtext("Heston Model: European Call Option Value", outer = T, cex = 1.25)
dev.off()


## Plot means w/ upper and lower confidence intervals as a function of N sims
pdf("./plots/q4/heston_call_CI_width.pdf", height = 3, width = 7)
par(mfrow = c(1,3),  oma = c(0, 0, 2, 0))
for (j in 1:length(unique(all_dat$rho))) {
  sub_dat <- subset(all_dat, rho == unique(all_dat$rho)[j])
  plot(sub_dat$diff[sub_dat$est == "Crude"] ~ unique(sub_dat$N), log = 'xy', type = 'b', ylim = c(ymin_diff, ymax_diff), 
       lwd = 1.5,  pch = 1,
       main = substitute(paste(rho, " = ", rh), list(rh = unique(all_dat$rho)[j])), xlab = 'N', ylab = 'CI Width')
  lines(sub_dat$diff[sub_dat$est == "CMC"] ~ unique(sub_dat$N), type = 'b', lwd = 1.5, pch = 2)
  legend('bottomleft', legend = c('Crude', 'CMC'), pch = c(1,2))
}
mtext("95% Confidence Interval Widths", outer = T, cex = 1.25)
dev.off()

## Plot time as a function of N sims
pdf("./plots/q4/heston_call_time.pdf", height = 3, width = 7)
par(mfrow = c(1,3),  oma = c(0, 0, 2, 0))
for (j in 1:length(unique(all_dat$rho))) {
  sub_dat <- subset(all_dat, rho == unique(all_dat$rho)[j])
  plot(sub_dat$Time[sub_dat$est == "Crude"] ~ unique(sub_dat$N), log = 'xy', type = 'b', ylim = c(ymin_time, ymax_time), 
       lwd = 1.5,  pch = 1,
       main = substitute(paste(rho, " = ", rh), list(rh = unique(all_dat$rho)[j])), xlab = 'N', ylab = 'Time')
  lines(sub_dat$Time[sub_dat$est == "CMC"] ~ unique(sub_dat$N), type = 'b', lwd = 1.5, pch = 2)
  legend('bottomright', legend = c('Crude', 'CMC'), pch = c(1,2))
}
mtext("Computation Time", outer = T, cex = 1.25)
dev.off()

## Plot mtime as a function of N sims
pdf("./plots/q4/heston_call_eff.pdf", height = 3, width = 7)
par(mfrow = c(1,3),  oma = c(0, 0, 2, 0))
for (j in 1:length(unique(all_dat$rho))) {
  sub_dat <- subset(all_dat, rho == unique(all_dat$rho)[j])
  plot(sub_dat$Eff.[sub_dat$est == "Crude"] ~ unique(sub_dat$N), log = 'xy', type = 'b', ylim = c(ymin_eff, ymax_eff), 
       lwd = 1.5,  pch = 1,
       main = substitute(paste(rho, " = ", rh), list(rh = unique(all_dat$rho)[j])), xlab = 'N', ylab = 'Time')
  lines(sub_dat$Eff.[sub_dat$est == "CMC"] ~ unique(sub_dat$N), type = 'b', lwd = 1.5, pch = 2)
  legend(3.8*1e4, 0.75, legend = c('Crude', 'CMC'), pch = c(1,2))
}
mtext("Efficiency without Computation Time", outer = T, cex = 1.25)
dev.off()

