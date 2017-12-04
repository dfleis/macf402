##
## Script for exploring for a suitable smoothing parameter h
##
setwd("~/drive/concordia/2015-2016/1_fall2015/macf402/assignment3/")

#==== DEFINE FXNS ====#
K <- function(h, x) { # kernel smoother
  1/(h * sqrt(2 * pi)) * exp(-1/2 * (x/h)^2)
}
nad_wat_1d <- function(h, x, x_vect, y_vect) { # Nadayara-Watson estimator
  sum(K(h, x - x_vect) * y_vect)/sum(K(h, x - x_vect))
}
smooth_vols <- function(spot, data, h, min, max, dx) { 
  # implmenetation of the Nadayara-Watson estimator
  unobs_moneyness <- seq(min, max, dx)
  moneyness <- spot/data$strike
  
  sapply(unobs_moneyness, FUN = function(z) {
    nad_wat_1d(h, z, moneyness, data$imp_vol)})
}

#==== LOAD DATA ====#
mar_dat <- read.csv("./data/imp_vols_mar.csv", stringsAsFactors = F)
apr_dat <- read.csv("./data/imp_vols_apr.csv", stringsAsFactors = F)
may_dat <- read.csv("./data/imp_vols_may.csv", stringsAsFactors = F)

#==== SET PARAMS ====#
h <- c(0.025,0.05,0.075,0.1) # smoothing parameters to test
spot <- 502.12 # market spot price on Feb 17 2012
min <- 0.5 # min unobs. moneyness to compute
max <- 1.5 # max unobs. moneyness to compute
dx <- 0.05 # granularity of unobs moneyness steps
unobs_moneyness <- seq(min, max, dx)

#==== MARCH DATA ====#
sigma <- mar_dat$imp_vol
obs_moneyness <- mar_dat$moneyness

# compute fitted vols
M_sigma_hat <- matrix(nrow = length(unobs_moneyness), ncol = length(h))
M_sigma_hat <- sapply(h, FUN = function(k) {
  smooth_vols(spot, mar_dat, k, min, max, dx)})
colnames(M_sigma_hat) <- h

# plot data
pdf("./plots/fitted_smile_mar.pdf", height = 6, width = 6)
par(mfrow = c(2,2), oma = c(0, 0, 2, 0)) 
plot(sigma~obs_moneyness,
     xlab = "Moneyness = S/K", ylab = "Implied Volatility",
     main = paste0("h = ", h[1]))
lines(M_sigma_hat[,1]~unobs_moneyness, lwd = 2)
for (i in 2:length(h)) {
  plot(sigma~obs_moneyness,
       xlab = "Moneyness = S/K", ylab = "Implied Volatility",
       main = paste0("h = ", h[i]))
  lines(M_sigma_hat[,i]~unobs_moneyness, lwd = 2)
}
mtext(paste0("Fitted Smile: AAPL Options with Maturity ", mar_dat$date[1]), outer = T, cex = 1.25)
dev.off()

#==== APRIL DATA ====#
sigma <- apr_dat$imp_vol
obs_moneyness <- apr_dat$moneyness

# compute fitted vols
M_sigma_hat <- matrix(nrow = length(unobs_moneyness), ncol = length(h))
M_sigma_hat <- sapply(h, FUN = function(k) {
  smooth_vols(spot, apr_dat, k, min, max, dx)})
colnames(M_sigma_hat) <- h

# plot data
pdf("./plots/fitted_smile_apr.pdf", height = 6, width = 6)
par(mfrow = c(2,2), oma = c(0, 0, 2, 0)) 
plot(sigma~obs_moneyness,
     xlab = "Moneyness = S/K", ylab = "Implied Volatility",
     main = paste0("h = ", h[1]))
lines(M_sigma_hat[,1]~unobs_moneyness, lwd = 2)
for (i in 2:length(h)) {
  plot(sigma~obs_moneyness,
       xlab = "Moneyness = S/K", ylab = "Implied Volatility",
       main = paste0("h = ", h[i]))
  lines(M_sigma_hat[,i]~unobs_moneyness, lwd = 2)
}
mtext(paste0("Fitted Smile: AAPL Options with Maturity ", apr_dat$date[1]), outer = T, cex = 1.25)
dev.off()

#==== MAY DATA ====#
sigma <- may_dat$imp_vol
obs_moneyness <- may_dat$moneyness

# compute fitted vols
M_sigma_hat <- matrix(nrow = length(unobs_moneyness), ncol = length(h))
M_sigma_hat <- sapply(h, FUN = function(k) {
  smooth_vols(spot, may_dat, k, min, max, dx)})
colnames(M_sigma_hat) <- h

# plot data
pdf("./plots/fitted_smile_may.pdf", height = 6, width = 6)
par(mfrow = c(2,2), oma = c(0, 0, 2, 0)) 
plot(sigma~obs_moneyness,
     xlab = "Moneyness = S/K", ylab = "Implied Volatility",
     main = paste0("h = ", h[1]))
lines(M_sigma_hat[,1]~unobs_moneyness, lwd = 2)
for (i in 2:length(h)) {
  plot(sigma~obs_moneyness,
       xlab = "Moneyness = S/K", ylab = "Implied Volatility",
       main = paste0("h = ", h[i]))
  lines(M_sigma_hat[,i]~unobs_moneyness, lwd = 2)
}
mtext(paste0("Fitted Smile: AAPL Options with Maturity ", may_dat$date[1]), outer = T, cex = 1.25)
dev.off()















