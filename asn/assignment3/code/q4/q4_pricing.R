##
## Script for pricing 
##
setwd("~/drive/concordia/2015-2016/1_fall2015/macf402/assignment3/")

fit_2d_data <- read.csv("./data/output/q2_fitted_2d_data.csv", stringsAsFactors = F)

spot <- 502.12
r <- 0.0175
K <- 572.50

maturity_wanted <- 60/365
moneyness_wanted <- spot/K

moneyness_err <- abs(fit_2d_data$moneyness - moneyness_wanted)
maturity_err <- abs(fit_2d_data$maturity - maturity_wanted)
idx <- which( (moneyness_err == min(moneyness_err)) & (maturity_err == min(maturity_err)))
sigma_hat <- fit_2d_data[idx,]$fitted_vol

bs_price <- function(S, tau, K, r, sigma, type) {
  d1 <- 1/(sigma * sqrt(tau)) * (log(S/K) + (r + 1/2 * sigma^2) * tau)
  d2 <- d1 - sigma * sqrt(tau)
  
  if (type == 'c') {
    pnorm(d1) * S - pnorm(d2) * K * exp(-r * tau) 
  } else {
    pnorm(-d2) * K * exp(-r * tau) - pnorm(-d1) * S
  }
}
bs_delta <- function(S, tau, K, r, sigma, type) {
  d1 <- 1/(sigma * sqrt(tau)) * (log(S/K) + (r + 1/2 * sigma^2) * tau)
  
  if (type == 'c') {
    pnorm(d1)
  } else {
    -pnorm(-d1)
  }
}


bs_price(spot, maturity_wanted, K, r, sigma_hat, 'c')
bs_delta(spot, maturity_wanted, K, r, sigma_hat, 'c')

















