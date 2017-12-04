##
## Script for testing arbitrage
##
setwd("~/drive/concordia/2015-2016/1_fall2015/macf402/assignment3/")

fit_2d_data <- read.csv("./data/output/q2_fitted_2d_data.csv", stringsAsFactors = F)

spot <- 502.12
r <- 0.0175

moneyness <- unique(fit_2d_data$moneyness)
maturity <- unique(fit_2d_data$maturity)
implied_vol <- matrix(fit_2d_data$fitted_vol, ncol = length(moneyness), nrow = length(maturity), byrow = T)


bs_price <- function(S,tau,K,r,sigma, type) {
  d1 <- 1/(sigma * sqrt(tau)) * (log(S/K) + (r + 1/2 * sigma^2) * tau)
  d2 <- d1 - sigma * sqrt(tau)
  
  if (type == 'c') {
    pnorm(d1) * S - pnorm(d2) * K * exp(-r * tau) 
  } else {
    pnorm(-d2) * K * exp(-r * tau) - pnorm(-d1) * S
  }
}

bs_strike_delta <- function(S, tau, K, r, sigma, type) {
  d2 <- 1/(sigma * sqrt(tau)) * (log(S/K) + (r - 1/2 * sigma^2) * tau)
  
  if (type == 'c') {
    -pnorm(d2) * exp(-r * tau)
  } else {
    (1 - pnorm(d2)) * exp(-r * tau)
  }
}

bs_strike_gamma <-  function(S, tau, K, r, sigma) {
  d2 <- 1/(sigma * sqrt(tau)) * (log(S/K) + (r - 1/2 * sigma^2) * tau)
  
  dnorm(d2) * exp(-r * tau) / (K * sigma * sqrt(tau))
}

bs_T <-  function(S, tau, K, r, sigma) {
  d2 <- 1/(sigma * sqrt(tau)) * (log(S/K) + (r - 1/2 * sigma^2) * tau)
  
  if (type == 'c') {
    -pnorm(d2) * exp(-r * tau)
  } else {
    (1 - pnorm(d2)) * exp(-r * tau)
  }
}



call_price <- matrix(nrow = nrow(implied_vol), ncol = ncol(implied_vol))
call_sd <- matrix(nrow = nrow(implied_vol), ncol = ncol(implied_vol))
put_price <- matrix(nrow = nrow(implied_vol), ncol = ncol(implied_vol))
put_sd <- matrix(nrow = nrow(implied_vol), ncol = ncol(implied_vol))
opt_sg <- matrix(nrow = nrow(implied_vol), ncol = ncol(implied_vol))

for (i in 1:length(moneyness)) {
  for (j in 1:length(maturity)) {
    call_price[j,i] <- bs_price(spot, maturity[j], spot/moneyness[i], r, implied_vol[j,i], 'c') 
    call_sd[j,i] <- bs_strike_delta(spot, maturity[j], spot/moneyness[i], r, implied_vol[j,i], 'c') 
    put_price[j,i] <- bs_price(spot, maturity[j], spot/moneyness[i], r, implied_vol[j,i], 'p')
    put_sd[j,i] <- bs_strike_delta(spot, maturity[j], spot/moneyness[i], r, implied_vol[j,i], 'p') 
    opt_sg[j,i] <- bs_strike_gamma(spot, maturity[j], spot/moneyness[i], r, implied_vol[j,i]) 
  }
}


par(mfrow = c(1,2),  oma = c(0, 0, 2, 0))
p_call <- persp(moneyness, maturity, call_price, phi = 15, theta = -30,
           col = "gray80", border = "gray30", expand = 0.5, ticktype = "detailed",
           xlab = "Moneyness = S/K", ylab = "Maturity", zlab = "Implied Volatility",
           cex.lab = 0.75, cex.axis = 0.75,
           main = "Call Options")

p_put <- persp(moneyness, maturity, put_price, phi = 15, theta = 150,
           col = "gray80", border = "gray30", expand = 0.5, ticktype = "detailed",
           xlab = "Moneyness = S/K", ylab = "Maturity", zlab = "Implied Volatility",
           cex.lab = 0.75, cex.axis = 0.75,
           main = "Put Options")
mtext("Price Surfaces given Fitted Volatility with h1 = 0.05, h2 = 0.05", outer = T, cex = 1.25)

pdf("./plots/q3/call_strike_delta.pdf", height = 6, width = 6)
par(mfrow = c(1,1))
p_call_sd <- persp(moneyness, maturity, call_sd, phi = 15, theta = -130,
                col = "gray80", border = "gray30", expand = 0.5, ticktype = "detailed",
                xlab = "Moneyness = S/K", ylab = "Maturity", zlab = "Derivative wrt Strike",
                cex.lab = 0.75, cex.axis = 0.75, zlim = c(min(call_sd), 0),
                main = "Derivative of Call Price wrt Strike")
dev.off()
pdf("./plots/q3/put_strike_delta.pdf", height = 6, width = 6)
par(mfrow = c(1,1))
p_put_sd <- persp(moneyness, maturity, put_sd, phi = 15, theta = -130,
                   col = "gray80", border = "gray30", expand = 0.5, ticktype = "detailed",
                   xlab = "Moneyness = S/K", ylab = "Maturity", zlab = "Derivative wrt Strike",
                   cex.lab = 0.75, cex.axis = 0.75, zlim = c(0, max(put_sd)),
                   main = "Derivative of Put Price wrt Strike")
dev.off()

pdf("./plots/q3/strike_gamma.pdf", height = 6, width = 6)
par(mfrow = c(1,1),  oma = c(0, 0, 2, 0))
p_opt_sg <- persp(moneyness, maturity, opt_sg, phi = 15, theta = -130,
                   col = "gray80", border = "gray30", expand = 0.5, ticktype = "detailed",
                   xlab = "Moneyness = S/K", ylab = "Maturity", zlab = "Second Derivative wrt Strike",
                   cex.lab = 0.75, cex.axis = 0.75, zlim = c(0, max(opt_sg)),
                   main = "Second Derivative of Option Price wrt Strike")
dev.off()







