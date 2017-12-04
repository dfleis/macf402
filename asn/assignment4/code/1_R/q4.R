## Set working directory
setwd("~/drive/concordia/2015-2016/1_fall2015/macf402/assignment4/")

set.seed(402)
N <- 1e4
nplot <- 20
CI <- 0.95

S0 <- 100
v0 <- 0.1
Tend <- 0.5
n <- 125

r <- 0
rho_v <- c(-0.7, 0, 0.7)
kappa <- 2
theta <- 0.03
sigma <- 0.1

dt <- Tend/n
times <- seq(0, Tend, dt)

Z1 <- matrix(rnorm(N * n), nrow = n, ncol = N)
Z2 <- matrix(rnorm(N * n), nrow = n, ncol = N)

### COMPUTE VOLS
v <- matrix(nrow = (n + 1), ncol = N)
v[1,] <- v0

for (j in 1:N) {
  for (i in 2:(n + 1)) {
    v[i,j] <- v[i - 1,j] + kappa * (theta - v[i - 1,j]) * dt + sigma * sqrt(v[i - 1,j] * dt) * Z1[i - 1,j]
    if (v[i] < 0) v[i] <- 0
  }
}
v_pctls <- apply(v, 1, quantile, probs = c(1 - CI, CI))

pdf("./plots/q4/heston_sample_vols.pdf", height = 4, width = 3.5)
par(mfrow = c(1,1), oma = c(0, 0, 2, 0))
plot(v[,1] ~ times, type = 'l', ylim = c(min(v[,1:nplot]), max(v[,1:nplot])), 
     col = 'gray70', xlab = 'Time', ylab = 'Volatility', 
     main = substitute(paste(kappa, " = ", k, ", ", theta, " = ", th, ", ", sigma, " = ", s), 
                       list(k = kappa, th = theta, s = sigma)))
for (i in 2:nplot) {
  lines(v[,i] ~ times, col = 'gray70')
}
lines(rowMeans(v) ~ times, lwd = 2)
lines(v_pctls[1,] ~ times, lwd = 2)
lines(v_pctls[2,] ~ times, lwd = 2)
mtext('Heston Model: Volatility Paths', outer = T, cex = 1.25)
dev.off()

### COMPUTE PRICES
pdf("./plots/q4/heston_sample_prices.pdf", height = 3.25, width = 8.5)
par(mfrow = c(1,3), oma = c(0, 0, 2, 0))

for (k in 1:length(rho_v)) {
  rho <- rho_v[k]
  x <- matrix(nrow = (n + 1), ncol = N)
  x[1,] <- log(S0)
  
  for (j in 1:N) {
    for (i in 2:(n + 1)) {
      x[i,j] <- x[i - 1,j] + (r - 0.5 * v[i - 1,j]) * dt + 
        sqrt(v[i - 1,j]) * (rho * sqrt(dt) * Z1[i - 1,j] + sqrt(1 - rho^2) * sqrt(dt) * Z2[i - 1,j])
    }
  }
  
  prices <- exp(x)
  price_pctls <- apply(prices, 1, quantile, probs = c(1 - CI, CI))
  
  plot(prices[,1] ~ times, type = 'l',
       ylim = c(65,155),
       col = 'gray70', xlab = 'Time', ylab = 'Price',
       main = substitute(paste(r, " = ", ir, ", ", rho, " = ", rh), 
                         list(ir = r, rh = rho_v[k])))
  for (i in 2:nplot) {
    lines(prices[,i] ~ times, col = 'gray70')
  }
  lines(rowMeans(prices) ~ times, lwd = 2)
  lines(price_pctls[1,] ~ times, lwd = 2)
  lines(price_pctls[2,] ~ times, lwd = 2)
}
mtext("Heston Model: Price Paths", outer = T, cex = 1.25)
dev.off()



























