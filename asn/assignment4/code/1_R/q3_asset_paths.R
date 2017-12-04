setwd("~/drive/concordia/2015-2016/1_fall2015/macf402/assignment4/")
set.seed(402)
N <- 1e5
nplot <- 20

S0 <- 500
K <- 500
r <- 0.0175
sigma <- 0.25
Tend <- 1
n <- 52
CI <- 0.95

dt <- Tend/n
times <- seq(0, 1, dt)

paths <- matrix(nrow = n, ncol = N)
paths <- replicate(N, {
  Z <- rnorm(n)
  c(S0, S0 * cumprod( exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z) ))
})

pctls <- apply(paths, 1, quantile, probs = c(1 - CI, CI))

pdf("./plots/q3/path_sim.pdf", height = 4, width = 7)
plot(paths[,1] ~ times, type = 'l', ylim = c(min(paths[,1:nplot]), max(paths[,1:nplot])), 
     col = 'gray70', main = 'Lognormal Asset Simulation', xlab = 'Time', ylab = 'Price')
for (i in 2:nplot) {
  lines(paths[,i] ~ times, col = 'gray70')
}
lines(rowMeans(paths) ~ times, lwd = 2)
lines(pctls[1,] ~ times, lwd = 2)
lines(pctls[2,] ~ times, lwd = 2)
dev.off()















