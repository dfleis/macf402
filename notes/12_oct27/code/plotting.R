set.seed(1238)
setwd('~/Drive/concordia/2015-2016/1_fall2015/macf402/notes/12_oct27/')

####
#### PRICING A VANILLA PUT
####
powmin <- 1; powmax <- 5
N <- 10^seq(powmin, powmax, 1/3)
pctl <- 0.95

dat <- data.frame(matrix(nrow = length(N), ncol = 3))
rownames(dat) <- N
colnames(dat) <- c('mean','CI_lo','CI_hi')

S_0 <- 10
K <- 9
sigma <- 0.1
r <- 0.06
T_end <- 1

Z_alpha <- qnorm(1 - (1 - pctl)/2)
df <- exp(-r * T_end)

for (i in 1:length(N)) {
  W_T <- sqrt(T_end) * rnorm(N[i])
  S_T <- S_0 * exp((r - 1/2 * sigma^2) * T_end + sigma * W_T)
  d_payoff <- df * pmax(K - S_T, 0)
  dat$mean[i] <- mean(d_payoff)
  dat$CI_hi[i] <- dat$mean[i] + Z_alpha * sd(d_payoff)/sqrt(N[i])
  dat$CI_lo[i] <- dat$mean[i] - Z_alpha * sd(d_payoff)/sqrt(N[i])
}
wt_avg <- sum(dat$mean * N)/sum(N)

pdf(file = './plots/putapprox.pdf', height = 5, width = 6.5) 
plot(x = N, dat$mean, pch = 4, col = 'blue',
     ylim = c(min(dat),max(dat)),
     xaxt = "n", xlab = "Num Samples", ylab = 'Put Option Value',
     log = 'x')
ticksat <- as.vector(sapply(seq(powmin, powmax), function(p) (1:10)*10^p))
axis(1, 10^seq(powmin, powmax), tcl = 0.5, lwd.ticks = 1)
axis(1, ticksat, labels = NA, tcl = 0.25, lwd.ticks = 1)
abline(h = wt_avg, col = 'red', lty = 'longdash')

for (i in 1:length(N)) {
  lines(x = c(N[i], N[i]), y = c(dat$CI_lo[i], dat$CI_hi[i]), col = 'blue')
}
dev.off()

####
#### PRICING A DOWN & OUT CALL
####
N <- 10^seq(2,6)

S_0 <- 100
K <- 110
b <- 95
r <- 0.05
sigma <- 0.2
T_end <- 1
s <- 52 # presumably number of steps for GBM path
pctl <- 0.95

dt <- 1/s
Z_alpha <- qnorm(1 - (1 - pctl)/2)

dat <- data.frame(matrix(nrow = length(N), ncol = 3))
colnames(dat) <- c('mean', 'CI_lo', 'CI_hi')
rownames(dat) <- N

for (i in 1:length(N)) {
  sim <- replicate(N[i], {
    Z <- rnorm(s)
    B <- sqrt(dt) * Z
    S <- S_0 * c(1,cumprod(exp(( r - 1/2 * sigma^2) * dt + sigma * B)))
    if (min(S) < b) {
      return (0)
    } else {
      return (max(S[s + 1] - K, 0))
    }
  })

  dat$mean[i] <- mean(sim)
  dat$CI_hi[i] <- dat$mean[i] + Z_alpha * sd(sim)/sqrt(N[i])
  dat$CI_lo[i] <- dat$mean[i] - Z_alpha * sd(sim)/sqrt(N[i])
} 
wt_avg <- sum(dat$mean * N)/sum(N)

pdf(file = './plots/d_o_callapprox.pdf', height = 5, width = 6.5) 
plot(x = N, dat$mean, pch = 4, col = 'blue',
     ylim = c(min(dat),max(dat)),
     xlab = "Num Samples", ylab = 'Down & Out Call Option Value',
     log = 'x')
abline(h = wt_avg, col = 'red', lty = 'longdash')

for (i in 1:length(N)) {
  lines(x = c(N[i], N[i]), y = c(dat$CI_lo[i], dat$CI_hi[i]), col = 'blue')
}
dev.off()

####
#### Showing discretization error
####
lambda <- 2
mu <- 1
X_0 <- 1
T_end <- 1
delta_t_max <- 2^-16
delta_t <- c(4*2^-8, 2*2^-8, 2^-8)

N_max <- T_end/delta_t_max
tau_max <- seq(0,1,delta_t_max)

Z_max <- rnorm(N_max)
W_max <- sqrt(delta_t_max) * Z_max
X_t_max <- X_0 * c(1, cumprod(exp( (lambda - 1/2 * mu^2) * delta_t_max + mu * W_max)))

i <- 3

x_sample <- c(1,seq(0,N_max,delta_t[i]/delta_t_max)[-1])
tau <- x_sample/max(x_sample)

Z <- Z_max[x_sample][-1]
W <- sqrt(delta_t[i]) * Z
X_t <- X_0 * c(1,cumprod(exp( (lambda - 1/2 * mu^2) * delta_t[i] + mu * W)))

plot(x = tau_max, y = X_t_max, type = 'l', col = 'purple',
     ylim = c(min(X_t,X_t_max), max(X_t, X_t_max)))
lines(x = tau , y = X_t, col = 'red', type = 'l', lty = 'dashed')

X_t
  
  
  
  
  
  
