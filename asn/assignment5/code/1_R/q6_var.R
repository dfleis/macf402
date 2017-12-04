#### functions
compute_d1 <- function(S, K, r, sigma, dt) {
  ( log(S / K) + (r + 0.5 * sigma^2) * dt ) / ( sigma * sqrt(dt) )
}
compute_d2 <- function(S, K, r, sigma, dt) {
  compute_d1(S, K, r, sigma, dt) - sigma * sqrt(dt)
}
bs_price_call <- function(S, K, r, sigma, dt) {
  if (dt == 0) {
    return( pmax(S - K, 0) )
  }
  d1 <- compute_d1(S, K, r, sigma, dt)
  d2 <- compute_d2(S, K, r, sigma, dt)
  pnorm(d1) * S - pnorm(d2) * K * exp(-r * dt) 
}

#### parameters
set.seed(26)
N <- 10^(2:6)

S0 <- 100
K <- 100
mu <- 0.02
sigma <- 0.5
t_max <- 1
h <- 5/365
z_alpha <- 1.645

short <- vector(mode = 'numeric', length = length(N))
short_ci_lo <- vector(mode = 'numeric', length = length(N))
short_ci_hi <- vector(mode = 'numeric', length = length(N))
long <- vector(mode = 'numeric', length = length(N))
long_ci_lo <- vector(mode = 'numeric', length = length(N))
long_ci_hi <- vector(mode = 'numeric', length = length(N))

C0 <- bs_price_call(S0, K, mu, sigma, t_max)
d1 <- compute_d1(S0, K, mu, sigma, t_max)
delta0 <- pnorm(d1) * S0
var <- z_alpha * sigma * sqrt(h) * delta0

for (i in 1:length(N)) {
  n <- N[i]
  print(n)
  Z <- rnorm(n)
  St <- S0 * exp((mu - 0.5 * sigma^2) * h + sigma * sqrt(h) * Z)
  Ct <- bs_price_call(St, K, mu, sigma, t_max - h)
  deltaC <- Ct - C0
  
  p_s <- sum(deltaC > var)/n
  p_l <- sum(deltaC < -var)/n
  short[i] <- p_s 
  short_ci_lo[i] <- p_s - z_alpha * sqrt(p_s * (1 - p_s) / n)
  short_ci_hi[i] <- p_s + z_alpha * sqrt(p_s * (1 - p_s) / n)
  long[i] <- p_l
  long_ci_lo[i] <- p_l - z_alpha * sqrt(p_l * (1 - p_l) / n)
  long_ci_hi[i] <- p_l + z_alpha * sqrt(p_l * (1 - p_l) / n)
}

breaks <- as.data.frame(cbind(N, short, short_ci_lo, short_ci_hi, long, long_ci_lo, long_ci_hi))

ymin <- min(breaks[,-1])
ymax <- max(breaks[,-1])

plot(breaks$short~breaks$N, log = 'x', type = 'b', ylim = c(ymin, ymax), lwd = 2)
segments(x0 = N, y0 = breaks$short_ci_lo, x1 = N, y1 = breaks$short_ci_hi, lwd = 2)
lines(breaks$long~breaks$N, type = 'b', pch = 2, lwd = 2)
segments(x0 = N, y0 = breaks$long_ci_lo, x1 = N, y1 = breaks$long_ci_hi, lwd = 2)
legend('topright', c("Short", "Long"), pch = c(1, 2), lwd = 2)










