#### functions
compute_d1 <- function(S, K, r, sigma, dt) {
  ( log(S / K) + (r + 0.5 * sigma^2) * dt ) / ( sigma * sqrt(dt) )
}
compute_d2 <- function(S, K, r, sigma, dt) {
  compute_d1(S, K, r, sigma, dt) - sigma * sqrt(dt)
}
bs_price_put <- function(S, K, r, sigma, dt) {
  d1 <- compute_d1(S, K, r, sigma, dt)
  d2 <- compute_d2(S, K, r, sigma, dt)
  pnorm(-d2) * K * exp(-r * dt) - pnorm(-d1) * S
}
tridiag <- function(a, b, c, n) {
  v <- c(a,b,c)
  B <- matrix(nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      B[i,j] <- 0
    }
  }
  
  B[1,1] <- v[2]
  B[1,2] <- v[3]
  
  for (i in 2:(n - 1)) {
    B[i, i - 1] <- v[1]
    B[i, i] <- v[2]
    B[i, i + 1] <- v[3]
  }
  
  B[n, n - 1] <- v[1]
  B[n, n] <- v[2]
  return (B)
}

#### parameters
S0 <- 100
K <- 100
sigma <- 0.2
r <- 0.03
t_max <- 10
S_min <- 0 # min spot price
S_max <- 3 * K # max spot price
rho <- -0.5 * (2  * r / sigma^2 - 1) 
xi <- -0.25 * (2 * r / sigma^2 + 1)^2

N <- 10 # nb of t steps
J <- 10 # nb of x steps

dS <- (S_max - S_min)/J
dt <- t_max/N
nu <- dt/dS^2

# define U on the interior of the space domain, but the entire time domain
Uint <- matrix(0, nrow = (J - 1), ncol = (N + 1))

# add the initial condition along time t = 0 over the interior of the space domain
for (j in 1:(J - 1)) {
  Uint[j,1] <- exp(-rho * log(j * dS)) * max(K - j * dS, 0)
}
# add the boundary conditions over the entire time domain
alpha <- vector(mode = 'numeric', length = N + 1)
for (n in 1:(N + 1)) {
  alpha[n] <- K * exp(-r * (n - 1) * dt )
}

# create the tridiagonal matrix
A <- tridiag(nu, 1 - 2 * nu, nu, J - 1)

# apply the FTCS 
for (n in 1:N) { # time march
  p <- c(nu * alpha[n], rep(0, J - 2))
  Uint[,n + 1] <- A %*% Uint[,n] + p
}

# transform interior back to V
Vint <- matrix(nrow = nrow(Uint), ncol = ncol(Uint))
for (n in 1:(N + 1)) {
  tau <- (n - 1) * dt * sigma^2 / 2
  for (j in 1:(J - 1)) {
    Vint[j,n] <- exp(rho * log(j * dS) + xi * tau) * Uint[j,n]
  }
}
# add the boundary conditions at S_min and S_max to create the full grid
V <- rbind(alpha, Vint, 0)

# create coordinate vectors
taxis <- seq(0, t_max, dt)
Saxis <- seq(S_min, S_max, dS)

# construct the true solution
V_true <- matrix(nrow = nrow(V), ncol = ncol(V))
for (n in 1:ncol(V_true)) { # time
  for (j in 1:nrow(V_true)) { # space 
    V_true[j,n] <- bs_price_put(Saxis[j], K, r, sigma, t_max - taxis[n])
  }
}
# flip time axis of true solution to match time vector
for (j in 1:nrow(V_true)) {
  V_true[j,] <- rev(V_true[j,])
}

# plot solution
persp(Saxis, taxis, V, theta = 45, phi = 15, ticktype = "detailed")
# plot true solution
 persp(Saxis, taxis, V_true, theta = 45, phi = 15, ticktype = "detailed")
# Compare with Black-Scholes price
# persp(Saxis, taxis, V_true - V, theta = 45, phi = 15, ticktype = "detailed")

# Calculate the maximum error
max(abs(V-V_true))













