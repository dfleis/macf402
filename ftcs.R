# FTCS Example for canonical heat equation
# u_t = u_xx
# u(0,x) = lambda(x)
# u(t,x_L) = alpha(t)
# u(t,x_R) = beta(t)
#
# where we assume that x_L = 0, x_R = pi
# lambda(x) = sin(x)
# alpha(t) = 0, beta(t) = 0
#
#
# can show analytically that the solution to the above problem is
# u(x,t) = e^(-t)sin(x)
#
# implement the FTCS method and compare to the analytic (true) solution
xL <- 0
xR <- pi
T <- 3


N <- 476
J <- 28

dt <- T/N
dx <- (xR - xL)/J
v <- dt/(dx^2)

# We define U on the interior of the space domain
Uint <- matrix(0, nrow = (J - 1), ncol = (N + 1))

# We add the initial condition at time t=0
for (i in 1:(J - 1)) {
    Uint[i, 1] <- sin(i*dx)
}

aa <- v
bb <- 1 - 2*v
cc <- v
J2 <- J - 1

# Create the tridiagonal matrix tridiag(aa,bb,cc)
# Note that the A here is the (I-vA) in the notes
A <- matrix(nrow = J2, ncol = J2, data = 0)
diag(A) <- bb
A[cbind(1:(J2 - 1), 2:J2)] <- cc
A[cbind(2:J2, 1:(J2 - 1))] <- aa

## Apply the FTCS
for (n in 1:N){
    Uint[,n + 1] <- A %*% Uint[,n]
}

# We add the boundary conditions at x_L and x_R
C <- rbind(0, Uint, 0)

# Plot the approximate solution
t <- seq(0, T, by = dt)
x <- seq(xL, xR, by = dx)
persp(x, t, C, theta = 45, phi = 15)

# Construct the true solution
N <- length(t)
J <- length(x)
C_true <- matrix(nrow = J, ncol = N, data = 0)

for (n in 1:N){
    for (j in 1:J){
    	C_true[j, n] <- exp(-t[n])*sin(x[j])
	}
}

# Compare the approximate solution C from FTCS and the true solution C_true 
persp(x, t, C - C_true, theta = 45, phi = 15)

# Calculate the maximum error
max(abs(C - C_true))

## For J=14 and N=199 we find the maximum error is
## 0.0012344
## For J=28 and N=199 we find that the maximum error is
## 1.89022e+97
## because the stability condition on dx and dt is not satisfied so the errors explode
## In order to double the number of space points we need to ensure v <= 0.5
## with J=28 this means we must have N >= 476.6148
## If we take J=28 and N=477 then we find that the maximum error is 
## 0.0007727066

