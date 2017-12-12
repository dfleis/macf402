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

N <- 100
Jmin <- 10
Jmax <- 30
Js <- Jmin:Jmax
errvec_J <- vector(mode = 'numeric', length = Jmax - Jmin + 1)

for (J in Js) {
  
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
  
  # Construct the true solution
  N <- length(t)
  J3 <- length(x)
  C_true <- matrix(nrow = J3, ncol = N, data = 0)
  
  for (n in 1:N){
    for (j in 1:J){
      C_true[j, n] <- exp(-t[n])*sin(x[j])
    }
  }
  
  # Calculate the maximum error
  errvec_J[J - Jmin + 1] <- max(abs(C - C_true))
}

plot(errvec_J ~ Js, type = 'o', pch = 19)
plot(log10(errvec_J) ~ Js, type = 'o', pch = 19)
