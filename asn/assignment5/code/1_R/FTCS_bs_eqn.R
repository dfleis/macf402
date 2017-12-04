#### set working directory
setwd("~/drive/concordia/2015-2016/1_fall2015/macf402/assignment5/")

#### functions
compute_d1 <- function(S, K, r, sigma, dt) {
  ( log(S / K) + (r + 0.5 * sigma^2) * dt ) / ( sigma * sqrt(dt) )
}
compute_d2 <- function(S, K, r, sigma, dt) {
  compute_d1(S, K, r, sigma, dt) - sigma * sqrt(dt)
}
bs_price_put <- function(S, K, r, sigma, dt) {
  if (dt == 0) {
    return( pmax(K - S, 0) )
  }
  d1 <- compute_d1(S, K, r, sigma, dt)
  d2 <- compute_d2(S, K, r, sigma, dt)
  pnorm(-d2) * K * exp(-r * dt) - pnorm(-d1) * S
}

#### parameters
S0 <- 100
K <- 100
sigma <- 0.2
r <- 0.03
t_max <- 1
S_max <- 3 * K 

ftcs_method <- function(N, J) {
  dS <- S_max/J
  dt <- t_max/N
  nu <- dt/dS^2
  
  # Create the tridiagonal matrix
  j_vec <- 1:(J - 1)
  aa <- (-0.5 * r * dt * j_vec + 0.5 * sigma^2 * dt * j_vec^2) / (1 + r * dt)
  bb <- (1 - sigma^2 * dt * j_vec^2) / (1 + r * dt)
  cc <- (0.5 * r * dt * j_vec + 0.5 * sigma^2 * dt * j_vec^2) / (1 + r * dt)
  A <- diag(bb) + cbind(rbind(0, diag(aa[-1])),0) + cbind(0,rbind(diag(head(cc,length(cc)-1)),0))
  
  V <- matrix(0, nrow = N + 1, ncol = J + 1)
  V[N + 1,] <- pmax(K - seq(0, S_max, dS), 0) # initial condition
  V[, 1] <- exp(-r * seq(t_max, 0, -dt)) * K # boundary condition for S -> 0
  
  
  for (i in N:1) { # FTCS
    p <- rep(0, J - 1)
    p[1] <- aa[1] * V[i + 1, 1] # boundary condition for S -> 0
    V[i, 2:J] <- A %*% V[(i + 1), 2:J] + p
  }
  
  # create coordinate vectors
  taxt <- seq(0, t_max, dt)
  Saxt <- seq(0, S_max, dS)
  
  # Construct the true solution
  V_true <- matrix(nrow = nrow(V), ncol = ncol(V))
  for (n in 1:nrow(V_true)) { # time
    V_true[n,] <- bs_price_put(Saxt, K, r, sigma, t_max - taxt[n])
  }
  
  #### PLOTS 
  plot_filepath <- paste0("./plots/q5_J", J, "_N", N, "_")
  ncols <- 128
  colfunc <- colorRampPalette(c('black','gray80'))
  V_cols <- colfunc(ncols)[as.numeric(cut(V, breaks = ncols))]
  V_true_cols <- colfunc(ncols)[as.numeric(cut(V_true, breaks = ncols))]
  border_col <- NA
  
  # plot FTCS solution
  pdf(paste0(plot_filepath, "ftcs.pdf"), height = 8, width = 8)
  par(mfrow = c(1,1), oma=c(0,0,2,0))
  persp(taxt, Saxt, V, theta = 140, phi = 10, ticktype = "detailed", col = V_cols,
        ylab = "Spot", xlab = "Time", zlab = "Call Price", border = border_col)
  mtext(paste("FTCS Solution\n N Spot Steps = ", N, "\nJ Time Steps = ", J), cex = 1.25)
  dev.off()
  
  ## plot black-scholes solution
  #pdf(paste0(plot_filepath, "bs.pdf"), height = 8, width = 8)
  #par(mfrow = c(1,1), oma=c(0,0,2,0))
  #persp(taxt, Saxt, V_true, theta = 140, phi = 10, ticktype = "detailed", col = V_true_cols,
  #      ylab = "Spot", xlab = "Time", zlab = "Call Price", border = border_col)
  #mtext("Black-Scholes Solution\n\n", cex = 1.25)
  #dev.off()
  
  ## plot errors
  #pdf(paste0(plot_filepath, "err1.pdf"), height = 8, width = 8)
  #par(mfrow = c(1,1), oma=c(0,0,2,0))
  #persp(taxt, Saxt, V_true - V, theta = 120, phi = 25, ticktype = "detailed", col = "gray80",
  #      ylab = "Spot", xlab = "Time", zlab = "Error", border = "gray40")
  #mtext(paste("FTCS Error\nN Spot Steps = ", N, "\nJ Time Steps = ", J), cex = 1.25)
  #dev.off()
  #pdf(paste0(plot_filepath, "err2.pdf"), height = 8, width = 8)
  #par(mfrow = c(1,1), oma=c(0,0,2,0))
  #persp(taxt, Saxt, V_true - V, theta = -120, phi = 25, ticktype = "detailed", col = "gray80",
  #      ylab = "Spot", xlab = "Time", zlab = "Error", border = "gray40")
  # mtext(paste("FTCS Error\nN Spot Steps = ", N, "\nJ Time Steps = ", J), cex = 1.25)
  #dev.off()
  
  mx_abs_err <- max(abs(V - V_true))
  return (mx_abs_err)
}

N <- seq(10, 450, 5) # nb of t steps
J <- seq(10, 120, 2) # nb of x steps
E <- matrix(nrow = length(N), ncol = length(J)) # rows are time, cols are space
for (n in 1:length(N)) {
  for (j in 1:length(J)) {
    E[n,j] <- ftcs_method(N[n], J[j])
    print(c(N[n], "---", J[j]))
  }
}

## plot errors
pdf("./plots/q5_max_abs_err.pdf", height = 8, width = 9)
par(mfrow = c(1,1), oma=c(0,0,2,0))
persp(N, J, log10(E), theta = -50, phi = 35, ticktype = "detailed",
      zlab = "Log10 Error", col = "gray80", border = "gray40")
mtext("Max Absolute Error\nFTCS - Black-Scholes", cex = 1.25)
dev.off()

dat_log10_err <- as.data.frame(log10(E))
dat_log10_err <- round(dat_log10_err, 2)

colnames(dat_log10_err) <- paste0("J=", J)
rownames(dat_log10_err) <- paste0("N=", N)
write.csv(x = dat_log10_err[,1:floor(length(J)/2)], file = "./data/q5_dat_err_a.csv", quote = F)
write.csv(x = dat_log10_err[,(floor(length(J)/2) + 1):length(J)], file = "./data/q5_dat_err_b.csv", quote = F)

## plot only region of stability
h <- log10(E) 
h[h > 0] <- NA # replace all errors > 1 => log10 error > 0 with NA

pdf("./plots/q5_max_abs_err_trunc.pdf", height = 8, width = 9)
par(mfrow = c(1,1), oma=c(0,0,2,0))
persp(N, J, h, theta = 220, phi = 20, ticktype = "detailed",
      zlab = "Log10 Error", col = "gray80", border = "gray40")
mtext("Max Absolute Error\nFTCS - Black-Scholes\nLog10 Error < 0", cex = 1.25)
dev.off()

## find boundary of stability wrt J
boundary <- apply(E, 2, FUN = function(x) { min(which(x <= 1)) })
boundary[boundary == Inf] <- NA

# fit to linear model
model <- lm(sqrt(boundary)~J)
int <- coef(model)["(Intercept)"]
slope <- coef(model)["J"]
model_sq <- (int + slope * J)^2

pdf("./plots/q5_model.pdf", height = 6, width = 5.5)
par(mfrow = c(1,1), oma=c(0,0,2,0))
plot(boundary~J, type = 'p', lwd = 2,
     ylab = "N")
lines(model_sq~J, lwd = 3, col = 'gray30')
text(x = 40, y = 60, substitute(paste(sqrt(hat(N)), " = ", b0, " + ", b1, "J"), 
                                list(b0 = round(int,4), b1 = round(slope,4))))
mtext("Least N for Stability\n(Log10 Error < 0)\n", cex = 1.25)
dev.off()






