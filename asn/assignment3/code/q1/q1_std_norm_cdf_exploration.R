##
## Script for plotting the normal CDF approximation used in the C++ 
## implementation for the implied volatility calculator
##

setwd("~/drive/concordia/2015-2016/1_fall2015/macf402/assignment3/")

std_norm_pdf <- function(x) {
  1/sqrt(2 * pi) * exp(-1/2 * x^2)
}
std_norm_cdf <- function(x) {
  b0 <- 0.2316419
  b1 <- 0.31938530
  b2 <- -0.356563782
  b3 <- 1.781477937
  b4 <- -1.821255978
  b5 <- 1.330274429
  t <- 1/(1 + b0 * x)
  
  y <- 1 - std_norm_pdf(x) * (b1*t + b2*t^2 + b3*t^3 + b4*t^4 + b5*t^5)
  return (y)
}

pdf("./plots/q1/naive_approx1.pdf", height = 4.75, width = 9.75)
par(mfrow = c(1,2))
x <- seq(-3,3,0.1)
plot(std_norm_cdf(x)~x, type = 'p', lwd = 2, pch = 1,
     ylab = 'Standard Normal CDF Output')
legend('bottomright', legend = c(expression('NaiveApprox', paste(Phi,"(x)"))),
       lty = c('dashed','solid'), pch = c(1, 46), lwd = c(2,2))
lines(pnorm(x)~x, lwd = 2)

x <- seq(-4,-3,0.01)
plot(std_norm_cdf(x)~x, type = 'p', lwd = 2, pch = 1,
     ylab = 'Standard Normal CDF Output')
legend('bottomright', legend = c(expression('NaiveApprox', paste(Phi,"(x)"))),
       lty = c('dashed', 'solid'), pch = c(1, 46), lwd = c(2,2))
lines(pnorm(x)~x, lwd = 2)
dev.off()

pdf("./plots/q1/naive_approx_error.pdf", height = 4.75, width = 4.75)
par(mfrow = c(1,1))
x <- seq(-7,7,0.01)
norm_err <- abs(std_norm_cdf(x) - pnorm(x))
plot(log10(norm_err)~x, type = 'l', lwd = 2,
     ylab = 'Log10 Absolute Error')
dev.off()

std_norm_cdf_TWO <- function(x_vect) { # vectorized/R optimized solution
  # we notice that the approximation is much better for positive x
  # we use symmetry of the normal CDF to compute the negative x values
  
  x <- c(-x_vect[which(x_vect < 0)], x_vect[which(x_vect >= 0)])
  
  b0 <- 0.2316419
  b1 <- 0.31938530
  b2 <- -0.356563782
  b3 <- 1.781477937
  b4 <- -1.821255978
  b5 <- 1.330274429
  t <- 1/(1 + b0 * x)
  
  y <- 1 - std_norm_pdf(x) * (b1*t + b2*t^2 + b3*t^3 + b4*t^4 + b5*t^5)
  y <- c(1 - y[which(x_vect < 0)], y[which(x_vect >= 0)])
  
  return (y)
}

pdf("./plots/q1/lessnaive_approx_error.pdf", height = 4.75, width = 9.75)
par(mfrow = c(1,2))
x <- seq(-7,7,0.1)
plot(std_norm_cdf_TWO(x)~x, type = 'p', pch = 1,
     ylab = 'Standard Normal CDF Output')
legend('topleft', legend = c(expression('LessNaiveApprox', paste(Phi,"(x)"))),
       lty = c('dashed', 'solid'), pch = c(1, 46), lwd = c(2,2))
lines(pnorm(x)~x, lwd = 2)

norm_err_TWO <- abs(std_norm_cdf_TWO(x) - pnorm(x))
plot(log10(norm_err_TWO)~x, type = 'l', lwd = 2,
     ylab = 'Log10 Absolute Error')
dev.off()



