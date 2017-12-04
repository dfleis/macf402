n <- 3

C <- c(12,4,NA)
B <- c(2,5,1)
A <- c(NA,5,7)
D <- c(4, 9, 99)

M <- matrix(nrow = n, ncol = n)
M[1,] <- c(B[1], C[1], 0)
M[2,] <- c(A[2], B[2], C[2])
M[3,] <- c(0, A[3], B[3])

M
x <- vector(length = n, mode = 'numeric')

for (k in 2:n) {
  m <- A[k] / B[k - 1]
  print(m)
  B[k] <- B[k] - m * C[k - 1]
  D[k] <- D[k] - m * D[k - 1]
}

x[n] <- D[n] / B[n]

for (k in (n - 1):1) {
  x[k] <- (D[k] - C[k] * x[k + 1]) / B[k]
}

x
M %*% x
