rm(list = ls())
setwd('/project/jointage/matt/DM/paper/')
source('/project/jointage/matt/DM/paper/Functions.R')
ToeplitzInv <- function(n, c) {
  #compute the inverse of the tri-diagonal Toeplitz matrix with constant diagonal 
  #entries c and off-diagonals -1 using the formula based on roots r and s 
  #of x^2 - c*x + 1 = 0. note, c must be greater than 2.
  
  zz <- c^2 - 4 
  
  r <- (c + sqrt(zz)) / 2
  s <- (c - sqrt(zz)) / 2
  
  denom <- r^(n+1) - s^(n+1)
  
  if (denom==0) {
    stop("singular matrix")
  }
  
  T_inv <- matrix(0, nrow=n, ncol=n)
  for (i in 1:n) {
    for (j in 1:n) {
      m <- min(i, j)
      M <- max(i, j)
      num <- (r^m - s^m) * (r^(n+1-M) - s^(n+1-M))
      T_inv[i, j] <- num / ((r - s) * denom)
    }
  }
  
  return(T_inv)
}

#example for inverting D matrix
n <- 10
sigma_sq <- 1
sigma_sq_ep <- 1
psi <- sigma_sq /sigma_sq_ep
c <- 2 + 1/psi
T_inv_nonrec <- ToeplitzInv(n, c)
print(T_inv_nonrec)s
#compare to reference
ref <- toeplitz(c(c, -1, rep(0, n-2)))
print(solve(ref))

