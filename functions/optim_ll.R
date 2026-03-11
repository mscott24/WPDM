#estimation of degradation models using joint loglikelihood
run_optim_ll <- function(df, modtype, params0, ...) { 
  
  if (!all(c("V", "tau", "Patient") %in% names(df))) {
    stop('Missing columns (V, tau, or Patient) in df')
  }
  
  if (modtype=='m3') {
    
    if (missing(params0)) {
      params0 <- c(0, 0.1, 0.1, 0.1)
    }
    
    result <- optim(fn = ll,
                    par = params0,
                    df = df, 
                    method = "L-BFGS-B",
                    lower = c(-Inf, 1e-4, 1e-4, 1e-4), 
                    upper = c(Inf, Inf, Inf, Inf),
                    control = list(fnscale = -1))
    stats <- list()
    stats$mu <- result$par[1]
    stats$sigma_sq <- result$par[2]
    stats$sigma_sq_mu <- result$par[3]
    stats$sigma_sq_ep <- result$par[4]
  } else if (modtype=='m2') {
    
    if (missing(params0)) {
      params0 <- c(0, 0.1, 0.1)
    }
    
    result <- optim(fn = ll,
                    par = params0,
                    df = df, 
                    method = "L-BFGS-B",
                    lower = c(-Inf, 1e-4, 1e-4), upper = c(Inf, Inf, Inf),
                    control = list(fnscale = -1))
    stats <- list()
    stats$mu <- result$par[1]
    stats$sigma_sq <- result$par[2]
    stats$sigma_sq_mu <- result$par[3]
  } else if (modtype=='m1') {
    
    if (missing(params0)) {
      params0 <- c(0, 0.1)
    }
    
    result <- optim(fn = ll,
                    par = params0,
                    df = df, 
                    method = "L-BFGS-B",
                    lower = c(-Inf, 1e-4), upper = c(Inf, Inf),
                    control = list(fnscale = -1))
    stats <- list()
    stats$mu <- result$par[1]
    stats$sigma_sq <- result$par[2]
  } else {
    stop('model type not specified.')
  }
  
  return(stats)
}
ll <- function(params, df, ...) {
  
  if (length(params) < 4) {params[(length(params)+1):4] <- 0}
  
  patient_df <- split(df, df$Patient, drop = TRUE)
  func <- lapply(patient_df, function(data) {
    tau <- as.numeric(data$tau)
    V <- as.numeric(data$V)
    n <- length(tau)
    S <- diag(tau, nrow = n, ncol = n)
    A_first_row <- c(2, -1, rep(0, max(0, n-2)))
    A <- toeplitz(A_first_row[1:n])
    
    Sigma <- params[2] * S + 
      params[3] * tcrossprod(tau) + 
      params[4] * A
    
    mu <- params[1] * tau
    e <- V - mu
    
    ll <- as.numeric(
      -(n/2) * log(2*pi) -
        0.5 * log(det(Sigma)) -
        0.5 * t(e) %*% solve(Sigma) %*% e
    )
    
    ll
  })
  
  ll <- sum(sapply(func, function(x) x))
  
  return(ll)
}