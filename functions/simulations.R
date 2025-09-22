#data generation tools
DefaultSetup <- function(...) {
  setup <- list()
  setup$n_sims <- 1000  #number of simulations     
  setup$n_subjs <- 100 #number of subjects
  setup$alpha <- 0.05 #coverage probability
  setup$mu <- -1 #slope mean
  setup$sigma_sq <- 0.25 #brownian motion
  setup$sigma_sq_mu <- 0.25 #slope variance
  setup$sigma_sq_ep <- 0.25 #measurement error
  setup$psi <- setup$sigma_sq_ep / setup$sigma_sq #variance ratio for profile likelihood
  setup$fu_lb <- 3 #lower bound for number of follow-up visits in unequal case
  setup$fu_ub <- 10 #upper bound for number of follow-up visits  in unequal case
  setup$n <- 10 #equal FU case
  setup$home_dir <- '/project/jointage/matt/DM/WPDM/data'
  setup$batch_size <- min(100, ceiling(setup$n_sims / (2 * parallel::detectCores())))
  setup$u_lb <- 0.5
  setup$u_ub <- 1.5
  setup$p <- 0
  setup$tau <- 0.05
  return(setup)
}
GenTimes <- function(bounds, lambda, p, tau, ...) {
  if (missing(p)) {p=0}
  if (missing(tau)) {tau=0.05}
  
  if (missing(lambda) & length(bounds)==1) {
    #equal FU case
    return(0:(bounds[1])) 
  }
  
  else {
    #different number of FU visits and duration
    n_samp <- sample((bounds[1]+1):(bounds[2]+1), size = 1)
    t <- cumsum(runif(n_samp, min = 0.5, max = 1.5))
    if (p>0) {
      ind_close <- which(runif(n_samp) < p)
      for (i in ind_close) {
        t[i + 1] <- t[i] + runif(1, 0, tau)
      }
    }
    t <- t - t[1] 
    return(t)
  }
}
GenData <- function(setup, sim_id, opt, ...) {
  if (missing(opt)) {opt <- 'all'}
  
  #equal FU visits and duration
  ff_eq <- paste(sim_id, '_eq.RData', sep='')
  if (opt %in% c('all', 'eq')) {
    print('Generating eq')
    set.seed(2024)
    time_eq <- matrix(rep(t(GenTimes(bounds=setup$n)), setup$n_subjs), 
                      nrow=(setup$n+1), ncol=setup$n_subjs)
    eq <- list()
    eq$setup <- setup
    sims <- MakeSim(setup=setup, time=time_eq)
    eq$sim_m1 <- sims$sim_m1
    eq$sim_m2 <- sims$sim_m2
    eq$sim_m3 <- sims$sim_m3
    eq$sim_m1$dfs <- lapply(eq$sim_m1$dfs, as.data.table)
    eq$sim_m2$dfs <- lapply(eq$sim_m2$dfs, as.data.table)
    eq$sim_m3$dfs <- lapply(eq$sim_m3$dfs, as.data.table)
    save(eq, file=ff_eq)
    rm(eq)
  }
  
  #unequal FU time and number of visits
  ff_ue <- paste(sim_id, '_ue.RData', sep='')
  if (opt %in% c('all', 'ue')) {
    print('Generating ue')
    set.seed(2024)
    time_ue <- PadList(replicate(setup$n_subjs, 
                                 GenTimes(bounds=c(setup$fu_lb, setup$fu_ub), lambda=setup$lambda),
                                 simplify = FALSE))
    ue <- list()
    ue$setup <- setup
    sims <- MakeSim(setup=setup, time=time_ue)
    ue$sim_m1 <- sims$sim_m1
    ue$sim_m2 <- sims$sim_m2
    ue$sim_m3 <- sims$sim_m3
    ue$sim_m1$dfs <- lapply(ue$sim_m1$dfs, as.data.table)
    ue$sim_m2$dfs <- lapply(ue$sim_m2$dfs, as.data.table)
    ue$sim_m3$dfs <- lapply(ue$sim_m3$dfs, as.data.table)
    save(ue, file=ff_ue)
    rm(ue)
  }
}
RunSim <- function(params, t, ...) {
  n <- length(t) - 1 
  q <- outer(t, t, FUN = pmin)
  iden <- diag(n+1) 
  Sigma <- params[2] * q +
    params[3] * (t %*% t(t))+ 
    params[4] * iden
  mu <- params[1] * t
  return(mvrnorm(mu=mu, Sigma=Sigma))
}
MakeSim <- function(setup, time, ...) {
  
  sims <- list()
  sims$sim_m3 <- list()
  sims$sim_m2 <- list()
  sims$sim_m1 <- list()
  sims$sim_m3$dfs <- list()
  sims$sim_m2$dfs <- list()
  sims$sim_m1$dfs <- list()
  
  params_m3 <- c(setup$mu,
                 setup$sigma_sq,
                 setup$sigma_sq_mu,
                 setup$sigma_sq_ep)
  params_m2 <- c(setup$mu,
                 setup$sigma_sq,
                 setup$sigma_sq_mu,
                 0)
  params_m1 <- c(setup$mu,
                 setup$sigma_sq,
                 0,
                 0)
  
  max_fu <- nrow(time)
  for (k in 1:setup$n_sims) {
    Y_k_m3 <- matrix(nrow=max_fu, ncol = setup$n_subjs)
    Y_k_m2 <- matrix(nrow=max_fu, ncol = setup$n_subjs)
    Y_k_m1 <- matrix(nrow=max_fu, ncol = setup$n_subjs)
    for (i in 1:setup$n_subjs) {
      time_i <- as.numeric(na.omit(time[,i]))
      Y_k_m3[1:length(time_i),i] <- RunSim(params = params_m3, t=time_i)
      Y_k_m2[1:length(time_i),i] <- RunSim(params = params_m2, t=time_i)
      Y_k_m1[1:length(time_i),i] <- RunSim(params = params_m1, t=time_i)
    }
    sims$sim_m3$dfs[[k]] <- FormatDiffDF(MakeLong(Y_k_m3, time), Yvar = 'Y', idvar = 'Patient', timevar = 't')
    sims$sim_m2$dfs[[k]] <- FormatDiffDF(MakeLong(Y_k_m2, time), Yvar = 'Y', idvar = 'Patient', timevar = 't')
    sims$sim_m1$dfs[[k]] <- FormatDiffDF(MakeLong(Y_k_m1, time), Yvar = 'Y', idvar = 'Patient', timevar = 't')
  }
  
  return(sims)
}

#joint likelihood estimators using LGBF-S from optim
mc_optim_ll  <- function(simdat, setup, modtype, ...) {
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  clusterExport(cl, c("ll", "run_optim_ll", 
                      "simdat", "modtype"), envir = environment())
  batch_size <- setup$batch_size
  ests <- NULL
  print(paste('simulations', modtype, ' optim running'))
  for (i in seq(1, setup$n_sims, by = batch_size)) {
    print(i)
    indices <- i:min(i + batch_size - 1, setup$n_sims)
    results <- pblapply(indices, function(idx) {
      run_optim_ll(simdat$dfs[[idx]], modtype)
    }, cl = cl)
    ests <- rbind(ests, do.call(rbind, results))
    gc()
  }
  stopCluster(cl)
  gc()
  
  stats <- list()
  stats$mu <- unlist(ests[,1])
  stats$sigma_sq <- unlist(ests[,2])
  
  if (modtype=='m2') {
    stats$sigma_sq_mu <- unlist(ests[,3])
  }
  
  if (modtype=='m3') {
    stats$sigma_sq_mu <- unlist(ests[,3])
    stats$sigma_sq_ep <- unlist(ests[,4])
  }
  return(stats)
}

#profile likelihood estimators - MLEs
mc_mle_m1 <- function(sim_m1, ...) {
  print('m1 mle simulations running')
  ests <- pblapply(sim_m1$df, mle_m1)
  stats_m1 <- list()
  stats_m1$mu <- sapply(ests, function(x) x$mu)
  stats_m1$sigma_sq <- sapply(ests, function(x) x$sigma_sq)
  return(stats_m1)
}
mc_mle_m2 <- function(sim_m2, setup, ...) {

  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  clusterExport(cl, c("mle_m2", "mu_mle_m2", 
                      "sim_m2", "profileL_mle_m2"), envir = environment())
  batch_size <- setup$batch_size
  ests <- NULL
  print('m2 mle simulations running')
  for (i in seq(1, setup$n_sims, by = batch_size)) {
    print(i)
    indices <- i:min(i + batch_size - 1, setup$n_sims)
    results <- pblapply(indices, function(idx) {
      mle_m2(df=sim_m2$dfs[[idx]])
    }, cl = cl)
    ests <- rbind(ests, do.call(rbind, results))
    gc()
  }
  stopCluster(cl)
  gc()
  stats_m2 <- list()
  stats_m2$mu <- unlist(ests[,1])
  stats_m2$sigma_sq <- unlist(ests[,2])
  stats_m2$sigma_sq_mu <- unlist(ests[,3])
  return(stats_m2)
}
mc_mle_m3 <- function(sim_m3, setup, ...) {
  
  print('m3 mle simulations running')
  
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  clusterExport(cl, c("mle_m3", "psi_ll_mle_m3", 
                      "sigma_sq_mu_mle_m3", "sigma_sq_mle_m3", "mu_mle_m3",
                      "sim_m3"), envir = environment())
  
  batch_size <- setup$batch_size
  ests <- list() 
  counter <- 1
  print(paste("batch size:", batch_size))
  for (i in seq(1, setup$n_sims, by = batch_size)) {
    print(paste("running batch starting at", i))
    indices <- i:min(i + batch_size - 1, setup$n_sims)
    
    results <- pblapply(indices, function(idx) {
      tryCatch({
        mle_m3(df = sim_m3$dfs[[idx]])
      }, error = function(e) {
        warning(sprintf("mle_m3 failed at idx=%d: %s", idx, e$message))
        return(rep(NA, 5))
      })
    }, cl = cl)
    
    for (j in seq_along(results)) {
      ests[[counter]] <- results[[j]]
      counter <- counter + 1
    }
    gc()
  }
  
  stopCluster(cl)
  gc()
  
  ests_mat <- do.call(rbind, ests)
  
  if (nrow(ests_mat) != setup$n_sims) {
    warning(sprintf("expected %d simulations, but got %d rows", setup$n_sims, nrow(ests_mat)))
  }
  
  stats_m3 <- list()
  stats_m3$psi  <- unlist(ests_mat[, 1])
  stats_m3$mu <- unlist(ests_mat[, 2])
  stats_m3$sigma_sq <- unlist(ests_mat[, 3])
  stats_m3$sigma_sq_mu <- unlist(ests_mat[, 4])
  stats_m3$sigma_sq_ep <- unlist(ests_mat[, 5])
  
  return(stats_m3)
}

#profile likelihood estimators - MLEs
mc_u_m2 <- function(sim_m2, ...) {
  print('m2 u simulations running')
  tmp <- lapply(sim_m2$df, u_m2)
  stats_m2 <- list()
  stats_m2$mu <- sapply(tmp, function(x) x$mu)
  stats_m2$sigma_sq <- sapply(tmp, function(x) x$sigma_sq)
  stats_m2$sigma_sq_mu <- sapply(tmp, function(x) x$sigma_sq_mu)
  return(stats_m2)
}
mc_u_m3 <- function(sim_m3, setup, ...) {
  print('m3 u simulations running')
  
  # Setup parallel cluster
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  clusterExport(cl, c("u_m3", "psi_ll_u_m3", 
                      "sigma_sq_mu_u_m3", "sigma_sq_u_m3", "mu_u_m3",
                      "sim_m3"), envir = environment())
  batch_size <- setup$batch_size
  ests <- list()
  counter <- 1
  print(paste("Batch size:", batch_size))
  for (i in seq(1, setup$n_sims, by = batch_size)) {
    print(paste("running batch starting at", i))
    indices <- i:min(i + batch_size - 1, setup$n_sims)
    
    results <- pblapply(indices, function(idx) {
      tryCatch({
        u_m3(sim_m3$df[[idx]])
      }, error = function(e) {
        warning(sprintf("u_m3 failed at idx=%d: %s", idx, e$message))
        return(rep(NA, 5))
      })
    }, cl = cl)
    
    for (j in seq_along(results)) {
      ests[[counter]] <- results[[j]]
      counter <- counter + 1
    }
    
    gc()
  }
  
  stopCluster(cl)
  gc()
  
  ests_mat <- do.call(rbind, ests)
  if (nrow(ests_mat) != setup$n_sims) {
    warning(sprintf("expected %d simulations, but got %d rows", setup$n_sims, nrow(ests_mat)))
  }
  
  stats_m3 <- list()
  stats_m3$psi <- unlist(ests_mat[, 1])
  stats_m3$mu <- unlist(ests_mat[, 2])
  stats_m3$sigma_sq <- unlist(ests_mat[, 3])
  stats_m3$sigma_sq_mu <- unlist(ests_mat[, 4])
  stats_m3$sigma_sq_ep <- unlist(ests_mat[, 5])
  
  return(stats_m3)
}

runMC <- function(simdat) {
  #optim
  mc_optim_output <- list()
  mc_optim_output$m1 <- mc_optim_ll(simdat = simdat$sim_m1, setup = simdat$setup, modtype = 'm1')
  mc_optim_output$m2 <- mc_optim_ll(simdat = simdat$sim_m2, setup = simdat$setup, modtype = 'm2')
  mc_optim_output$m3 <- mc_optim_ll(simdat = simdat$sim_m3, setup = simdat$setup, modtype = 'm3')
  #MLEs
  mc_mle_output <- list()
  mc_mle_output$m1 <- mc_mle_m1(sim_m1 = simdat$sim_m1)
  mc_mle_output$m2 <- mc_mle_m2(sim_m2 = simdat$sim_m2, setup = simdat$setup)
  mc_mle_output$m3 <- mc_mle_m3(sim_m3 = simdat$sim_m3, setup = simdat$setup)
  #U
  mc_u_output <- list()
  mc_u_output$m2 <- mc_u_m2(sim_m2 = simdat$sim_m2)
  mc_u_output$m3 <- mc_u_m3(sim_m3 = simdat$sim_m3, setup = simdat$setup)
  
  simdat$mc_optim_output <- mc_optim_output
  simdat$mc_mle_output <- mc_mle_output
  simdat$mc_u_output <- mc_u_output
  return(simdat)
}


AnalyzeMCs <- function(simdat, opt, ...) {
  
  if (missing(opt)) {opt=0}
  if (opt == 1) {
    simdat$mc_optim_output$m1$results <- NULL
    simdat$mc_optim_output$m2$results <- NULL
    simdat$mc_optim_output$m3$results <- NULL
    simdat$mc_mle_output$m1$results <- NULL
    simdat$mc_mle_output$m2$results <- NULL
    simdat$mc_mle_output$m3$results <- NULL
    simdat$mc_u_output$m2$results <- NULL
    simdat$mc_u_output$m3$results <- NULL
  }
  
  estimators <- c('mc_optim_output', 'mc_mle_output', 'mc_u_output')
  modtypes <- c('m1', 'm2', 'm3')
  datasets <- c('sim_m1', 'sim_m2', 'sim_m3')
  params <- c('mu', 'sigma_sq', 'sigma_sq_mu', 'sigma_sq_ep', 'psi')
  
  for (e in 1:length(estimators)) {
    if (is.null(simdat[[estimators[e]]])) {next}
    mc_output <- simdat[[estimators[e]]]
    
    for (m in 1:length(modtypes)) {
      if (is.null(mc_output[[modtypes[m]]])) {next}
      if ('results' %in% names(mc_output[[modtypes[m]]])) {next}
      
      nparams <- length(mc_output[[modtypes[m]]])
      params_em <- params[1:nparams]
      vals <- as.numeric(do.call(cbind, simdat$setup[params_em]))
      stats <- mc_output[[modtypes[m]]][params_em]
      if (("sigma_sq" %in% params_em) & any(stats$sigma_sq < 0)) {
        print('sigma_sq has a value less than 0')
        i <- which(stats$sigma_sq < 0)
        stats$sigma_sq[which(stats$sigma_sq < 0)] <- 1e-4
      }
      if (("sigma_sq_mu" %in% params_em) & any(stats$sigma_sq_mu < 0)) {
        print('sigma_sq_mu has a value less than 0')
        which(stats$sigma_sq_mu < 0)
        stats$sigma_sq_mu[which(stats$sigma_sq_mu < 0)] <- 1e-4
      }
      if (("sigma_sq_ep" %in% params_em) & any(stats$sigma_sq_ep < 0)) {
        print('sigma_sq_ep has a value less than 0')
        stats$sigma_sq_ep[which(stats$sigma_sq_ep < 0)] <- 1e-4
      }
      
      mat_row0 <- as.matrix(do.call(cbind, stats))
      i_na <- which(apply(mat_row0, 1, function(x) all(is.na(x))))
      if (length(i_na) > 0) {
        mat_row <- mat_row0[-i_na, , drop = FALSE]
      } else {
        mat_row <- mat_row0
      }
      nsims <- nrow(mat_row) 
      
      mn <- as.numeric(colMeans(mat_row))
      std <- as.numeric(apply(mat_row, 2, sd))
      bias <- as.numeric(colMeans(mat_row - matrix(rep(vals, each = nsims), nrow = nsims)))
      rmse <- as.numeric(sqrt(colMeans((mat_row - matrix(rep(vals, each = nsims), nrow = nsims))^2)))
      mae <- as.numeric(colMeans(abs(mat_row - matrix(rep(vals, each = nsims), nrow = nsims))))
      mape <- as.numeric(100 * colMeans(abs(mat_row - matrix(rep(vals, each = nsims), nrow = nsims)) / 
                                          abs(matrix(rep(vals, each = nsims), nrow = nsims))))
      re <- as.numeric(100 * bias / vals)
      
      nparams_cp <- min(nparams, 4)
      i_int <- matrix(0, ncol=nparams_cp)
      indices <- setdiff(1:simdat$setup$n_sims, i_na)
      for (k in indices) {
        df <- tidyr::drop_na(simdat[[datasets[m]]]$dfs[[k]], V, tau) 
        cis <- calcCIs(df = df, results = mat_row0[k,1:nparams_cp])
        for (p in 1:nparams_cp) {
          if ((cis$CI_lower[p] < vals[p]) &  (vals[p] < cis$CI_upper[p])) {i_int[p] <- i_int[p] + 1}
        }
      }
      cp <- as.numeric(i_int / nsims)
      if (nparams>4) {cp <- c(cp, NA)}
      
      simdat[[estimators[e]]][[modtypes[m]]][['results']] <- data.frame(
        estimator = rep(estimators[e], nparams),
        modtype = rep(modtypes[m], nparams),
        param = params_em,
        mean = mn,
        sd = std,
        bias = bias,
        rmse = rmse,
        mae = mae,
        mape = mape,
        re = re,
        cp = cp)
      
    }
  }
  return(simdat)
}