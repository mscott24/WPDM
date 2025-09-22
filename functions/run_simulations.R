rm(list = ls())
setwd('/project/jointage/matt/DM/WPDM/')
source('functions/WPDM.R')
setwd('/project/jointage/matt/DM/WPDM/data/')

####sim 1####
setup <- DefaultSetup()
setup$n_sims <- 1000
setup$n_subjs <- 100
setup$mu <- -1
setup$sigma_sq <- 0.25
setup$sigma_sq_mu <- 0.25
setup$sigma_sq_ep <- 0.25
GenData(setup=setup, sim_id='sim1')

load('sim1_eq.RData')
eq <- runMC(simdat = eq)
eq <- AnalyzeMCs(simdat=eq)
save(eq, file='sim1_eq.RData')
rm(eq)

load('sim1_ue.RData')
ue <- runMC(simdat = ue)
ue <- AnalyzeMCs(simdat=ue)
save(ue, file='sim1_ue.RData')
rm(ue)

# ####sim 2####
setup <- DefaultSetup()
setup$n_sims <- 1000
setup$n_subjs <- 100
setup$mu <- -1
setup$sigma_sq <- 0.25
setup$sigma_sq_mu <- 0.05
setup$sigma_sq_ep <- 0.49
GenData(setup=setup, sim_id='sim2')
#
load('sim2_eq.RData')
eq <- runMC(simdat = eq)
eq <- AnalyzeMCs(simdat=eq)
save(eq, file='sim2_eq.RData')
rm(eq)
load('sim2_ue.RData')
ue <- runMC(simdat = ue)
ue <- AnalyzeMCs(simdat=ue)
save(ue, file='sim2_ue.RData')
rm(ue)

####supplemental simulations####
setwd('/project/jointage/matt/DM/WPDM/data/param_checks')
mu_vals <- c(-1, 0.1, 1)
var_vals <- c(0.05, 0.25, 1)
sim_num <- 3
for (mu in mu_vals) {
  for (var in var_vals) {
    setup <- DefaultSetup()
    setup$n_sims <- 1000
    setup$n_subjs <- 100
    setup$mu <- mu
    setup$sigma_sq <- var
    setup$sigma_sq_mu <- var
    setup$sigma_sq_ep <- var
    ff <- paste0('sim', sim_num)
    GenData(setup = setup, sim_id = ff)
    ff_eq <- paste0(ff, '_eq.RData')
    load(ff_eq)
    eq <- runMC(simdat = eq)
    eq <- AnalyzeMCs(simdat = eq)
    save(eq, file = ff_eq)
    rm(eq)

    ff_ue <- paste0(ff, '_ue.RData')
    load(ff_ue)
    ue <- runMC(simdat = ue)
    ue <- AnalyzeMCs(simdat = ue)
    save(ue, file=ff_ue)
    rm(ue)

    sim_num <- sim_num + 1
  }
}

####number of follow-up visits - figure 1####
setwd('/project/jointage/matt/DM/WPDM/data/nfu_checks/')
n_vals <- 3:10
for (n in n_vals) {
    setup <- DefaultSetup()
    setup$n <- n
    setup$n_sims <- 1000
    setup$n_subjs <- 100
    setup$mu <- -1
    setup$sigma_sq <- 0.25
    setup$sigma_sq_mu <- 0.25
    setup$sigma_sq_ep <- 0.25 
    ff <- paste0('nfu_', n)
    GenData(setup = setup, sim_id = ff, opt = 'eq')
    
    ff_eq <- paste0(ff, '_eq.RData')
    load(ff_eq)
    eq <- runMC(simdat = eq)
    eq <- AnalyzeMCs(simdat = eq)
    save(eq, file = ff_eq)
    rm(eq)
}

###time interval proportion analysis###
rm(list = ls())
setwd('/project/jointage/matt/DM/WPDM/')
source('functions/WPDM.R')
setwd('/project/jointage/matt/DM/WPDM/data/')

setup <- DefaultSetup()
setup$n_sims <- 500
setup$n_subjs <- 100
setup$mu <- -1
setup$sigma_sq <- 0.25
setup$sigma_sq_mu <- 0.25
setup$sigma_sq_ep <- 0.25
setup$n <- 10

t <- 0:10
center <- 5
prop_bunched <- 0.5
sims <- list()
sims$setup <- setup
sims$sim_m3 <- list()
sims$sim_m3$dfs <- list()
params <- c(sims$setup$mu,sims$setup$sigma_sq,sims$setup$sigma_sq_mu,sims$setup$sigma_sq_ep)

set.seed(2024)
n_subj <- sims$setup$n_subjs
prop_bunched <- prop_bunched 
center <- center  
t <- 0:10 
max_fu <- length(t)
deltas <- c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25,  0.30)
idx_list <- replicate(n_subj,
                      sample(seq_along(t), size = ceiling(length(t) * prop_bunched), replace = FALSE),
                      simplify = FALSE)
for (delta in deltas) {
  if (file.exists(paste0('followup_time_checks/delta', delta, '.RData'))) next
  print(delta)
  
  time <- matrix(nrow = length(t), ncol = n_subj)
  for (i in seq_len(n_subj)) {
    idx <- idx_list[[i]]
    tnew <- t
    tnew[idx] <- (1 - delta) * t[idx] + delta * center
    time[, i] <- sort(tnew) - min(tnew)
  }
  
  for (k in 1:sims$setup$n_sims) {
    Y_k_m3 <- matrix(nrow = max_fu, ncol = n_subj)
    for (i in seq_len(n_subj)) {
      Y_k_m3[1:length(t), i] <- RunSim(params = params, t = time[, i])
    }
    sims$sim_m3$dfs[[k]] <- FormatDiffDF(MakeLong(Y_k_m3, time), Yvar = 'Y', idvar = 'Patient', timevar = 't'
    )
  }
  sims$mc_mle_output$m3  <- mc_mle_m3(sim_m3 = sims$sim_m3, setup = sims$setup)
  sims$mc_optim_output$m3 <- mc_optim_ll(simdat = sims$sim_m3, setup = sims$setup, modtype = 'm3')
  sims <- AnalyzeMCs(simdat = sims)
  sims$time <- time
  save(sims, file = paste0('followup_time_checks/delta', delta, '.RData'))
}


