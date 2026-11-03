#### Code to run the primary simulations

#Reproducibility note: seeds are set to 2024 before each data-generating step 
#in the 'GenData' function of the 'simulations.R' file.

#Set your path in the WPDM folder here
#setwd('path to WPDM')
source('functions/WPDM.R')
setwd('data')

####sim 1####
setup <- DefaultSetup()
setup$home_dir <- getwd()
setup$n_sims <- 1000
setup$n_subjs <- 100
setup$mu <- -1
setup$sigma_sq <- 0.25
setup$sigma_sq_mu <- 0.25
setup$sigma_sq_ep <- 0.25
setup$psi <- setup$sigma_sq_ep / setup$sigma_sq
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

#####sim 2####
setup <- DefaultSetup()
setup$home_dir <- getwd()
setup$n_sims <- 1000
setup$n_subjs <- 100
setup$mu <- -1
setup$sigma_sq <- 0.25
setup$sigma_sq_mu <- 0.05
setup$sigma_sq_ep <- 0.49
setup$psi <- setup$sigma_sq_ep / setup$sigma_sq
GenData(setup=setup, sim_id='sim2')

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

####varying number of follow-up visits####
setwd('~/Desktop/Dissertation/WPDM/data/nfu_checks/')
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
    setup$psi <- setup$sigma_sq_ep / setup$sigma_sq
    ff <- paste0('nfu_', n)
    GenData(setup = setup, sim_id = ff, opt = 'eq')
    
    ff_eq <- paste0(ff, '_eq.RData')
    load(ff_eq)
    eq <- runMC(simdat = eq)
    eq <- AnalyzeMCs(simdat = eq)
    save(eq, file = ff_eq)
    rm(eq)
}