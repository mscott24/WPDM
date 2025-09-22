####robustness
rm(list = ls())
setwd('/project/jointage/matt/DM/WPDM/')
source('functions/WPDM.R')
setwd('data/')
library(xtable)

genTables <- function(simdat) {
  #optim
  m1 <- data.frame(t(simdat$mc_optim_output$m1$results))
  m1$params <- rownames(m1)
  m1$model <- 'm1'
  m1 <- m1[-(1:3), ]
  m2 <- data.frame(t(simdat$mc_optim_output$m2$results))
  m2$params <- rownames(m2)
  m2$model <- 'm2'
  m2 <- m2[-(1:3), ]
  m3 <- data.frame(t(simdat$mc_optim_output$m3$results))
  m3$params <- rownames(m3)
  m3$model <- 'm3'
  m3 <- m3[-(1:3), ]
  opt <- bind_rows(m1, m2, m3)
  rownames(opt) <- NULL
  opt <- opt %>% select(model, params, X1, X2, X3, X4)
  colnames(opt) <- c('model', 'params', 'mu', 'sigma_sq', 'sigma_sq_mu', 'sigma_sq_ep')
  
  #mle
  m1 <- data.frame(t(simdat$mc_mle_output$m1$results))
  m1$params <- rownames(m1)
  m1$model <- 'm1'
  m1 <- m1[-(1:3), ]
  m2 <- data.frame(t(simdat$mc_mle_output$m2$results))
  m2$params <- rownames(m2)
  m2$model <- 'm2'
  m2 <- m2[-(1:3), ]
  m3 <- data.frame(t(simdat$mc_mle_output$m3$results))
  m3$params <- rownames(m3)
  m3$model <- 'm3'
  m3 <- m3[-(1:3), ]
  mle <- bind_rows(m1, m2, m3)
  rownames(mle) <- NULL
  mle <- mle %>% select(model, params, X1, X2, X3, X4, X5)
  colnames(mle) <- c('model', 'params', 'mu', 'sigma_sq', 'sigma_sq_mu', 'sigma_sq_ep', 'psi')
  
  #u
  m2 <- data.frame(t(simdat$mc_u_output$m2$results))
  m2$params <- rownames(m2)
  m2$model <- 'm2'
  m2 <- m2[-(1:3), ]
  m3 <- data.frame(t(simdat$mc_u_output$m3$results))
  m3$params <- rownames(m3)
  m3$model <- 'm3'
  m3 <- m3[-(1:3), ]
  u <- bind_rows(m2, m3)
  rownames(u) <- NULL
  u <- u %>% select(model, params, X1, X2, X3, X4, X5)
  colnames(u) <- c('model', 'params', 'mu', 'sigma_sq', 'sigma_sq_mu', 'sigma_sq_ep', 'psi')
  
  opt <- opt %>% rename_with(~ paste0(., ".optim"), -c(model, params))
  mle <- mle %>% rename_with(~ paste0(., ".mle"),   -c(model, params))
  u <- u %>% rename_with(~ paste0(., ".u"),     -c(model, params))
  
  result <- opt %>%
    left_join(mle, by = c("model", "params")) %>%
    left_join(u, by = c("model", "params"))
  return(result)
}
#table c1
load('sim1_eq.RData')
tbl_c1_a <- genTables(simdat=eq)
tbl_c1_a$simulation <- 'sim 1'
load('sim2_eq.RData')
tbl_c1_b <- genTables(simdat=eq)
tbl_c1_b$simulation <- 'sim 2'
out <- rbind(tbl_c1_a, tbl_c1_b)
xtable(out)

#table c2
load('sim1_ue.RData')
tbl_c2_a <- genTables(simdat=ue)
tbl_c2_a$simulation <- 'sim 1'
load('sim2_ue.RData')
tbl_c2_b <- genTables(simdat=ue)
tbl_c2_b$simulation <- 'sim 2'
out <- rbind(tbl_c2_a, tbl_c2_b)
xtable(out)

#supplemental simulations
#'sim3': mu = -1.00, vars = 0.05
load('param_checks/sim3_eq.RData')
View(genTables(simdat=eq))
print(eq$setup)
#'sim4': mu = -1.00, vars = 0.25
load('param_checks/sim4_eq.RData')
View(genTables(simdat=eq))
print(eq$setup)
#'sim5': mu = -1.00, vars = 1.00
load('param_checks/sim5_eq.RData')
View(genTables(simdat=eq))
print(eq$setup)
#'sim6': mu = 0.01, vars = 0.05
load('param_checks/sim6_eq.RData')
View(genTables(simdat=eq))
print(eq$setup)
#'sim7': mu = 0.01, vars = 0.25
load('param_checks/sim7_eq.RData')
View(genTables(simdat=eq))
print(eq$setup)
#'sim8': mu = 0.01, vars = 1.00
load('param_checks/sim8_eq.RData')
View(genTables(simdat=eq))
print(eq$setup)
#'sim9': mu = -1.00, vars = 0.05
load('param_checks/sim9_eq.RData')
View(genTables(simdat=eq))
print(eq$setup)
#'sim10': mu = -1.00, vars = 0.25
load('param_checks/sim10_eq.RData')
View(genTables(simdat=eq))
print(eq$setup)
#'sim11': mu = -1.00, vars = 1.00
load('param_checks/sim11_eq.RData')
View(genTables(simdat=eq))
print(eq$setup)