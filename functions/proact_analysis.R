rm(list = ls())
setwd('/project/jointage/matt/DM/WPDM/')
source('functions/WPDM.R')
setwd('data/')

library(purrr)
library(rsample) 
library(rlang) 
library(yardstick) 
library(covBM)

############# DATA CREATION #############
setwd('/project/jointage/matt/DM/paper/proact_data/')
frs <- read.csv('PROACT_ALSFRS.csv')
dth <- read.csv('PROACT_DEATHDATA.csv')
demog <- read.csv('PROACT_DEMOGRAPHICS.csv')
trt <- read.csv('PROACT_TREATMENT.csv')
rilu <- read.csv('PROACT_RILUZOLE.csv')
hist <- read.csv('PROACT_ALSHISTORY.csv')

# remove duplicate row for death
dth <- subset(dth, subject_id!=442984)

# remove duplicate onset time
hist <- hist %>% 
  dplyr::select(subject_id, Onset_Delta, Diagnosis_Delta) %>% 
  distinct() %>% 
  filter(!is.na(Onset_Delta) | !is.na(Diagnosis_Delta))

# join data.frames to FRS as base
frs <- frs %>%
  arrange(subject_id, ALSFRS_Delta) %>%
  left_join(demog, by = "subject_id") %>%
  left_join(dth, by = "subject_id", ) %>%
  left_join(trt, by = "subject_id") %>%
  left_join(rilu, by = "subject_id") %>%
  left_join(hist, by = "subject_id") 

frs$subject_id <- as.factor(frs$subject_id)
frs$ALSFRS_Delta <- as.numeric(frs$ALSFRS_Delta)
frs$ALSFRS_Total <- as.numeric(frs$ALSFRS_Total)
frs$ALSFRS_R_Total <- as.numeric(frs$ALSFRS_R_Total)
frs$Onset_Delta <- as.numeric(frs$Onset_Delta)
frs$Sex <- as.factor(frs$Sex)
frs$Ethnicity <- as.factor(frs$Ethnicity)
frs$Subject_used_Riluzole <- as.factor(frs$Subject_used_Riluzole)

# define ethnicity variable
i <- which(frs$Ethnicity=='')
frs$Ethnicity[i] <- 'Unknown'
frs$Ethnicity <- droplevels(frs$Ethnicity)

# create three-category death variable
frs$Subject_Died[is.na(frs$Subject_Died)] <- 'Unknown'
frs$Subject_Died <- as.factor(frs$Subject_Died)

# convert diagnostic time from days to months (age is still in years)
frs$ALSFRS_Delta <- as.numeric(frs$ALSFRS_Delta / 30.44)
frs$Onset_Delta <- as.numeric(frs$Onset_Delta / 30.44)
frs$Diagnosis_Delta <- as.numeric(frs$Diagnosis_Delta / 30.44)

# remove NA FRS and time entries
frs <- frs %>% filter(!is.na(ALSFRS_R_Total) & !is.na(ALSFRS_Delta))

# remove duplicate rows
frs <- frs %>% distinct()

# define additional measure
frs <- frs %>% 
  arrange(subject_id, ALSFRS_Delta) %>%
  group_by(subject_id) %>% 
  mutate(delta_frs = last(ALSFRS_R_Total) - first(ALSFRS_R_Total), 
         delta_time = last(ALSFRS_Delta) - first(ALSFRS_Delta),
         time0 = first(ALSFRS_Delta),
         timeLast = last(ALSFRS_Delta),
         time = ALSFRS_Delta - first(ALSFRS_Delta),
         frs_bl = first(ALSFRS_R_Total),
         slope = (last(ALSFRS_R_Total) - first(ALSFRS_R_Total))/(last(ALSFRS_Delta) - first(ALSFRS_Delta)),
         n_fu = n())

# only in patients with 3 or more follow-up visits
frs <- subset(frs, n_fu>3)

# remove instances where there are multiple FRS values on the same day 
# (take first one) and add in difference measures
frs <- AddDiff(df = frs, Yvar ='ALSFRS_R_Total', timevar = 'time', idvar = 'subject_id')
frs <- frs %>% filter(tau!=0 | is.na(tau))

final <- subset(frs, Subject_Died!='Unknown' &
                  Study_Arm=='Placebo' & 
                  Onset_Delta>=-36 &
                  Age < 80)

############# ANALYSIS #############
# m1_optim <- WPDM(df=final, modtype = 'm1', estimator = 'optim')
# m1_mle <- WPDM(df=final, modtype = 'm1', estimator = 'mle')
# m2_optim <- WPDM(df=final, modtype = 'm2', estimator = 'optim')
# m2_mle <- WPDM(df=final, modtype = 'm2', estimator = 'mle')
# m2_u <- WPDM(df=final, modtype = 'm2', estimator = 'u')
# m3_optim <- WPDM(df=final, modtype = 'm2', estimator = 'optim')
# m3_mle <- WPDM(df=final, modtype = 'm2', estimator = 'mle')
# m3_u <- WPDM(df=final, modtype = 'm2', estimator = 'u')
# l <- lmeBM(fixed = Y ~ t,
#            data=final,
#            random=~t|Patient,
#            covariance=covBM(form=~t|Patient),
#            method="ML", 
#            control=lmeControl(maxIter = 100000, tolerance = 1e-3))
# final_allcovs <- subset(final, !is.na(Subject_used_Riluzole) & !is.null(Sex))
# l_adj <- lmeBM(fixed = Y ~ t*Subject_Died + Age*t + Sex*t + Subject_used_Riluzole*t,
#                data=final_allcovs,
#                random=~t|Patient,
#                covariance=covBM(form=~t|Patient),
#                method="ML", 
#                control=lmeControl(maxIter = 100000, tolerance = 1e-3))

###run cross validation
setwd('/project/jointage/matt/DM/WPDM/data')

# covBM package source functions 
lmeBM <-function(fixed,
                 data,					
                 random,					
                 covariance = NULL,
                 method = c("REML", "ML"),
                 control = list(),
                 keep.data = TRUE){
  origCall<-match.call()
  CallList<-as.list(match.call())[-1]
  
  if(missing(covariance)) stop("lmeBM used without specification of covariance structure")
  if(!missing(control) && !is.list(control)) stop("'control' argument is not provided as a list")
  
  index<-which(names(CallList)=="covariance")
  covObject<-eval(CallList[[index]])
  if(!inherits(covObject, c("covBM", "covFracBM", "covIOU"))) stop("lmeBM used without supported covariance structure")
  names(CallList)[index]<-"correlation"		### covariance argument is fed into 'lme' function as correlation argument, hence renaming here
  
  lme_fit<-do.call(lme, CallList)
  if(!is.null(lme_fit$modelStruct$corStruct) && inherits(lme_fit$modelStruct$corStruct, c("covBM", "covFracBM", "covIOU"))){
    attr(lme_fit$modelStruct$corStruct, "sigma") <-	lme_fit$sigma
  }
  
  if(!is.character(lme_fit$apVar)){lme_fit$apVar<-refactor_apVar(lme_fit$apVar, class(lme_fit$modelStruct$corStruct)[1])}
  
  lme_fit$call<-origCall
  
  return(lme_fit)
}
refactor_apVar<-function(apVar, cov_class){
  
  npars<-dim(apVar)[1]
  
  if(cov_class=="covBM"){
    #browser()
    for(i in 1:(npars-2)){
      apVar[npars-1, i]<- 2*apVar[npars,i] + apVar[npars-1,i]
    }
    apVar[npars-1, npars-1]<- 4*apVar[npars,npars] + 4*apVar[npars,npars-1] + apVar[npars-1,npars-1]
    apVar[npars-1, npars]<- 2*apVar[npars,npars] + apVar[npars-1,npars]
    for(i in 1:npars){
      if(i!=npars-1) {apVar[i, npars-1]<-apVar[npars-1, i]}
    }
    attr(apVar,"Pars")["corStruct"]<- as.numeric(attr(apVar,"Pars")["corStruct"] + 2*attr(apVar,"Pars")["lSigma"])
  }
  
  if(cov_class=="covFracBM" || cov_class=="covIOU"){
    for(i in 1:(npars-3)){
      apVar[npars-2, i]<- 2*apVar[npars,i] + apVar[npars-2,i]
    }
    apVar[npars-2, npars-2]<- 4*apVar[npars,npars] + 4*apVar[npars,npars-2] + apVar[npars-2,npars-2]
    apVar[npars-2, npars-1]<- 2*apVar[npars,npars-1] + apVar[npars-2,npars-1]
    apVar[npars-2, npars]<- 2*apVar[npars,npars] + apVar[npars-2,npars]
    for(i in 1:npars){
      if(i!=npars-2) {apVar[i, npars-2]<-apVar[npars-2, i]}
    }
    attr(apVar,"Pars")["corStruct1"]<- as.numeric(attr(apVar,"Pars")["corStruct1"] + 2*attr(apVar,"Pars")["lSigma"])
  }		
  
  return(apVar)
}

predictWPDM <- function(y0, mu, t) {
  return(y0 + mu*t)
}

calcPreds <- function(pred, y, label) {
  mae  <- mean(abs(y - pred))
  rmse <- sqrt(mean((y - pred)^2))
  mase <- mae / mean(abs(y - mean(y)))
  return(data.frame(mae, rmse, mase, label))
}

####Cross validation####
seed <- 2024
k <- 5
df <- final

set.seed(seed)

unique_patients <- unique(df$Patient)
patient_folds <- sample(rep(1:k, length.out = length(unique_patients)))
names(patient_folds) <- unique_patients
folds <- patient_folds[as.character(df$Patient)]

cv <- list()
cv$m1_optim_pred <- matrix(nrow=nrow(df), ncol=1)
cv$m2_optim_pred <- matrix(nrow=nrow(df), ncol=1)
cv$m3_optim_pred <- matrix(nrow=nrow(df), ncol=1)
cv$m1_mle_pred <- matrix(nrow=nrow(df), ncol=1)
cv$m2_mle_pred <- matrix(nrow=nrow(df), ncol=1)
cv$m3_mle_pred <- matrix(nrow=nrow(df), ncol=1)
cv$m2_u_pred <- matrix(nrow=nrow(df), ncol=1)
cv$m3_u_pred <- matrix(nrow=nrow(df), ncol=1)
cv$m1_optim_ests <- matrix(nrow=k, ncol=2)
cv$m2_optim_ests <- matrix(nrow=k, ncol=3)
cv$m3_optim_ests <- matrix(nrow=k, ncol=4)
cv$m1_mle_ests <- matrix(nrow=k, ncol=2)
cv$m2_mle_ests <- matrix(nrow=k, ncol=3)
cv$m3_mle_ests <- matrix(nrow=k, ncol=4)
cv$m2_u_ests <- matrix(nrow=k, ncol=3)
cv$m3_u_ests <- matrix(nrow=k, ncol=4)
cv$lme_pred <- matrix(nrow=nrow(df), ncol=1)
cv$lme_ests <- matrix(nrow=k, ncol=6)
cv$lme_bl_pred <- matrix(nrow=nrow(df), ncol=1)
cv$lme_bl_ests <- matrix(nrow=k, ncol=6)
cv$bm_pred <- matrix(nrow=nrow(df), ncol=1)
cv$bm_ests <- matrix(nrow=k, ncol=7)
cv$bm_bl_pred <- matrix(nrow=nrow(df), ncol=1)
cv$bm_bl_ests <- matrix(nrow=k, ncol=7)
for (i in seq_len(k)) {
  cat("---- Fold", i, "----\n")
  test_idx <- which(folds == i)
  train_idx <- which(folds != i)
  train <- df[train_idx, ]
  test  <- df[test_idx, ]
  
  cv$m1_optim_ests[i,] <- WPDM(df = train, modtype = 'm1', estimator = 'optim')$stats$Estimate
  cv$m1_optim_pred[test_idx] <- predictWPDM(y0=test$frs_bl, mu=cv$m1_optim_ests[i,1], t=test$t)    
  cv$m2_optim_ests[i,] <- WPDM(df = train, modtype = 'm2', estimator = 'optim')$stats$Estimate
  cv$m2_optim_pred[test_idx] <- predictWPDM(y0=test$frs_bl, mu=cv$m2_optim_ests[i,1], t=test$t)    
  cv$m3_optim_ests[i,] <- WPDM(df = train, modtype = 'm3', estimator = 'optim')$stats$Estimate
  cv$m3_optim_pred[test_idx] <- predictWPDM(y0=test$frs_bl, mu=cv$m3_optim_ests[i,1], t=test$t)     
  
  cv$m1_mle_ests[i,] <- WPDM(df = train, modtype = 'm1', estimator = 'mle')$stats$Estimate
  cv$m1_mle_pred[test_idx] <- predictWPDM(y0=test$frs_bl, mu=cv$m1_mle_ests[i,1], t=test$t)    
  cv$m2_mle_ests[i,] <- WPDM(df = train, modtype = 'm2', estimator = 'mle')$stats$Estimate
  cv$m2_mle_pred[test_idx] <- predictWPDM(y0=test$frs_bl, mu=cv$m2_mle_ests[i,1], t=test$t)    
  cv$m3_mle_ests[i,] <- WPDM(df = train, modtype = 'm3', estimator = 'mle')$stats$Estimate
  cv$m3_mle_pred[test_idx] <- predictWPDM(y0=test$frs_bl, mu=cv$m3_mle_ests[i,1], t=test$t)   
  
  cv$m2_u_ests[i,] <- WPDM(df = train, modtype = 'm2', estimator = 'u')$stats$Estimate
  cv$m2_u_pred[test_idx] <- predictWPDM(y0=test$frs_bl, mu=cv$m2_u_ests[i,1], t=test$t)    
  cv$m3_u_ests[i,] <- WPDM(df = train, modtype = 'm3', estimator = 'u')$stats$Estimate
  cv$m3_u_pred[test_idx] <- predictWPDM(y0=test$frs_bl, mu=cv$m3_u_ests[i,1], t=test$t)   
  
  fit_bm <- lmeBM(
    fixed = Y ~ t,
    data = train,
    random = ~ t | Patient,
    covariance = covBM(form = ~ t | Patient),
    method = "ML",
    control = lmeControl(maxIter = 1e5, tolerance = 1e-3))
  cv$bm_pred[test_idx] <- as.matrix(predict(fit_bm, newdata = test, level = 0))
  
  b0 <- fixef(fit_bm)["(Intercept)"]
  b1 <- fixef(fit_bm)["t"]
  vc  <- VarCorr(fit_bm)
  var_int <- as.numeric(vc["(Intercept)", "Variance"])
  var_slope <- as.numeric(vc["t", "Variance"])
  cor_slope_int <- as.numeric(vc["t", "Corr"])
  resid_var <- as.numeric(vc["Residual", "Variance"])
  var_bm <- intervals(fit_bm)$corStruct[2]
  cv$bm_ests[i,] <- as.numeric(c(b0, b1, var_int, var_slope, cor_slope_int, resid_var, var_bm))
  
  fit_lme <- lme(
    fixed = Y ~ t,
    data = train,
    random = ~ t | Patient,
    method = "ML",
    control = lmeControl(opt = "optim", maxIter = 1e5, msMaxIter = 1e5, tolerance = 1e-6))
  cv$lme_pred[test_idx] <- as.matrix(predict(fit_lme, newdata = test, level = 0))
  
  b0 <- fixef(fit_lme)["(Intercept)"]
  b1 <- fixef(fit_lme)["t"]
  vc  <- VarCorr(fit_lme)
  var_int <- as.numeric(vc["(Intercept)", "Variance"])
  var_slope <- as.numeric(vc["t", "Variance"])
  cor_slope_int <- as.numeric(vc["t", "Corr"])
  resid_var <- as.numeric(vc["Residual", "Variance"])
  cv$lme_ests[i,] <- as.numeric(c(b0, b1, var_int, var_slope, cor_slope_int, resid_var))
  
  #lmes with bl Y as covariate
  train_star <- subset(train, t!=0)
  fit_lme_bl <- lme(
    fixed = Y ~ t + frs_bl,
    data = train,
    random = ~ t | Patient,
    method = "ML",
    control = lmeControl(opt = "optim", maxIter = 1e5, msMaxIter = 1e5, tolerance = 1e-6))
  cv$lme_bl_pred[test_idx] <- as.matrix(predict(fit_lme_bl, newdata = test, level = 0))
  
  b0 <- fixef(fit_lme_bl)["(Intercept)"]
  b1 <- fixef(fit_lme_bl)["t"]
  vc  <- VarCorr(fit_lme_bl)
  var_int <- as.numeric(vc["(Intercept)", "Variance"])
  var_slope <- as.numeric(vc["t", "Variance"])
  cor_slope_int <- as.numeric(vc["t", "Corr"])
  resid_var <- as.numeric(vc["Residual", "Variance"])
  cv$lme_bl_ests[i,] <- as.numeric(c(b0, b1, var_int, var_slope, cor_slope_int, resid_var))
}

save(cv, file ='cv_results.RData')
stop()

load('cv_results.RData')
vals_m1_optim <- calcPreds(pred=as.numeric(cv$m1_optim_pred), y=final$Y, label='M1 optim')
vals_m2_optim <- calcPreds(pred=as.numeric(cv$m2_optim_pred), y=final$Y, label='M2 optim')
vals_m3_optim <- calcPreds(pred=as.numeric(cv$m3_optim_pred), y=final$Y, label='M3 optim')
vals_m1_mle <- calcPreds(pred=as.numeric(cv$m1_mle_pred), y=final$Y, label='M1 mle')
vals_m2_mle <- calcPreds(pred=as.numeric(cv$m2_mle_pred), y=final$Y, label='M2 mle')
vals_m3_mle <- calcPreds(pred=as.numeric(cv$m3_mle_pred), y=final$Y, label='M3 mle')
vals_m2_u <- calcPreds(pred=as.numeric(cv$m2_u_pred), y=final$Y, label='M2 u')
vals_m3_u <- calcPreds(pred=as.numeric(cv$m3_u_pred), y=final$Y, label='M3 u')
vals_lme <- calcPreds(pred=as.numeric(cv$lme_pred), y=final$Y, label='LME')
vals_bm <- calcPreds(pred=as.numeric(cv$bm_pred), y=final$Y, label='LME+BM')

summary_table <- rbind(vals_m1_optim, vals_m2_optim, vals_m3_optim,
                       vals_m1_mle, vals_m2_mle, vals_m3_mle,
                       vals_m2_u, vals_m3_u, 
                       vals_lme, vals_bm)
print(summary_table, digits = 3)