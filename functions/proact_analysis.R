#### PRO-ACT data set generation and analysis

library(purrr)
library(rsample) 
library(rlang) 
library(yardstick) 
library(covBM)

#set your working directory to WPDM 
#setwd('path to WPDM')
source('functions/WPDM.R')

#You need to define these two functions before running covBM or else you will 
#encounter a version error. 
lmeBM <-function(fixed, data, random,	covariance = NULL, method = c("REML", "ML"),
                 control = list(), keep.data = TRUE){
  origCall<-match.call()
  CallList<-as.list(match.call())[-1]
  
  if(missing(covariance)) stop("lmeBM used without specification of covariance structure")
  if(!missing(control) && !is.list(control)) stop("'control' argument is not provided as a list")
  
  index<-which(names(CallList)=="covariance")
  covObject<-eval(CallList[[index]])
  if(!inherits(covObject, c("covBM", "covFracBM", "covIOU"))) stop("lmeBM used without supported covariance structure")
  names(CallList)[index]<-"correlation"		### covariance argument is fed into 'lme' function as correlation argument, hence renaming here
  
  lme_fit <- do.call(lme, CallList, envir = parent.frame())
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

############# DATA CREATION #############
#path to your pro-act data here
setwd('your path to proact_data')

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
frs$Subject_used_Riluzole <- as.factor(frs$Subject_used_Riluzole)

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
         Onset_Duration = time0 - Onset_Delta,
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
                  Onset_Duration>=0 &
                  Onset_Duration<=36 &
                  Age < 80)

############# ANALYSIS #############
#degradation, LME, and covBM model estimates

runPROACT <- function(dat) {
  
  m1_optim <- WPDM(df=dat, modtype = 'm1', estimator = 'optim')
  m1_mle <- WPDM(df=dat, modtype = 'm1', estimator = 'mle')
  m2_optim <- WPDM(df=dat, modtype = 'm2', estimator = 'optim')
  m2_mle <- WPDM(df=dat, modtype = 'm2', estimator = 'mle')
  m2_u <- WPDM(df=dat, modtype = 'm2', estimator = 'u')
  m3_optim <- WPDM(df=dat, modtype = 'm3', estimator = 'optim')
  m3_mle <- WPDM(df=dat, modtype = 'm3', estimator = 'mle')
  m3_u <- WPDM(df=dat, modtype = 'm3', estimator = 'u')
  
  # unadjusted LME
  l_lme <- nlme::lme(fixed = Y ~ t,
                     data=dat,
                     random=~t|Patient,
                     method="ML",
                     control=lmeControl(maxIter = 100000, tolerance = 1e-3))
  
  # unadjusted covBM model
  l_covbm_unadjusted <- lmeBM(fixed = Y ~ t,
                              data=dat,
                              random=~t|Patient,
                              covariance=covBM(form=~t|Patient),
                              method="ML",
                              control=lmeControl(maxIter = 100000, tolerance = 1e-3))
  
  
  stats <- list()
  stats$m1_optim <- m1_optim
  stats$m1_mle <- m1_mle
  stats$m2_optim <- m2_optim
  stats$m2_mle <- m2_mle
  stats$m2_u <- m2_u
  stats$m3_optim <- m3_optim
  stats$m3_mle <- m3_mle
  stats$m3_u <- m3_u
  stats$l_lme <- l_lme
  stats$l_covbm_unadjusted <- l_covbm_unadjusted
  return(stats)
  
}

all <- final
alive <- subset(final, Subject_Died=='No')
died <- subset(final, Subject_Died=='Yes')

stats_all <- runPROACT(all)
stats_alive <- runPROACT(alive)
stats_died <- runPROACT(died)

# adjusted covBM model
all_wcovs <- subset(all, !is.na(Subject_used_Riluzole) & !is.na(Sex))
age_mean <- mean(all_wcovs$Age, na.rm = TRUE)
all_wcovs$Age_c <- all_wcovs$Age - age_mean
l_covbm_adjusted <- lmeBM(fixed = Y ~ t*Subject_Died + Age_c*t + Sex*t + Subject_used_Riluzole*t,
                          data=all_wcovs,
                          random=~t|Patient,
                          covariance=covBM(form=~t|Patient),
                          method="ML",
                          control=lmeControl(maxIter = 100000, tolerance = 1e-3))
summary(l_covbm_adjusted)
stop()
save(stats_all, stats_alive, stats_died, l_covbm_adjusted, file = "proact_results.RData")