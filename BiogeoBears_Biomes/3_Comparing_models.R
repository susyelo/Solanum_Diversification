### Libraries ----
library(optimx)  # (either 2012 or 2013 version, as of January 2014)
library(FD)        # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; prob. better than library(parallel))
library(parallel)
library(minqa)
library(BioGeoBEARS)
library(plyr)
library(magrittr)

##cladoRcpp.R needs to be source first
source("./Functions/cladoRcpp.R")
sourceall("Functions/")
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)
# slight speedup hopefully

#####################
# Functions ----
#####################

load_obj <- function(f)
{
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}


resDEC_M0<-load_obj("./Output/Solanum_DEC_M0_unconstrained_v1.Rdata")
resDECj_M0<-load_obj("./Output/Solanum_DEC+J_M0_unconstrained_v1.Rdata")
resDEC_M1<-load_obj("./Output/Solanum_DEC_M1_2_v1.Rdata")
resDECj_M1<-load_obj("./Output/Solanum_DEC+J_M1_2_v1.Rdata")

resDIVALIKE_M0<-load_obj("./Output/Solanum_DIVALIKE_M0_unconstrained_v1.Rdata")
resDIVALIKEj_M0<-load_obj("./Output/Solanum_DIVALIKE+J_M0_unconstrained_v1.Rdata")
resDIVALIKE_M1<-load_obj("./Output/Solanum_DIVALIKE_M1_2_v1.Rdata")
resDIVALIKEj_M1<-load_obj("./Output/Solanum_DIVALIKE+J_M1_2_v1.Rdata")

resBAYAREALIKE_M0<-load_obj("./Output/Solanum_BAYEAREALIKE_M0_unconstrained_v1.Rdata")
resBAYAREALIKEj_M0<-load_obj("./Output/Solanum_BAYAREALIKE+J_M0_unconstrained_v1.Rdata")
resBAYAREALIKE_M1<-load_obj("./Output/Solanum_BAYEAREALIKE_M1_2_v1.Rdata")
resBAYAREALIKEj_M1<-load_obj("./Output/Solanum_BAYAREALIKE+J_M1_2_v1.Rdata")

################################
#### Creating AIC table ----
################################
models<-list(resDEC_M0,resDECj_M0,resDEC_M1,resDECj_M1,
             resDIVALIKE_M0,resDIVALIKEj_M0,resDIVALIKE_M1,resDIVALIKEj_M1,
             resBAYAREALIKE_M0,resBAYAREALIKEj_M0,resBAYAREALIKE_M1,resBAYAREALIKEj_M1)

Models_table<-
  models %>% 
  lapply(. %>% extract_params_from_BioGeoBEARS_results_object)%>%
  ldply

Models_table$AIC<-2*Models_table$numparams-2*Models_table$LnL
Models_table$deltaAIC<-Models_table$AIC-min(Models_table$AIC)
Models_table<-round(Models_table,3)
Models_table$Models<-c("DEC_M0", "DEC+J_M0", "DEC_M1", "DEC+J_M1",
                       "DIVALIKE_M0", "DIVALIKE+J_M0", "DIVALIKE_M1", "DIVALIKE+J_M1",
                       "BAYAREALIKE_M0", "BAYAREALIKE+J_M0","BAYAREALIKE_M1", "BAYAREALIKE+J_M1")

Models_table<-Models_table[order(Models_table$AIC),]
write.csv(Models_table, "./Output/Models_AIC_table.csv")

