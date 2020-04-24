#### Libraries --
#install.packages("BioGeoBEARS", dependencies=TRUE, repos="http://cran.rstudio.com")
#install.packages("snow",dependencies=TRUE)
library(optimx)  # (either 2012 or 2013 version, as of January 2014)
library(FD)        # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; prob. better than library(parallel))
library(parallel)
library(minqa)
library(BioGeoBEARS)
library(plyr)
library(magrittr)

########################################################
# TO GET THE OPTIMX/OPTIM FIX, AND THE UPPASS FIX, 
# SOURCE THE REVISED FUNCTIONS WITH THESE COMMANDS
#
# CRUCIAL CRUCIAL CRUCIAL: 
# YOU HAVE TO RUN THE SOURCE COMMANDS AFTER 
# *EVERY TIME* YOU DO library(BioGeoBEARS). THE CHANGES ARE NOT "PERMANENT", 
# THEY HAVE TO BE MADE EACH TIME.  IF YOU ARE GOING TO BE OFFLINE, 
# YOU CAN DOWNLOAD EACH .R FILE TO YOUR HARD DRIVE AND REFER THE source()
# COMMANDS TO THE FULL PATH AND FILENAME OF EACH FILE ON YOUR
# LOCAL SYSTEM INSTEAD.
########################################################

##cladoRcpp.R needs to be source first
source("./Functions/cladoRcpp.R")
sourceall("Functions/")
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)
# slight speedup hopefully


#########################################################
# General data ----
########################################################
tr<-read.tree("./Data/Solanum_biomes_tree.newick")
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn="./Data/biomes_presence_matrix.data")
#tipranges

# Set the maximum number of areas any species may occupy; this cannot be larger 
# than the number of areas you set up, but it can be smaller.
max_range_size = 4



#######################################################
## functions ---
########################################################

## Preparing BGB run object
BGB_object_default<-function(tree_path, geo_data_path){
  BGB_object = define_BioGeoBEARS_run()
  BGB_object$trfn =tree_path
  BGB_object$geogfn = geo_data_path
  BGB_object$max_range_size = max_range_size
  BGB_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BGB_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  BGB_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BGB_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
  BGB_object$num_cores_to_use = 4
  BGB_object$force_sparse=FALSE # sparse=FALSE causes pathology & isn't much faster at this scale
  BGB_object = readfiles_BioGeoBEARS_run(BGB_object)
  BGB_object$return_condlikes_table = TRUE
  BGB_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BGB_object$calc_ancprobs = TRUE    # get ancestral states from optim run
  check_BioGeoBEARS_run(BGB_object)
  return(BGB_object)
}

## BGB Run function 
BGB_run_default<-function(runslow = TRUE,name_file=" ", BGB_object){
  resfn = name_file
  if (runslow)
  {
    res = bears_optim_run(BGB_object)
    res    
    save(res, file=resfn)
  
  } else {
    # Loads to "res"
    load(resfn)
  }
  return(res)
}


## BGB plot function
plot_BGB_two_models<-function(Title= "", BGB_result){
  analysis_titletxt =Title
  # Setup
  results_object = BGB_result
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  
  # States
  res = plot_BioGeoBEARS_results(results_object, 
                                  analysis_titletxt, 
                                  addl_params=list("j"), 
                                  plotwhat="text", 
                                  label.offset=0.45, 
                                  tipcex=0.7, 
                                  statecex=0.7, 
                                  splitcex=0.6, titlecex=0.8, plotsplits=TRUE,
                                  cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  
  # Pie chart
  plot_BioGeoBEARS_results(results_object, 
                           analysis_titletxt,
                           addl_params=list("j"), 
                           plotwhat="pie", 
                           label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, 
                           plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  
return(res)
}

#######################################################
# Run DEC ----
#######################################################
BGB_object<-BGB_object_default(tree_path = "./Data/Solanum_biomes_tree.newick",
                               geo_data_path = "./Data/biomes_presence_matrix.data")

#BGB_object$areas_adjacency_fn = "./Data/Areas_adjacency.txt"
BGB_object<-readfiles_BioGeoBEARS_run(BGB_object)
check_BioGeoBEARS_run(BGB_object)



res_DECM0<-BGB_run_default(runslow = TRUE,
                       name_file="./Output/Solanum_DEC_M0_unconstrained_v1.Rdata", 
                       BGB_object=BGB_object)


#######################################################
# Run DEC+J ----
#######################################################
BGB_object<-BGB_object_default(tree_path = "./Data/Solanum_biomes_tree.newick",
                               geo_data_path = "./Data/biomes_presence_matrix.data")

#BGB_object$areas_adjacency_fn = "./Data/Areas_adjacency.txt"
BGB_object<-readfiles_BioGeoBEARS_run(BGB_object)




# Set up DEC+J model
#********************
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = res_DECM0$outputs@params_table["d","est"]
estart = res_DECM0$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BGB_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BGB_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BGB_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BGB_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add j as a free parameters
BGB_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BGB_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BGB_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

check_BioGeoBEARS_run(BGB_object)


res_DECjM0<-BGB_run_default(runslow = TRUE,
                            name_file="./Output/Solanum_DEC+J_M0_unconstrained_v1.Rdata", 
                            BGB_object=BGB_object)


#######################################################
# Run DIVALIKE ----
#######################################################
BGB_object<-BGB_object_default(tree_path = "./Data/Solanum_corrected.newick",
                               geo_data_path = "./Data/Solanum_geo_2.data")

#BGB_object$areas_adjacency_fn = "./Data/Areas_adjacency.txt"
BGB_object<-readfiles_BioGeoBEARS_run(BGB_object)

# Set up DIVALIKE model
# Remove subset-sympatry
BGB_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BGB_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BGB_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BGB_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BGB_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BGB_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BGB_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BGB_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BGB_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BGB_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5


check_BioGeoBEARS_run(BGB_object)

resDIVALIKE_M0<-BGB_run_default(runslow = TRUE,
                                name_file="./Output/Solanum_DIVALIKE_M0_unconstrained_v1.Rdata", 
                                BGB_object=BGB_object)


#######################################################
# Run DIVALIKE+J ----
#######################################################
BGB_object<-BGB_object_default(tree_path = "./Data/Solanum_biomes_tree.newick",
                               geo_data_path = "./Data/biomes_presence_matrix.data")

#BGB_object$areas_adjacency_fn = "./Data/Areas_adjacency.txt"
BGB_object<-readfiles_BioGeoBEARS_run(BGB_object)

# Set up DIVALIKE+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDIVALIKE_M0$outputs@params_table["d","est"]
estart = resDIVALIKE_M0$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BGB_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BGB_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BGB_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BGB_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Remove subset-sympatry
BGB_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BGB_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BGB_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BGB_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BGB_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BGB_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BGB_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BGB_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BGB_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BGB_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# Add jump dispersal/founder-event speciation
BGB_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BGB_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BGB_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
BGB_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BGB_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999

check_BioGeoBEARS_run(BGB_object)

resDIVALIKEj_M0<-BGB_run_default(runslow = TRUE,
                                 name_file="./Output/Solanum_DIVALIKE+J_M0_unconstrained_v1.Rdata", 
                                 BGB_object=BGB_object)



#######################################################
# BAYAREALIKE AND BAYAREALIKE+J ANALYSIS
#######################################################

# NOTE: As with DIVA, the BioGeoBEARS BayArea-like model is 
# not identical with the full Bayesian model implemented 
# in the "BayArea" program of Landis et al. (2013). 
#
# Instead, this is a simplified likelihood interpretation
# of the model.  Basically, in BayArea and BioGeoBEARS-BAYAREALIKE, 
# "d" and "e" work like they do in the DEC model of Lagrange 
# (and BioGeoBEARS), and then BayArea's cladogenesis assumption
# (which is that nothing in particular happens at cladogenesis) is 
# replicated by BioGeoBEARS.
#
# This leaves out 3 important things that are in BayArea:
# 1. Distance dependence (you can add this with a distances 
#    matrix + the "x" parameter in BioGeoBEARS, however)
# 2. A correction for disallowing "e" events that drive
#    a species extinct (a null geographic range)
# 3. The neat Bayesian sampling of histories, which allows
#    analyses on large numbers of areas.
#
# The main purpose of having a "BAYAREALIKE" model is 
# to test the importance of the cladogenesis model on 
# particular datasets. Does it help or hurt the data 
# likelihood if there is no special cladogenesis process?
# 
# I thus now call the model "BAYAREALIKE", and you should also. ;-)
# 

#######################################################
# Run BAYAREALIKE -----
#######################################################
BGB_object<-BGB_object_default(tree_path = "./Data/Solanum_biomes_tree.newick",
                               geo_data_path = "./Data/biomes_presence_matrix.data")

#BGB_object$areas_adjacency_fn = "./Data/Areas_adjacency.txt"
BGB_object<-readfiles_BioGeoBEARS_run(BGB_object)

# Set up BAYAREALIKE model
# No subset sympatry
BGB_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BGB_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BGB_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BGB_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BGB_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BGB_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# Adjust linkage between parameters
BGB_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BGB_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BGB_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BGB_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BGB_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BGB_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# Check the inputs
check_BioGeoBEARS_run(BGB_object)

resBAYAREALIKE_M0<-BGB_run_default(runslow = TRUE,
                                   name_file="./Output/Solanum_BAYEAREALIKE_M0_unconstrained_v1.Rdata", 
                                   BGB_object=BGB_object)



#######################################################
# Run BAYAREALIKE+J -----
#######################################################
BGB_object<-BGB_object_default(tree_path = "./Data/Solanum_biomes_tree.newick",
                               geo_data_path = "./Data/biomes_presence_matrix.data")

#BGB_object$areas_adjacency_fn = "./Data/Areas_adjacency.txt"
BGB_object<-readfiles_BioGeoBEARS_run(BGB_object)

# Set up BAYAREALIKE+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resBAYAREALIKE_M0$outputs@params_table["d","est"]
estart = resBAYAREALIKE_M0$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BGB_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BGB_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BGB_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BGB_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# No subset sympatry
BGB_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BGB_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BGB_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BGB_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BGB_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BGB_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
BGB_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BGB_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BGB_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
BGB_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

# Adjust linkage between parameters
BGB_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BGB_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BGB_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BGB_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BGB_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BGB_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
# machines. I can't replicate this on my Mac machines, but it is almost certainly
# just some precision under-run issue, when optim/optimx tries some parameter value 
# just below zero.  The "min" and "max" options on each parameter are supposed to
# prevent this, but apparently optim/optimx sometimes go slightly beyond 
# these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
# slightly for each parameter:
BGB_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
BGB_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999

BGB_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
BGB_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

BGB_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BGB_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

check_BioGeoBEARS_run(BGB_object)


resBAYAREALIKEj_M0<-BGB_run_default(runslow = TRUE,
                                    name_file="./Output/Solanum_BAYAREALIKE+J_M0_unconstrained_v1.Rdata", 
                                    BGB_object=BGB_object)

