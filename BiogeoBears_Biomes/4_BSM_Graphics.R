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


### Functions ----
## Loading RData function
load_obj <- function(f)
{
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
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


#########################################################
# General data ----
########################################################
tr<-read.tree("./Data/Solanum_corrected.newick")
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn="./Data/Solanum_geo_2.data")


## Open the best model
resDECj_M1<-load_obj("./Output/Solanum_DEC+J_M1_2_v1.Rdata")

#######################################################
# PDF plots------
#######################################################
pdffn = "./Figures/Solanum_DECJ_M1.pdf"
pdf(pdffn, width=20, height=100)
plot_DEC_M1<-plot_BGB_two_models(resDECj_M1,Title = "BioGeoBEARS DEC Solanum constrain dispersal matrix")
dev.off()  # Turn off PDF



# ---------- Biogeographical Stochastic Mapping -------------

#######################################################
# Stochastic mapping on DEC
#######################################################
clado_events_tables = NULL
ana_events_tables = NULL
lnum = 0

#######################################################
# Get the inputs for Biogeographical Stochastic Mapping
# Note: this can be slow for large state spaces and trees, since 
# the independent likelihoods for each branch are being pre-calculated
# E.g., for 10 areas, this requires calculation of a 1024x1024 matrix
# for each branch.  On a tree with ~800 tips and thus ~1600 branches, this was about 1.6 gigs
# for storage of "BSM_inputs_file.Rdata".
# Update: 2015-09-23 -- now, if you used multicore functionality for the ML analysis,
# the same settings will be used for get_inputs_for_stochastic_mapping().
#######################################################
BSM_inputs_fn = "./Output/BSM_inputs_file.Rdata"
runInputsSlow = TRUE

##Changing the areas to one letter names

resDEC_M1$inputs$geogfn<-"./Data/Solanum_geo_2.data"


if (runInputsSlow)
{
  stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=resDECj_M1)
  save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
} else {
  # Loads to "stochastic_mapping_inputs_list"
  load(BSM_inputs_fn)
} # END if (runInputsSlow)


# Check inputs (doesn't work the same on unconstr)
names(stochastic_mapping_inputs_list)
stochastic_mapping_inputs_list$phy2
stochastic_mapping_inputs_list$COO_weights_columnar
stochastic_mapping_inputs_list$unconstr
set.seed(seed=as.numeric(Sys.time()))




# Extract BSM output
clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables
head(clado_events_tables[[2]])
head(ana_events_tables[[1]])
length(clado_events_tables)
length(ana_events_tables)


#######################################################
# Plot one stochastic map, manual method
#######################################################
# (we have to convert the stochastic maps into event
#  maps for plotting)

######################
# Get the color scheme
######################
include_null_range = TRUE
areanames = names(tipranges@df)
areas = areanames
max_range_size = 3

# Note: If you did something to change the states_list from the default given the number of areas, you would
# have to manually make that change here as well! (e.g., areas_allowed matrix, or manual reduction of the states_list)
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)

colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, 
                                                           states_list_0based=states_list_0based, 
                                                           max_range_size=max_range_size, 
                                                           plot_null_range=TRUE)

############################################
# Setup for painting a single stochastic map
############################################
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
stratified=FALSE
clado_events_table = clado_events_tables[[1]]
ana_events_table = ana_events_tables[[1]]

# cols_to_get = names(clado_events_table[,-ncol(clado_events_table)])
# colnums = match(cols_to_get, names(ana_events_table))
# ana_events_table_cols_to_add = ana_events_table[,colnums]
# anagenetic_events_txt_below_node = rep("none", nrow(ana_events_table_cols_to_add))
# ana_events_table_cols_to_add = cbind(ana_events_table_cols_to_add, anagenetic_events_txt_below_node)
# rows_to_get_TF = ana_events_table_cols_to_add$node <= length(tr$tip.label)
# master_table_cladogenetic_events = rbind(ana_events_table_cols_to_add[rows_to_get_TF,], clado_events_table)

############################################
# Open a PDF
############################################
model_name = "DECjM1"
pdffn = paste0("./Figures/",model_name, "_single_stochastic_map_n1.pdf")
pdf(file=pdffn,width=20, height=100)

# Convert the BSM into a modified res object
master_table_cladogenetic_events = clado_events_tables[[1]]
resmod = stochastic_map_states_into_res(res=resDECj_M1, master_table_cladogenetic_events=master_table_cladogenetic_events)

plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", 
                         addl_params=list("j"), label.offset=0.5, 
                         plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, 
                         colors_list_for_states=colors_list_for_states, skiptree=FALSE, show.tip.label=TRUE)

# Paint on the branch states
paint_stochastic_map_branches(res=resmod, master_table_cladogenetic_events=master_table_cladogenetic_events, 
                              colors_list_for_states=colors_list_for_states, lwd=5, lty=par("lty"), root.edge=TRUE, stratified=FALSE)

plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), 
                         plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, 
                         colors_list_for_states=colors_list_for_states, skiptree=TRUE, show.tip.label=TRUE)

############################################
# Close PDF
############################################
dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)



#######################################################
# Plot all 50 stochastic maps to PDF
#######################################################
# Setup
include_null_range = include_null_range
areanames = areanames
areas = areanames
max_range_size = max_range_size
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)
colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=TRUE)
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
stratified = FALSE

# Loop through the maps and plot to PDF
pdffn = paste0("./Figures/",model_name, "_", length(clado_events_tables), "BSMs_v1.pdf")
pdf(file=pdffn, width=20, height=100)

nummaps_goal = 50
for (i in 1:nummaps_goal)
{
  clado_events_table = clado_events_tables[[i]]
  analysis_titletxt = paste0(model_name, " - Stochastic Map #", i, "/", nummaps_goal)
  plot_BSM(results_object=res, clado_events_table=clado_events_table, stratified=FALSE, analysis_titletxt=analysis_titletxt, addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, show.tip.label=TRUE, include_null_range=include_null_range)
} # END for (i in 1:nummaps_goal)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)


#######################################################
# Summarize stochastic map tables
#######################################################
length(clado_events_tables)
length(ana_events_tables)

head(clado_events_tables[[1]][,-20])
tail(clado_events_tables[[1]][,-20])

head(ana_events_tables[[1]])
tail(ana_events_tables[[1]])

areanames = names(tipranges@df)
actual_names = areanames
actual_names

# Get the dmat and times (if any)
dmat_times = get_dmat_times_from_res(res=resDECj_M1, numstates=NULL)
dmat_times

# Simulate the source areas
BSMs_w_sourceAreas = simulate_source_areas_ana_clado(resDECj_M1, clado_events_tables, ana_events_tables, areanames)
clado_events_tables = BSMs_w_sourceAreas$clado_events_tables
ana_events_tables = BSMs_w_sourceAreas$ana_events_tables

# Count all anagenetic and cladogenetic events
counts_list = count_ana_clado_events(clado_events_tables, ana_events_tables, areanames, actual_names)

save(counts_list,file="./Output/Event_count_list_DECjM1.RData")

##To find the cladogenetic events in an specifiic node (in this case 20) in an specific table (in this case 23)
### Checking the cladogenetic events in the old world clade ----
Clades_in_phylo<-read.csv("./Data/Clades_in_phylo.csv")

sp_old_world<-Clades_in_phylo$species[which(Clades_in_phylo$clades=="Old_world")]

ow_tree_tips<-tr$tip.label[which(tr$tip.label%in%sp_old_world)]
ow_node<-getMRCA(tr,ow_tree_tips)

rootnode_row_TF = clado_events_tables[[1]]$node == ow_node
clado_events_tables[[45]][rootnode_row_TF,]
# print the actual event (several equivalent methods):
clado_events_tables[[70]]$clado_event_txt[rootnode_row_TF]


summary_counts_BSMs = counts_list$summary_counts_BSMs
print(conditional_format_table(summary_counts_BSMs))

# Histogram of event counts
hist_event_counts(counts_list, pdffn=paste0("./Figures/",model_name, "_histograms_of_event_counts.pdf"))
#######################################################


