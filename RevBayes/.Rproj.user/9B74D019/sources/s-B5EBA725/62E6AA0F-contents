## LIB ----
library(coda)
library(ape)
library(BAMMtools)

## Functions ---
nodes_in_Sarkinen<-function(tree_runs, tree_sarkinen, t_nodes)
{
  node_shifts<-NULL
  for (j in 1:length(t_nodes))
  { 
    tips_shift<-getDescendants(as.phylo(tree_runs),t_nodes[j])
    tips_names<-na.omit(ed_runs[[i]]$tip.label[tips_shift])
    node_tmp<-getMRCA(tree_sarkinen,Sarkinen_tree$tip.label[which(Sarkinen_tree$tip.label%in%tips_names)])
    
    if (!is.null(node_tmp)){
      if(node_tmp!=getMRCA(tree_sarkinen,tree_sarkinen$tip.label)){
        node_shifts<-c(node_shifts,node_tmp)
      }
      
    }
  }
  return(node_shifts)
}


## FILES --- 

## Reading mcmc files and calculating the effective size of each run -- 
ESS_shifts<-NULL
ESS_loglik<-NULL
runs<-as.numeric(list.files("./Data/Input_BAMM2.5/"))

for(i in 1:length(runs)){
  mcmcout<-read.csv(file.path("./Data/Input_BAMM2.5/",runs[i], "mcmc_out.txt"))
  #print(i)
  #discard some burning (here 10%)
  burnstart<-floor(.25*nrow(mcmcout))
  postburn<-mcmcout[burnstart:nrow(mcmcout),]
  ESS_shifts[i]<-effectiveSize(postburn$N_shifts)
  ESS_loglik[i]<-effectiveSize(postburn$logLik)
}


Eff_size<-data.frame(Run=runs, ESS_loglik, ESS_shifts)

## Reading phylogenetic trees ---
trees_run <- NULL

for(i in 1:length(Eff_size$Run)){
  tree<-read.tree(file.path("./Data/Input_BAMM2.5/",Eff_size$Run[i], "tree.tre"))
  run_tmp<-paste("Run_", Eff_size$Run[i], sep="")
  trees_run[[run_tmp]]<-tree
}
class(trees_run) <- "multiPhylo"

## extract mcmc values post burn --
mcmc_out<-lapply(runs,
                 function(x)read.csv(file.path("./Data/Input_BAMM2.5/",x, "mcmc_out.txt")))

mcmc_post<-lapply(mcmc_out,
                  function(x)x[floor(.25*nrow(x)):nrow(x),])

## Creating pseudoposterior samples
library(plyr)
mcmc_combined<-ldply(mcmc_post, data.frame)

## Plotting prior vs posterior probabilities of the number of shifts along the tree
source("./Functions/plotPrior_mod.R")

pdf("./Figures/Priors_posteriors_BAMM.pdf")
plotPrior_mod(mcmc_combined, expectedNumberOfShifts=1)
legend("topright",legend=c("Prior","Posterior"),fill=c("red","lightblue"))
title(xlab="Number of shifts", ylab="Probabilities")
dev.off()


### Plot the prior vs posterior distribution for all the BAMM runs

pdf("./Figures/Priors_posteriors_BAMM_1_25.pdf", width = 10)
m<- matrix(1:25,ncol=5,nrow=5,byrow=TRUE)
layout(m)

for(i in 1:25) {
  par(mar = c(3, 3, 2, 0))
  plotPrior_mod(mcmc_post[[i]], expectedNumberOfShifts=1)
  title(main=paste("Tree", i),cex=0.5)
}
dev.off()

pdf("./Figures/Priors_posteriors_BAMM_26_50.pdf", width = 10)
m<- matrix(1:25,ncol=5,nrow=5,byrow=TRUE)
layout(m)

for(i in 26:50) {
  par(mar = c(3, 3, 2, 0))
  plotPrior_mod(mcmc_post[[i]], expectedNumberOfShifts=1)
  title(main=paste("Tree", i),cex=0.5)
}
dev.off()

pdf("./Figures/Priors_posteriors_BAMM_51_75.pdf", width = 10)
m<- matrix(1:25,ncol=5,nrow=5,byrow=TRUE)
layout(m)

for(i in 51:75) {
  par(mar = c(3, 3, 2, 0))
  plotPrior_mod(mcmc_post[[i]], expectedNumberOfShifts=1)
  title(main=paste("Tree", i),cex=0.5)
}
dev.off()

pdf("./Figures/Priors_posteriors_BAMM_76_100.pdf", width = 10)
m<- matrix(1:25,ncol=5,nrow=5,byrow=TRUE)
layout(m)

for(i in 76:100) {
  par(mar = c(3, 3, 2, 0))
  plotPrior_mod(mcmc_post[[i]], expectedNumberOfShifts=1)
  title(main=paste("Tree", i),cex=0.5)
}
dev.off()



#source("./Functions//computeBayesfactors_mod.R")
## Bayes factors cannot be calculated since the null model (no shifts) was never sampled in the posterior samples
#BF<-computeBayesFactors_mod(mcmc_combined, expectedNumberOfShifts=1, burnin=0.1)

## Reading event data 
events.sol<-lapply(runs,
                   function(x)read.csv(file.path("./Data/Input_BAMM2.5/",x, "event_data.txt")))


## calculate diversification variables  ----
ed_runs<-NULL

for(i in 1:length(runs)){
  print(paste("Processing run",runs[i]))
  ed <- getEventData(trees_run[[names(trees_run)[i]]], events.sol[[i]], burnin=0.2, nsample=1000)
  ed_runs[[names(trees_run)[i]]]<-ed
  save(ed, file=paste0(getwd(),'/Output/Tree', New_run[i], '.RData'))
}

## Read RData for all the event data ---
load_obj <- function(f)
{
env <- new.env()
nm <- load(f, env)[1]
env[[nm]]
}

##file_names=list.files("./Output/", pattern="*.RData", full.names = TRUE)
##ed_runs<-lapply(file_names, load_obj)
##names(ed_runs)<- names(trees_run)

## Get tip rates ---
get_tip_rates<-function(ed_runs){
  print(paste("Processing runs"))
  tips_rates<-getTipRates(ed_runs, returnNetDiv = TRUE, statistic = "mean")$netdiv.avg
  return(tips_rates)
}

tip_rates<-lapply(ed_runs, get_tip_rates)

## Sorting tip rates
re_order<-function(x=1){
  new_order<-tip_rates[[x]][names(tip_rates[[1]])]
  return(new_order)
}

sort_tip_rates<-lapply(1:100,re_order)

tip_rates_df<-as.data.frame(do.call("cbind", sort_tip_rates))
colnames(tip_rates_df)<-names(tip_rates)

## ploting distribution of rates between old world and new world ---
taxonomy<-read.csv("./Data/Clades_solanum.csv")
tip_rates_df$Clade<-taxonomy$Clade[match(rownames(tip_rates_df), taxonomy$Species)]
tip_rates_df$Cont<-ifelse(tip_rates_df$Clade=="old_world", "old_world", "Other")

tip_rates_df$mean<-apply(tip_rates_df[,1:100],1,mean)

## Saving Average tip rates
write.csv(tip_rates_df, file="./Output/Avg_TipsRates.csv")

###############################################
##### Plotting avg rates on the phylogeny 
#############################################

### Get branchshiftprior
branchShift_prior<-NULL

for(i in 1:length(Eff_size$Run)){
  branch_shift_prior <- getBranchShiftPriors(trees_run[[names(trees_run)[i]]], expectedNumberOfShifts=1)
  branchShift_prior[[names(trees_run)[i]]]<-branch_shift_prior
}

## Locating the node with the most-probable shift configuration
library(phytools)
Sarkinen_tree<-read.tree("./Data/Solanum_Sarkinen_non_negative.tre")


## Identifying the nodes supporting the credible shift set from each run
tips_name_shift<-NULL
for(i in 1:length(ed_runs)){
  
  print(paste("Processing run",i))
  css <- credibleShiftSet(ed_runs[[i]], expectedNumberOfShifts=1, threshold=50, set.limit = 0.95)
  node_tmp<-css$shiftnodes[[1]]
  tips_name_shift[[names(ed_runs[i])]]<-lapply(node_tmp, 
                                                 function(x)
                                                   na.omit(css$tip.label[getDescendants(trees_run[[i]], x)]))
}

## Finding nodes in Sarkinen with the most-probable shift configuration
node_bf<-NULL
for(i in 1:length(tips_name_shift)){
  
  node_bf[[names(tips_name_shift[i])]]<-lapply(tips_name_shift[[i]], 
                                               function(x)
                                                 getMRCA(Sarkinen_tree,as.character(x))
                                               )
}

nodes_info<-as.data.frame(table(unlist(node_bf)))
nodes_info$bf_freq<-round(nodes_info$Freq/100,2)
nodes_info$bf_freq_char<-paste(nodes_info$Freq, "/", length(ed_runs), sep="")
write.csv(nodes_info, "./Output/Nodes_with_shifts.csv")

## Data for the LTT ---
Rtt_matrix<-NULL
Rtt_ow<-NULL
Rtt_petota<-NULL
Rtt_Torva<-NULL
Rtt_Geminata<-NULL
Rtt_Tomato<-NULL


for(i in 1:length(ed_runs)){
  tips_clade<-trees_run[[i]]$tip.label[which(trees_run[[i]]$tip.label%in%taxonomy$Species[which(taxonomy$Clade=="old_world")])]
  node_clade<-getMRCA(trees_run[[i]], tips_clade)
  
  petota_tips<-trees_run[[i]]$tip.label[which(trees_run[[i]]$tip.label%in%taxonomy$Species[which(taxonomy$Clade=="Petota")])]
  petota_node<-getMRCA(trees_run[[i]], petota_tips)
  
  Torva_tips<-trees_run[[i]]$tip.label[which(trees_run[[i]]$tip.label%in%taxonomy$Species[which(taxonomy$Clade=="Torva")])]
  Torva_node<-getMRCA(trees_run[[i]], Torva_tips)
  
  Geminata_tips<-trees_run[[i]]$tip.label[which(trees_run[[i]]$tip.label%in%taxonomy$Species[which(taxonomy$Clade=="Geminata")])]
  Geminata_node<-getMRCA(trees_run[[i]], Geminata_tips)
  
  Tomato_tips<-trees_run[[i]]$tip.label[which(trees_run[[i]]$tip.label%in%taxonomy$Species[which(taxonomy$Clade=="Tomato")])]
  Tomato_node<-getMRCA(trees_run[[i]], Tomato_tips)
  
  
  print(names(ed_runs)[i])
  Rtt_ow[[names(ed_runs)[i]]] <- getRateThroughTimeMatrix(ed_runs[[i]], node = node_clade,nslices = 50)
  Rtt_petota[[names(ed_runs)[i]]] <- getRateThroughTimeMatrix(ed_runs[[i]], node = petota_node,nslices = 50)
  Rtt_Torva[[names(ed_runs)[i]]] <- getRateThroughTimeMatrix(ed_runs[[i]], node = Torva_node,nslices = 50)
  Rtt_Geminata[[names(ed_runs)[i]]] <- getRateThroughTimeMatrix(ed_runs[[i]], node = Geminata_node,nslices = 50)
  Rtt_Tomato[[names(ed_runs)[i]]] <- getRateThroughTimeMatrix(ed_runs[[i]], node = Tomato_node,nslices = 50)
  
  Rtt_matrix[[names(ed_runs)[i]]]<-getRateThroughTimeMatrix(ed_runs[[i]],nslices = 50)
}

save(Rtt_matrix,file="./Output/RatesThroughTime.Rdata")
save(Rtt_ow,file="./Output/RatesThroughTime_ow.Rdata")
save(Rtt_petota,file="./Output/RatesThroughTime_Petota.Rdata")
save(Rtt_Torva,file="./Output/RatesThroughTime_Torva.Rdata")
save(Rtt_Geminata,file="./Output/RatesThroughTime_Geminata.Rdata")
save(Rtt_Tomato,file="./Output/RatesThroughTime_Tomato.Rdata")



