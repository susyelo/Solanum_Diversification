################
## Libraries ---
################

library(OutbreakTools)
library(ape)
library(phytools)
library(magrittr)

################
## Functions ---
################
source("./Functions/annotated.tree.R")
source("./Functions/plotBranchbyTrait_mod.R")

################
## Data ---
################

taxonomy<-read.csv("./Data/Clades_solanum.csv")

## Sarkinen et al 2013 Solanum phylogeny
Sola_tree<-read.tree("Data/Solanum_Sarkinen_non_negative.tre")
tip_rates_all<-read.csv("./Output/Avg_rates_tips.csv")


################################################################################
## Calculate the mean of the div rate across all the runs for each species  ----
################################################################################
tip_rates_all$DivMean<-
  tip_rates_all%>%
  within(.,rm("X"))%>%
  rowMeans(.)


## Add clade information to the dataframe ---
tip_rates_all$Clade<-taxonomy$Clade[match(tip_rates_all$X,taxonomy$Species)]
tip_rates_all<-droplevels(tip_rates_all)

#############################################
## Plot the avg div rate in the phylogeny ---
#############################################

Sola_clades<-na.omit(unique(taxonomy$Clade))
Sola_clades<-Sola_clades[Sola_clades!="Jaltomata"]
Sola_clades<-droplevels(Sola_clades)

o<-match(Sola_tree$tip.label,tip_rates_all$X)
Trait_to_map<-tip_rates_all$DivMean[o]

Tips_nodes_rates<-c(Trait_to_map,rep(NA,nrow(Sola_tree$edge)/2))

for (i in 1:length(Sola_clades)){
  print(paste("Running", i, Sola_clades[i]))
  clade_tips<-tip_rates_all$X[which(tip_rates_all$Clade==Sola_clades[i])]
  tmp_node<-getMRCA(Sola_tree, as.character(clade_tips))
  if(!is.null(tmp_node)){
    clade_desc<-getDescendants(Sola_tree,node=tmp_node)
    Tips_nodes_rates[tmp_node]<-mean(tip_rates_all$DivMean[tip_rates_all$X%in%clade_tips])
    Tips_nodes_rates[clade_desc][is.na(Tips_nodes_rates[clade_desc])]<-mean(tip_rates_all$DivMean[tip_rates_all$X%in%clade_tips])
  }
  
}

Tips_nodes_rates[is.na(Tips_nodes_rates)]<-mean(tip_rates_all$DivMean[which(tip_rates_all$Clade!="old_world")])


## Creating table with the mean and sd div rates per each clade ---
extract_div_clade<-function(clade="old_world", runs=c(1:10)){
  tmp<-tip_rates_all[which(tip_rates_all$Clade==clade),match(paste("Run_",runs, sep=""),names(tip_rates_all))]
  mean_clade<-round(mean(colMeans(tmp)),4)
  sd_clade<-round(sd(colMeans(tmp)),4)
  clade_res<-data.frame(Clade=clade,mean_clade=mean_clade, sd_clade=sd_clade)
  return(clade_res)
}

DivRates_clades<-lapply(sort(unique(tip_rates_all$Clade)),extract_div_clade)
Clades_div_rates<-do.call(rbind,DivRates_clades)
write.csv(Clades_div_rates,"./Output/Clades_DivRates_RevBayes1_20.csv")
#write.csv(Clades_div_rates,"./Output/Clades_DivRates_RevBayes1_10.csv")

## Plotting rates ---
pdf("Figures/Avg_DivRates_RevBayes_Sarkinen_tree_20runs.pdf",  height = 9)
plotBranchbyTrait_mod(Sola_tree,Tips_nodes_rates,mode="predefined", legend=10,
                      show.tip.label = FALSE, palette="Zissou", title="Net Diversification rates",edge.width = 2)
dev.off()

