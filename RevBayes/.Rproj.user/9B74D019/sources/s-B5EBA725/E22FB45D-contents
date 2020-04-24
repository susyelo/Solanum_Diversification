## Libraries ---
library(BAMMtools)
library(ggtree)
library(magrittr)


##sources ---
source("./Functions/plotBranchbyTrait_mod.R")

## Data ---
#Reading tiprates files estimating using BAMMtools
Avg_tip_rates<-read.csv("./Output/Avg_TipsRates.csv")

#Reading the information about the nodes associated with rates of diversification in BAMMtools
nodes_info<-read.csv("./Output/Nodes_with_shifts.csv")

## Sarkinen et al 2013 Solanum phylogeny
Sola_tree<-read.tree("Data/Solanum_Sarkinen_non_negative.tre")

## Subset dataframe to the tips in Sarkinen et al phylogeny ---
tiprates_sarkinen<-Avg_tip_rates[which(Avg_tip_rates$X%in%Sola_tree$tip.label),]
variables<-grep("Run",colnames(tiprates_sarkinen))

##Calculating the median values of rates across all the 100 runs of BAMM
tiprates_sarkinen$median<-apply(tiprates_sarkinen[,variables],1,median)


## Reading trees for each run ---
runs<-as.numeric(list.files("./Data/Input_BAMM2.5/"))

trees_run <- NULL

for(i in 1:length(runs)){
  tree<-read.tree(file.path("./Data/Input_BAMM2.5/",runs[i], "tree.tre"))
  run_tmp<-paste("Run_", runs[i], sep="")
  trees_run[[run_tmp]]<-tree
}
class(trees_run) <- "multiPhylo"



## Ordering the rates values using the Sola_tree as reference
o<-match(Sola_tree$tip.label,tiprates_sarkinen$X)
Trait_to_map<-as.numeric(tiprates_sarkinen$mean[o])
names(Trait_to_map)<-tiprates_sarkinen$X[o]

## Nodes that appear more than 10% of the times
sig_nodes<-nodes_info[which(nodes_info$bf_freq>0.10),]

## Calculating the value of the rates at the nodes. This value will be the average rates for each of the clades.
#The rest of the edges will have the average value of diverisification rates of Solanum without taking into Old world clade
Tips_nodes_rates<-c(as.numeric(Trait_to_map),rep(NA,nrow(Sola_tree$edge)/2))

clades_sola<-unique(tiprates_sarkinen$Clade)

## For each clade calculates the average rates values and assings to the edges associated with the clades.
for (i in 1:length(clades_sola)){
  print(paste("Running", i, clades_sola[i]))
  clade_tips<-tiprates_sarkinen$X[which(tiprates_sarkinen$Clade==clades_sola[i])]
  tmp_node<-getMRCA(Sola_tree, as.character(clade_tips))
  if(!is.null(tmp_node)){
    clade_desc<-getDescendants(Sola_tree,node=tmp_node)
    Tips_nodes_rates[tmp_node]<-mean(Tips_nodes_rates[clade_desc],na.rm=TRUE)
    Tips_nodes_rates[clade_desc][is.na(Tips_nodes_rates[clade_desc])]<-mean(Tips_nodes_rates[clade_desc],na.rm=TRUE)
  }

}

Tips_nodes_rates[is.na(Tips_nodes_rates)]<-mean(tiprates_sarkinen$mean[which(tiprates_sarkinen$Clade!="old_world")])


pdf("./Figures/Avg_rates_Sarkinen_phylo_ladde.pdf", height = 9)
plotBranchbyTrait_mod(ladderize(Sola_tree),Tips_nodes_rates,mode="predefined", legend=7,
                      show.tip.label = FALSE, palette="Zissou", title="Net Diversification rates",edge.width = 2)

nodelabels(thermo=sig_nodes$Freq/100,node=as.numeric(as.character(sig_nodes$Var1)), piecol=c("grey","white"),width = 0.2)
dev.off()

## Reading pc2
pc2<-read.csv("~/Box Sync/Susy Echeverria-Londono/Solanaceae/Analysis/Environment_analysis/Output/pc2.csv")

Sola_tree$pc2<-pc2$pc2

layout(matrix(c(1,2), 1,2,byrow=T))
par(mar=c(0,0,0,0))
plotBranchbyTrait_mod(Sola_tree,Tips_nodes_rates,mode="predefined", legend=7,
                      show.tip.label = FALSE, palette="Zissou", title="Net Diversification rates",edge.width = 2)

nodelabels(thermo=sig_nodes$Freq/100,node=as.numeric(as.character(sig_nodes$Var1)), piecol=c("grey","white"),width = 0.2)
barplot(Sola_tree$pc2, horiz = TRUE, axisname=FALSE)



extract_div_clade<-function(clade="old_world", runs=c(1:100)){
  tmp<-Avg_tip_rates[which(Avg_tip_rates$Clade==clade),match(paste("Run_",runs, sep=""),names(Avg_tip_rates))]
  mean_clade<-round(mean(colMeans(tmp)),4)
  sd_clade<-round(sd(colMeans(tmp)),4)
  clade_res<-data.frame(Clade=clade,mean_clade=mean_clade, sd_clade=sd_clade)
  return(clade_res)
}

DivRates_clades<-lapply(sort(unique(Avg_tip_rates$Clade)),extract_div_clade)
Clades_div_rates<-do.call(rbind,DivRates_clades)
write.csv(Clades_div_rates,"./Output/Clades_DivRates_BAMM1_100.csv")


DivRates_clades1_10<-lapply(sort(unique(Avg_tip_rates$Clade)),function(x)extract_div_clade(x,runs=1:10))
Clades_div_rates1_10<-do.call(rbind,DivRates_clades1_10)
write.csv(Clades_div_rates1_10,"./Output/Clades_DivRates_BAMM1_10.csv")

DivRates_clades1_20<-lapply(sort(unique(Avg_tip_rates$Clade)),function(x)extract_div_clade(x,runs=1:20))
Clades_div_rates1_20<-do.call(rbind,DivRates_clades1_20)
write.csv(Clades_div_rates1_20,"./Output/Clades_DivRates_BAMM1_20.csv")



##Creating density plots
library(HistogramTools)
PlotRelativeFrequency(hist(Avg_tip_rates$mean, nclass = 100, xlim=c(0,0.8), ylim = c(0,1)))
