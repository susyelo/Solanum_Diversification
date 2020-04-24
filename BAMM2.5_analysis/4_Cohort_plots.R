####################
#### LIBRARIES ----
####################
library(BAMMtools)

####################
#### FUNCTIONS ----
####################

source("./Functions/barLegend.R")
source("./Functions/plotBranchbyTrait_mod.R")

load_obj <- function(f)
{
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

## Get cohort matrix per run ---
get_cohort<-function(ed_runs){
  print(paste("Processing runs"))
  cmat<-getCohortMatrix(ed_runs)
  return(cmat)
}

## sort each cohort matrix by species names
sort_tips<-function(cmat_obj){
  print(paste("Processing runs"))
  cmat_obj<-cmat_obj[order(rownames(cmat_obj)),]
  return(cmat_obj)
}


####################
#### DATA ----
####################
## Reading trees for each run ---
runs<-as.numeric(list.files("./Data/Input_BAMM2.5/"))

#Reading tiprates files estimating using BAMMtools
Avg_tip_rates<-read.csv("./Output/Avg_TipsRates.csv")

#Reading the information about the nodes associated with rates of diversification in BAMMtools
nodes_info<-read.csv("./Output/Nodes_with_shifts.csv")
Sola_tree<-read.tree("./Data/Solanum_Sarkinen_non_negative.tre")


trees_run <- NULL

for(i in 1:length(runs)){
  tree<-read.tree(file.path("./Data/Input_BAMM2.5/",runs[i], "tree.tre"))
  run_tmp<-paste("Run_", runs[i], sep="")
  trees_run[[run_tmp]]<-tree
}
class(trees_run) <- "multiPhylo"

## Reading event data from each of the pastis runs
### This will take a while ... 
file_names=list.files("./Output/", pattern="*.RData", full.names = TRUE)
ed_runs<-lapply(file_names, load_obj)
names(ed_runs)<- names(trees_run)




####################
#### MAIN ----
####################
## Getting macroevolutionary cohort information per PASTIS run
cmat<-lapply(ed_runs, get_cohort)

## Sorting the data for each run to make sure there is a standard comparison across runs
cmat<-lapply(cmat, sort_tips)

### Calculating the mean of the values across the matrices 
avg_cmat<-Reduce("+", cmat) / length(cmat)

### Subsetting to the Sarkinen phylogeny 
index_1 <- match(Sola_tree$tip.label,rownames(avg_cmat))
index_2 <- match(Sola_tree$tip.label,colnames(avg_cmat))
avg_cmat_sarkinen<- avg_cmat[index_1, index_2]



## Calculating the value of the rates at the nodes. This value will be the average rates for each of the clades. 
#The rest of the edges will have the average value of diverisification rates of Solanum without taking into Old world clade

## Subset dataframe to the tips in Sarkinen et al phylogeny ---
tiprates_sarkinen<-Avg_tip_rates[which(Avg_tip_rates$X%in%Sola_tree$tip.label),]
variables<-grep("Run",colnames(tiprates_sarkinen))

##Calculating the median values of rates across all the 100 runs of BAMM
tiprates_sarkinen$median<-apply(tiprates_sarkinen[,variables],1,median)

## Ordering the rates values using the Sola_tree as reference
o<-match(Sola_tree$tip.label,tiprates_sarkinen$X)
Trait_to_map<-as.numeric(tiprates_sarkinen$mean[o])
names(Trait_to_map)<-tiprates_sarkinen$X[o]


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



breaks <- seq(min(avg_cmat_sarkinen), max(avg_cmat_sarkinen),length.out=100)
#cols<-wes_palette("Zissou", 100, type = "continuous")
cols<-colorRampPalette(brewer.pal(5,"PuBuGn"))(100)

### Plotting the macroevolutionary cohort using 
pdf("Figures/Avg_cohort.pdf", width = 11)
ofs = 0
figs <- matrix(c(0, 0.2, 0.8, 1, 0.2, 0.95, 0.8 + ofs, 1, 
                 0, 0.2 - ofs, 0, 0.8, 0.2, 0.95, 0, 0.8, 0.98, 1, 0.25, 
                 0.75), byrow = TRUE, nrow = 5, ncol = 4)

par(fig = figs[2, ], new = FALSE, mar = c(0, 0, 1, 4))
plotBranchbyTrait_mod(Sola_tree,Tips_nodes_rates,mode="predefined", legend=FALSE,
                      show.tip.label = FALSE, palette="Zissou", title="Net Diversification rates",edge.width = 3,direction = "downwards")

par(fig = figs[3, ], new = TRUE, mar = c(5, 1, 0, 0))
plotBranchbyTrait_mod(Sola_tree,Tips_nodes_rates,mode="predefined", legend=FALSE,
                      show.tip.label = FALSE, palette="Zissou", title="Net Diversification rates",edge.width = 2,direction = "rightwards")


par(fig = figs[4, ], new = TRUE, mar = c(0, 0, 0, 0))
plot(0, 0, type = "n", axes = FALSE, ann = FALSE, xlim = c(0,1), ylim = c(0, 1))

image(avg_cmat_sarkinen, axes = FALSE, xlab = "", ylab = "", col = cols, xlim = c(0, 1), ylim = c(0, 1), add = TRUE, useRaster = TRUE)

barLegend(cols, breaks, fig = figs[5, ], side = 2)

dev.off()
