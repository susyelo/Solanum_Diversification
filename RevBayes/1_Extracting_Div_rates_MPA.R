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

################
## Data ---
################
taxonomy<-read.csv("./Data/Clades_solanum.csv")

## Sarkinen et al 2013 Solanum phylogeny
Sola_tree<-read.tree("Data/Solanum_Sarkinen_non_negative.tre")

runs<-list.files("./Data/Revbayes_runs/")


##########################################################################################################
## Creating a dataframe with the values of Div rates across all the runs per each Salanum species ---
##########################################################################################################
tip_rates_all<-NULL

for(i in 1:length(runs)){
tree_run<- read.annotated.nexus(file.path("./Data/RevBayes_runs/", runs[i], "/output/Solanum_BSBD_MAP.tree"))

### 1. Extracting avg rates values from all the branches ---
avg_index <- unlist(sapply(tree_run$annotations, function(e) e$index))
avg_mu_rates <- unlist(sapply(tree_run$annotations, function(e) e$avg_mu))
avg_lambda_rates <- unlist(sapply(tree_run$annotations, function(e) e$avg_lambda))
avg_div_rates<-avg_lambda_rates-avg_mu_rates

tip_index<-length(tree_run$tip.label)-1
## I don't understand why the 0 is including in the list and the order of the tips is length(tree1$tip.label)-1 to 0
names(avg_lambda_rates)<-avg_index[1:length(avg_lambda_rates)]

## 2. Extract the values only for the species ----
o<-match(as.character(tip_index:0),names(avg_lambda_rates))
tip_rates<-avg_lambda_rates[o]
names(tip_rates)<-tree_run$tip.label

## 3. Order by alphabetic order
tip_rates<-tip_rates[order(names(tip_rates))]

## 4. Add columns into the dataframe
tip_rates_all<-cbind(tip_rates_all,tip_rates)

}
colnames(tip_rates_all)<-runs

write.csv(tip_rates_all, "./Output/Avg_rates_tips.csv")
