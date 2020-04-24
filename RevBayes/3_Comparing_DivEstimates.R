################
## Libraries ---
################


################
## Functions ---
################
CI_DivRates<-function(data, data_mean, data_sd)
{
  data$upper<-data_mean+data_sd*1.96
  data$lower<-data_mean-data_sd*1.96
  return(data)
}


################
## Data ---
################
DivRates_BAMM1_10<-read.csv("./Data/Clades_DivRates_BAMM1_10.csv")
DivRates_BAMM1_100<-read.csv("./Data/Clades_DivRates_BAMM1_100.csv")
DivRates_RevBayes<-read.csv("Output/Clades_DivRates_RevBayes1_10.csv")

DivRates_BAMM1_10<-CI_DivRates(DivRates_BAMM1_10, DivRates_BAMM1_10$mean_clade, DivRates_BAMM1_10$sd_clade)
DivRates_BAMM1_100<-CI_DivRates(DivRates_BAMM1_100, DivRates_BAMM1_100$mean_clade, DivRates_BAMM1_100$sd_clade)
DivRates_RevBayes<-CI_DivRates(DivRates_RevBayes, DivRates_RevBayes$mean_clade, DivRates_RevBayes$sd_clade)


### Plotting estimates for specific clades: Old world, Petota, Torva, Thelopodium 
clades_to_plot<-c("old_world","Petota","Torva","Thelopodium","Geminata")
o<-match(clades_to_plot,as.character(DivRates_RevBayes$Clade))

plot(DivRates_RevBayes$mean_clade[o],1:5, xlim=c(0,1.4),ylim=c(1,6),xlab='', ylab='',yaxt='n')
axis(2,at=1:5,labels=clades_to_plot)
segments(DivRates_RevBayes$lower[o],1:5,DivRates_RevBayes$upper[o])

o<-match(clades_to_plot,as.character(DivRates_BAMM1_10$Clade))

points(DivRates_BAMM1_10$mean_clade[o],c(1.2,2.2,3.2,4.2,5.2))
segments(DivRates_BAMM1_10$lower[o],c(1.2,2.2,3.2,4.2,5.2),DivRates_BAMM1_10$upper[o],col="red")
