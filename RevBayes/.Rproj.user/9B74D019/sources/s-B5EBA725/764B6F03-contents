#### Plotting average rates of diversification ---

##############
### LIBRARIES ----
##############

library(ape)
library(phytools)
library(paleotree)
library(wesanderson)


##############
### FUNCTIONS ----
##############

### Diversification rates through time for all the runs ---
plot_DivRates_uncertainty<-function(Rates_matrix, ylim=c(0,1))
{
  Rtt_DIV<-lapply(Rates_matrix, function(x)(x$lambda-x$mu))
  meanDivRate <- lapply(Rtt_DIV, function(x)colMeans(x))
  
  MEAN_DIV<-apply(simplify2array(meanDivRate), 1, mean)
  SD_DIV<-apply(simplify2array(meanDivRate), 1, sd)
  
  Div_plus<-MEAN_DIV+SD_DIV*1.96
  Div_minus<-MEAN_DIV-SD_DIV*1.96
  
  time_rates<-Rates_matrix[[1]]$times-max(Rates_matrix[[1]]$times)
  
  plot(1, type="n", xlab="Time", ylab="Diversification rates", 
       xlim=c(-15.6,0),ylim=ylim)
  
  X.Vec<-c(time_rates, max(time_rates), 
           rev(time_rates), min(time_rates))
  
  Y.Vec <- c(((Div_minus)), 
             tail((Div_plus), 1), 
             rev((Div_plus)), 
             ((Div_minus))[1])
  
  polygon(X.Vec, Y.Vec, col = "grey", border = NA)
  points(time_rates,MEAN_DIV,type="l",lwd=1)
  
}

## Plot LTT 
plot_ltt<-function(trees_ltt, col_curve="grey"){
  X.Vec<-c(trees_ltt$int.times[,1], min(trees_ltt$int.times[,1]), 
           rev(trees_ltt$int.times[,1]), max(trees_ltt$int.times[,1]))
  
  Y.Vec <-log(c(((trees_ltt$median.curve[,2])), 
                tail((trees_ltt$median.curve[,3]), 1), 
                rev((trees_ltt$median.curve[,3])), 
                ((trees_ltt$median.curve[,2]))[1])+1)
  
  polygon(X.Vec, Y.Vec, col = col_curve, border = NA)
  points(trees_ltt$int.times[,1],log(trees_ltt$median.curve[,1]+1),type="l",lwd=1)
  
}


## Transparent colors
## Mark Gardener 2015
## www.dataanalytics.org.uk

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}
## END

##############
### DATA ----
##############

### Reading trees for each of the BAMM runs ----
runs<-as.numeric(list.files("./Data/Input_BAMM2.5/"))
trees_run <- NULL

for(i in 1:length(runs)){
  tree<-read.tree(file.path("./Data/Input_BAMM2.5/",runs[i], "tree.tre"))
  run_tmp<-paste("Run_", runs[i], sep="")
  trees_run[[run_tmp]]<-tree
}
class(trees_run) <- "multiPhylo"


##############
### MAIN ----
##############

### Lineages through time comparision among continents ----
Neo_trees<-NULL
AU_trees<-NULL
AF_trees<-NULL


Cont_info<-read.csv("./Data/Solanum_geo_realm.csv")

## Subsetting by continent
Neo<-Cont_info$species[which(Cont_info$Cox_area=="Neotropic")]
AU<-Cont_info$species[which(Cont_info$Cox_area=="Australian")]
AF<-Cont_info$species[which(Cont_info$Cox_area=="African")]

## subset each of the trees with the lineages of each of the continents
for (i in 1:length(trees_run)){
  Neo_trees[[i]]<-drop.tip(trees_run[[i]],trees_run[[i]]$tip[which(trees_run[[i]]$tip.label%in%Neo==FALSE)])
  AU_trees[[i]]<-drop.tip(trees_run[[i]],trees_run[[i]]$tip[which(trees_run[[i]]$tip.label%in%AU==FALSE)])
  AF_trees[[i]]<-drop.tip(trees_run[[i]],trees_run[[i]]$tip[which(trees_run[[i]]$tip.label%in%AF==FALSE)])
}

class(Neo_trees) <- "multiPhylo"
class(AU_trees) <- "multiPhylo"
class(AF_trees) <- "multiPhylo"


global_ltt<-multiDiv(trees_run, 
                     int.length = 0.3, 
                     plotLogRich=TRUE, 
                     plot = FALSE)
Neo_ltt<-multiDiv(Neo_trees, 
                  int.length = 0.3, 
                  plotMultCurves = FALSE, 
                  divPalette = cols,
                  plotLogRich=TRUE, 
                  plot = FALSE)
AU_ltt<-multiDiv(AU_trees, 
                 int.length = 0.3, 
                 plotMultCurves = FALSE, 
                 divPalette = cols,
                 plotLogRich=TRUE, 
                 plot = FALSE)
AF_ltt<-multiDiv(AF_trees, 
                 int.length = 0.3, 
                 plotMultCurves = FALSE, 
                 divPalette = cols,
                 plotLogRich=TRUE, 
                 plot = FALSE)

## Color palette
col_regions<-wes_palette("Darjeeling1", 7,  type = "continuous")
col_regions<-col_regions[c(6,5,7,4,2,1,3)]
col_regions2 <- unlist(lapply(col_regions, function(x)t_col(x, percent = 40)))


## Plotting PNG
png("./Figures/Ltt_per_continent.png")
plot(1, type="n", xlab="Time before present (Myr)", ylab="Log(lineages)", 
     xlim=c(15.7,0),ylim=c(0, 7))

plot_ltt(global_ltt, col_curve="grey")
plot_ltt(Neo_ltt, col_curve=col_regions2[3])
plot_ltt(AF_ltt, col_curve=col_regions2[7])
plot_ltt(AU_ltt, col_curve=col_regions2[1])

legend(16, 7, legend = c("Global", "Neotropics", "Australia","Africa"), 
       lwd=c(3,3,3,3),
       col=c("grey",col_regions2[c(3,1,7)]))
dev.off()


## Plotting PDF
pdf("./Figures/Ltt_per_continent.pdf")
plot(1, type="n", xlab="Time before present (Ma)", ylab="Log(lineages)", 
     xlim=c(15.7,0),ylim=c(0, 7), cex.lab = 1.5, cex.axis = 1.5)

plot_ltt(global_ltt, col_curve="grey")
plot_ltt(Neo_ltt, col_curve=col_regions2[3])
plot_ltt(AF_ltt, col_curve=col_regions2[7])
plot_ltt(AU_ltt, col_curve=col_regions2[1])

legend(16, 7, legend = c("Global", "Neotropics", "Australia","Africa"), 
       lwd=c(3,3,3,3),
       col=c("grey",col_regions2[c(3,1,7)]), 
       cex = 1.5)
dev.off()
