##Libraries
library(letsR)
library(maptools)
library(wesanderson)# Check https://github.com/karthik/wesanderson#wes-anderson-palettes for palettes references
library(RColorBrewer)
library(magrittr)


##Data ---
Solanum_geo<-read.csv("./Data/Solanum_BotRec_postClean.csv", fileEncoding="latin1")

### Dropping cultivated species ---
cultivated<-read.csv("./Data/Cultivated_species.csv")
cultivated<-gsub("S.", "Solanum", cultivated$Species)
cultivated<-gsub(" ", "_", cultivated)
Solanum_geo<-Solanum_geo[which(Solanum_geo$species%in%cultivated==FALSE),]

##Shapefiles
world<-readShapePoly("./Data/TM_WORLD_BORDERS-0.3/TM_WORLD_BORDERS-0.3.shp")
all_but_antartica<-world[which(world$NAME!="Antarctica"),]
spl_tmp<- list("sp.lines", as(all_but_antartica, "SpatialLines"), col="dimgrey")

Australia<-world[which(world$NAME=="Australia"),]
sp_australia<- list("sp.lines", as(Australia, "SpatialLines"), col="black")

tipRates<-read.csv("./Output/Avg_rates_tips.csv")
################################################################################
## Calculate the mean of the div rate across all the runs for each species  ----
################################################################################
tipRates$DivMean<-
  tipRates%>%
  within(.,rm("X"))%>%
  rowMeans(.)


## Adding diversification rates to the geographical data ---
Solanum_geo$Div<-tipRates$DivMean[match(Solanum_geo$species, tipRates$X)]
Solanum_geo<-na.omit(Solanum_geo) ## 62304 records  out of 65499 with Div rate.  1000 species with geographical info 

xy<-data.frame(x=Solanum_geo$longdec, y=Solanum_geo$latdec)


# Presence/abscence raster of Solanum ---
Solanum_grids<-lets.presab.points(xy, Solanum_geo$species, resol=1)
Solanum_grids$Richness_Raster[Solanum_grids$Richness_Raster==0]<-NA

#****************************************
# Mapping Diversification rates data ----
#****************************************
source("./Functions/lets_maplizer_mod.R")
Div_map<-lets.maplizer.mod(Solanum_grids,
                           Solanum_geo$Div, 
                           Solanum_geo$species, ras = TRUE)

##Set again the colors
colors=wes_palette(18, name = "Zissou", type = "continuous")
par(bg="transparent")

##General map ---
e <- extent(-170, 180, -60, 75)
p1<-spplot(crop(Div_map$Raster,e),
           sp.layout=spl_tmp,
           col.regions=colors,  
           colorkey=list(labels=list(cex=2)))


## Just Australia
e <- extent(109.90,159.10, -48.75, -10.05)
p1_aus<-spplot(crop(Div_map$Raster, e),sp.layout=sp_australia, colorkey=TRUE,col.regions=colors)

png("./Figures/Mean_diversification_rates_australia_RevBayes.png", width = 1000, height = 600)
par(mar=c(0,0,0,0))
p1_aus
dev.off()

png("./Figures/Mean_diversification_rates_RevBayes.png", width = 700, height = 300)
par(mar=c(0,0,0,0))
p1
dev.off()


##******************
## Weighted values
##******************
range_size<-lets.rangesize(Solanum_grids, units="squaremeter")## Using Sq metres is better than cell size! 

Solanum_geo$range<-range_size[match(Solanum_geo$species,rownames(range_size))]
Solanum_geo$logRange<-log(Solanum_geo$range+1)
Solanum_geo$sqrRange<-sqrt(Solanum_geo$range)
Solanum_geo$powerRange<-(Solanum_geo$range)^(2/3)


Solanum_geo$W_div<-as.vector(Solanum_geo$Div)/Solanum_geo$logRange

Div_weighted<- lets.maplizer.mod(Solanum_grids,
                                 Solanum_geo$W_div,
                                 Solanum_geo$species,
                                 ras = TRUE)

Range_size_map<- lets.maplizer.mod(Solanum_grids, 
                                   Solanum_geo$sqrRange,
                                   Solanum_geo$species,
                                   ras = TRUE)

e <- extent(-170, 180, -60, 75)
Div_weighted_map<-spplot(crop(Div_weighted$Raster,e),
                         sp.layout=spl_tmp,
                         col.regions=colors,
                         colorkey=list(labels=list(cex=2)))


Range_size_map_p<-spplot(Range_size_map$Raster,
                         sp.layout=spl_tmp,
                         col.regions=colors)

png("./Figures/Weighted_Mean_diversification_rates_revBayes.png", width = 700, height = 300)
par(mar=c(0,0,0,0))
Div_weighted_map
dev.off()

e <- extent(109.90,159.10, -48.75, -10.05)
p1_aus_we<-spplot(crop(Div_weighted$Raster, e),
                  sp.layout=sp_australia, 
                  col.regions=colors,
                  colorkey=list(labels=list(cex=2)))

png("./Figures/Weighted_Mean_div_Australia_revBayes.png", width = 700, height = 500)
par(mar=c(0,0,0,0))
p1_aus_we
dev.off()
