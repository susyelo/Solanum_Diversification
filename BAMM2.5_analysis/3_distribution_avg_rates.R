##Libraries
library(letsR)
library(maptools)
library(wesanderson)# Check https://github.com/karthik/wesanderson#wes-anderson-palettes for palettes references
library(RColorBrewer)
library(gridExtra)
library(lattice)

##Data ---
Solanum_geo<-read.csv("./Data/Solanum_BotRec_postClean.csv")

### Dropping cultivated species ---
cultivated<-read.csv("./Data/Cultivated_species.csv")
cultivated<-gsub("S.", "Solanum", cultivated$Species)
cultivated<-gsub(" ", "_", cultivated)
Solanum_geo<-Solanum_geo[which(Solanum_geo$species%in%cultivated==FALSE),]

##Shapefiles
world<-readShapePoly("./Data/TM_WORLD_BORDERS-0.3/TM_WORLD_BORDERS-0.3.shp")
Australia<-world[which(world$NAME=="Australia"),]
all_but_antartica<-world[which(world$NAME!="Antarctica"),]

spl<- list("sp.lines", as(world, "SpatialLines"), col="darkgray")
spl_tmp<- list("sp.lines", as(all_but_antartica, "SpatialLines"), col="dimgrey")
sp_australia<- list("sp.lines", as(Australia, "SpatialLines"), col="black")


#****************************************
# Adding  Diversification rates data ----
#****************************************
tipRates<-read.csv("./Output/Avg_TipsRates.csv")

## Adding diversification rates to the geographical data ---
Solanum_geo$Div<-tipRates$mean[match(Solanum_geo$species, tipRates$X)]
Solanum_geo<-na.omit(Solanum_geo) ## 62304 records  out of 65499 with Div rate.  1000 species with geographical info

xy<-data.frame(x=Solanum_geo$longdec, y=Solanum_geo$latdec)

# Presence/abscence raster of Solanum ---
Solanum_grids<-lets.presab.points(xy, Solanum_geo$species, resol=1)
Solanum_grids$Richness_Raster[Solanum_grids$Richness_Raster==0]<-NA

##Setting the color palette
colors=wes_palette("Zissou1",100, type=c("continuous"))
colors=c("#949494","#949494","#949494","#f0c020", "#ed541d", "#cd0a0a")
colors<-colorRampPalette(colors)(100)


## Adding some breaks to the richness map
plot(log(Solanum_grids$Richness_Raster), 
     breaks = c(0:5),
     col =c("#949494",wes_palette("Zissou1",5)[-1]))


## Richness map ----
colors=wes_palette(14, name = "Zissou1", type = "continuous")

e <- extent(-170, 180, -60, 75)
ras<-spplot(crop(log(Solanum_grids$Richness_Raster+1),e),
       sp.layout=spl_tmp,
       col.regions=colors, 
       colorkey=list(labels=list(cex=1),space = "bottom"),
       at = seq(0.5,5,0.35),
       main=expression(paste("Species richness ","(log"["e"],")")))



## Species richness map of Solanum ----
png("./Figures/Species_richness_Solanum_log.png", width = 700, height = 300)
par(mar=c(0,0,0,0),bg="transparent")
ras
dev.off()



#****************************************
# Mapping Diversification rates data ----
#****************************************

source("./Functions/lets_maplizer_mod.R")
Div_map<-lets.maplizer.mod(Solanum_grids,
                           Solanum_geo$Div,
                           Solanum_geo$species, ras = TRUE)

##Set again the colors
colors=wes_palette(14, name = "Zissou1", type = "continuous")
par(bg="transparent")

##General map ---
e <- extent(-170, 180, -60, 75)
p1<-spplot(crop(Div_map$Raster,e),
           sp.layout=spl_tmp,
           col.regions=colors,
           colorkey=list(labels=list(cex=2)),
           key.space=list(x=0.2,y=0.9,corner=c(0,1)),
           scales=list(draw=T))


Div_map_final<-spplot(crop(Div_map$Raster,e),
                         sp.layout=spl_tmp,
                         col.regions=colors,
                         colorkey=list(labels=list(cex=1),space = "bottom"),
                         at = seq(0.2,0.75,0.04),
                         main=expression(paste("Mean diversification rates ","(species ","Myr"^"-1",")")))


## Just Australia
e <- extent(109.90,159.10, -48.75, -10.05)
p1_aus<-spplot(crop(Div_map$Raster, e),
               sp.layout=sp_australia, 
               colorkey=TRUE,
               at = seq(0.2,0.75,0.04),
               col.regions=colors)

png("./Figures/Mean_diversification_rates_australia.png", width = 1000, height = 600)
par(mar=c(0,0,0,0))
p1_aus
dev.off()

png("./Figures/Mean_diversification_rates.png", width = 700, height = 300)
par(mar=c(0,0,0,0))
p1
dev.off()

##******************
## Weighted values
##******************
range_size<-lets.rangesize(Solanum_grids, units="squaremeter")## Using Sq metres is better than cell size!
range_size<-range_size*0.001


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
                         colorkey=list(labels=list(cex=1),space = "bottom"),
                         at = seq(0.01,0.045,0.004),
                         main=expression(paste("Weighted mean diversification rates ","(species ","Myr"^"-1","/","log"["e"],"m"^"2",")")))


Range_size_map_p<-spplot(Range_size_map$Raster,
                      sp.layout=spl,
                      col.regions=colors)

png("./Figures/Weighted_Mean_diversification_rates.png", width = 700, height = 300)
par(mar=c(0,0,0,0))
Div_weighted_map
dev.off()


## Plotting all maps at the same time ----
layout.height <- trellis.par.get("layout.heights") 
trellis.par.set(layout.heights = layout.height)
png("./Figures/Richness_div.png")
ras
dev.off()

layout.height <- trellis.par.get("layout.heights") 
trellis.par.set(layout.heights = layout.height)
pdf("./Figures/Richness_div.pdf")
ras
dev.off()


layout.height <- trellis.par.get("layout.heights") 
trellis.par.set(layout.heights = layout.height)
png("./Figures/Div_rates_maps.png")
grid.arrange(Div_map_final,
             Div_weighted_map,
             ncol=1,
             nrow=2)
dev.off()

layout.height <- trellis.par.get("layout.heights") 
trellis.par.set(layout.heights = layout.height)
pdf("./Figures/Div_rates_maps.pdf")
grid.arrange(Div_map_final,
             Div_weighted_map,
             ncol=1,
             nrow=2)
dev.off()


trellis.par.set(layout.heights = layout.height)
pdf("./Figures/Richness_div_rates.pdf")
grid.arrange(ras,
             Div_map_final,
             Div_weighted_map,
             ncol=1,
             nrow=3)
dev.off()


e <- extent(109.90,159.10, -48.75, -10.05)
p1_aus_we<-spplot(crop(Div_weighted$Raster, e),
                  sp.layout=list(sp_australia),
                  col.regions=colors,
                  colorkey=list(labels=list(cex=1.5),space = "bottom"),
                  #key.space=list(x=0.2,y=0.9,corner=c(0,1)),
                  scales=list(draw=T,cex=1.5), 
                  at = seq(0.01,0.045,0.004))

png("./Figures/Weighted_Mean_div_Australia.png", width = 900, height = 500)
par(mar=c(0,0,0,0))
p1_aus_we
dev.off()

pdf("./Figures/Weighted_Mean_div_Australia.pdf")
par(mar=c(0,0,0,0))
p1_aus_we
dev.off()


png("./Figures/Range_size_map_sqrt.png", width = 1000, height = 600)
par(mar=c(0,0,0,0))
Range_size_map_p
dev.off()

## Just old world
OW<-geo_total[which(geo_total$Clade=="old_world"),]
OW_xy<-data.frame(x=OW$LONG, y=OW$LAT)

# Presence/abscence matrix ---
OW_grids<-lets.presab.points(OW_xy, OW$species, resol=1)

OW_Div_map<-lets.maplizer.mod(OW_grids, OW$Div, OW$species, ras = TRUE)
Ow_div<-spplot(OW_Div_map$Raster,sp.layout=spl,col.regions=colors)

png("./Figures/Old_world_Div.png", width = 1000, height = 600)
par(mar=c(0,0,0,0))
Ow_div
dev.off()

## Just old world. The weighted version
OW_Div_map_W<-lets.maplizer.mod(OW_grids, OW$W_div, OW$species, ras = TRUE)
Ow_div_W<-spplot(OW_Div_map_W$Raster,sp.layout=spl,col.regions=colors)

png("./Figures/Old_world_Div_W.png", width = 1000, height = 600)
par(mar=c(0,0,0,0))
Ow_div_W
dev.off()
