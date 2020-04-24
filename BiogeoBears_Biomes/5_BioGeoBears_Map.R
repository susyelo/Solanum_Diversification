library(maps)
library(geosphere)
library(dplyr)
library(nycflights13)
library(maptools)
library(diagram)
library(ggplot2)
library(rgeos)
library(tidyr)
library(BSDA)
library(fields)

## Fucntions ----

Points2Regions = function(points, regions_shp, abbre=FALSE)
{  
  pointsSP = SpatialPoints(points, proj4string=CRS(proj4string(regions_shp)))  
  # use 'over' to get indices of the Polygons object containing each point 
  indices = over(pointsSP, regions_shp)
  # return the ADMIN names of each country
  
  if (abbre==TRUE){
    return(indices$WWF_REALM)
  }else
  {
    return(indices$WWF_REALM2)
  }
}



## Processes 
## Reading BioGeoBears data --- 

load("./Output/Event_count_list_DECjM1.RData")

d_counts_mean<-counts_list$d_counts_fromto_means
d_counts_sd<-counts_list$d_counts_fromto_sds

rownames(d_counts_mean)<-c("African","Australian","Indo_Pacific","Nearctic","Neotropic","Palearctic")
names(d_counts_mean)<-c("African","Australian","Indo_Pacific","Nearctic","Neotropic","Palearctic")
names(d_counts_sd)<-c("African","Australian","Indo_Pacific","Nearctic","Neotropic","Palearctic")

new_d_counts_mean<-d_counts_mean%>%
  gather(To,dRate_mean)

new_d_counts_sd<-d_counts_sd%>%
  gather(To,dRate_sd)


new_d_counts_mean$From<-rep(rownames(d_counts_mean),6)
new_d_counts_mean$sd<-new_d_counts_sd$dRate_sd
new_d_counts<-new_d_counts_mean[,c("From","To","dRate_mean","sd")]

new_d_counts<-new_d_counts[which(new_d_counts$dRate_mean!=0),]
new_d_counts$p_values<-tsum.test(new_d_counts$dRate_mean,new_d_counts$sd,98)$p.value

new_d_counts<-new_d_counts[which(new_d_counts$dRate_mean>=2&new_d_counts$p_values<0.001),]##values that are significant different from 0
  


Cox_maps<-readShapePoly("~/Box Sync/Susy Echeverria-Londono/Solanaceae/Analysis/Biogeography/Data/shapefiles/Cox_areas/Cox_areas.shp")
centroids <- getSpPPolygonsLabptSlots(Cox_maps)
#plot(Cox_maps)
#points(centroids, pch = 3, col = "Red")

## Creating datafrane with the centroids coordinates
cens_cox<-data.frame(name=Cox_maps$WWF_REALM2,long=centroids[,1],lat=centroids[,2])
edges <- merge(merge(new_d_counts, cens_cox[, c("name", "long", "lat")], by.x = "From", by.y = "name"), cens_cox[, c("name", "long", "lat")], by.x = "To", by.y = "name")
edges<-edges[,c(2,1,3:9)]


library(wesanderson)
cens_cox$colors=wes_palette("Zissou",7, type=c("continuous"))
edges$col <- as.character(cens_cox$colors[match(edges$to, cens_cox$name)])

Cox_maps_sim<-gSimplify(Cox_maps,tol=0.5,topologyPreserve=TRUE)
proj4string(Cox_maps_sim)<-"+proj=longlat +datum=WGS84"
##
plot(Cox_maps_sim, border="grey40", col = "lightgray", mar = rep(0, 4),xlim=c(-180,180),ylim = c(-60,70))

map(Cox_maps_sim, center=180, col="white",bg="gray",
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))


## Also check: http://phylo.wikidot.com/biogeonick
## Geographic data to count the number of species per region
Cox_maps<-readShapePoly("~/Box Sync/Susy Echeverria-Londono/Solanaceae/Analysis/Biogeography/Data/shapefiles/Cox_areas/Cox_areas.shp")
Solanum_geo<-read.csv("./Data/Solanum_BotRec_postClean.csv")##1074 with records
coordinates(Solanum_geo)<-c("longdec","latdec")
Solanum_geo$Cox_area<-Points2Regions(Solanum_geo,Cox_maps)

N_species_per_regions<-tapply(Solanum_geo$species, Solanum_geo$Cox_area, function(x)length(unique(x)))
reg_cols<-data.frame(N_sp=N_species_per_regions[order(N_species_per_regions)])


##Cox_maps_sim order = Australia, Antarctic,Africa,Indopacific, Neartic, Neotropics, Palearti 
reg_cols$cols<-c(rev(grey.colors(6)),"white")

col_breaks = c(seq(min(N_species_per_regions,na.rm=TRUE),max(N_species_per_regions,na.rm=TRUE),length=100))
               

pdf("./Figures/Dispersion_maps_DEC+J_moreThan2.pdf")
par(mar = c(0,0,0,0))

plot(Cox_maps_sim, col=reg_cols$cols[match(Cox_maps$WWF_REALM2,rownames(reg_cols))],border=NA,ann=FALSE, axes = FALSE)

for (i in 1:nrow(edges)){
  point1 = c(edges$long.x[i],edges$lat.x[i])
  point2 = c(edges$long.y[i],edges$lat.y[i])
  
  curvedarrow(from=point1,to=point2,  lcol=wes_palette("Zissou",1), curve=.1, lwd=edges$dRate_mean[i],
              arr.adj = 1, arr.pos = 0.5,arr.type = "triangle",arr.lwd = edges$dRate_mean[i]+1)
}

par(mar = c(5,5,5,5))

image.plot(zlim=c(min(N_species_per_regions,na.rm=T),max(N_species_per_regions,na.rm=T)),nlevel=50, legend.only=TRUE,
           horizontal=T, col=rev(grey.colors(50)),legend.mar =8)
dev.off()




color.bar(lut=rev(grey.colors(100)),min=min(N_species_per_regions,na.rm=TRUE),max=max(N_species_per_regions,na.rm=TRUE))
plot( 1:10, (1:10)*10, type="n", bty="n") 
colorbar.plot( 2.5, 30, col_breaks,col=grey.colors(100),xlim=min(N_species_per_regions,na.rm=TRUE))

image.plot(zlim=c(min(N_species_per_regions,na.rm=T),max(N_species_per_regions,na.rm=T)), nlevel=50, legend.only=TRUE,
           horizontal=F, col=rev(grey.colors(100)),legend.mar = 45)
