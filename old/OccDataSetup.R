##niche trackR##
#Building and cross-predicting ENM's for wintering and breeding season hummingbirds with PCA and Maxent 
library(ecodist);library(ggplot2);library(raster);library(dismo);library(rgeos);library(data.table)
setwd("~/Documents/nicheTracker/")
map <- crop(shapefile("~/Documents/worldclim/country_outlines/cntry06/cntry06.shp"),ext)


###################### Multispecies Read-in and processing 

print("reading in range maps")
setwd("~/Documents/nicheTracker/data/ebird/")

range.files <- list.files("~/Documents/Trochilidae_rangeMaps/")
range.files <- range.files[grep(".shp",range.files)]
range.file.index <- c(338,75,57,336,60,334,76,73,58,337,339,335)
range.files <- range.files[range.file.index]
ranges <- lapply(range.files,FUN=function(e) shapefile(paste("~/Documents/Trochilidae_rangeMaps/",e,sep="")))
ranges.winter <- lapply(ranges,function(e) e[which(e@data$SEASONAL == 1 | e@data$SEASONAL == 3),])
ranges.breeding <-  lapply(ranges,function(e) e[which(e@data$SEASONAL == 1 | e@data$SEASONAL == 2),])

print("reading in ebird reports")
occ.wnt <- list()
occ.brd <- list()
files <- list.files()
files <- files[grep(".txt",files)]
for (i in 1:12){
  a <- read.delim(files[i],header=T,quote="")
  a$yday <- as.POSIXlt(as.character(a$OBSERVATION.DATE),"%Y-%m-%d",tz="GMT")$yday
  a$month <- as.numeric(substr(a$OBSERVATION.DATE,6,7))
  a.brd <- a[which(a$month %in% brd.months[[i]]==T),]
  a.brd <- SpatialPoints(data.frame(a.brd$LONGITUDE,a.brd$LATITUDE))
  a.brd <- gIntersection(a.brd,ranges.breeding[[i]])
  a.wnt <- a[which(a$month %in% wnt.months[[i]]==T),]
  if(i %in% c(1,10) == T){
    a.wnt <- subset(a.wnt,a.wnt$COUNTRY== "Mexico" | a.wnt$STATE_PROVINCE == "California")
  }
  a.wnt <- SpatialPoints(data.frame(a.wnt$LONGITUDE,a.wnt$LATITUDE))
  if (i != 1){
  a.wnt <- gIntersection(a.wnt,ranges.winter[[i]])
  }
  plot(map)+points(a.wnt,col="blue")+points(a.brd,col="red")
  occ.wnt <- c(occ.wnt,a.wnt)
  occ.brd <- c(occ.brd,a.brd)
#   assign(paste(substr(files[i],5,9),".wnt",sep=""),a.wnt)
#   assign(paste(substr(files[i],5,9),".brd",sep=""),a.brd)
  print(paste(i,"/12",sep=""))
}
names(occ.wnt) <- names(wnt.months)
names(occ.brd) <- names(brd.months)

occ.wnt <- lapply(occ.wnt,function(e) SpatialPoints(gridSample(e,r,n=1)))
occ.brd <- lapply(occ.brd,function(e) SpatialPoints(gridSample(e,r,n=1)))

setwd("~/Documents/nicheTracker/")
