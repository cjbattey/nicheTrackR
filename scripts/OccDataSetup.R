##niche trackR##
#Building and cross-predicting ENM's for wintering and breeding season hummingbirds with PCA and Maxent 
library(ecodist);library(ggplot2);library(raster);library(dismo);library(rgeos)

print("reading in occurrence data for S. rufus")
#read in rufous hummingbird data
ruhu <- read.delim("~/Documents/ruhu/sampling/ebd_rufhum_relNov-2014/ebd_rufhum_relNov-2014.txt")
ruhu$yday <- as.POSIXlt(as.character(ruhu$OBSERVATION.DATE),"%m/%d/%Y",tz="GMT")$yday
ruhu$date <- as.Date(ruhu$OBSERVATION.DATE,"%m/%d/%Y")
ruhu$month <- as.numeric(substr(ruhu$date,6,7))
ruhu$monthName <- months(ruhu$date)
#subset out birds from expanded wintering range in SE US
ruhuW <- subset(ruhu,ruhu$LONGITUDE < -96 & ruhu$LATITUDE < 25)
ext <- extent(c(-180,0,0,90))
map <- crop(shapefile("~/Documents/worldclim/country_outlines/cntry06/cntry06.shp"),ext)
ruhu.range <- shapefile("~/Documents/Trochilidae_rangeMaps/Selasphorus_rufus_22688296.shp")
ruhu.wint.range <- ruhu.range[which(ruhu.range@data$SEASONAL == 3),]
ruhu.brd.range <- ruhu.range[which(ruhu.range@data$SEASONAL == 2),]

print("subsetting to breeding and wintering ranges")
#subset to wintering and breeding months and filter by natureserve range maps
ruhu.wint <- ruhuW[which(ruhuW$month %in% c(11,12,1) == T),]
ruhu.wint <- SpatialPoints(data.frame(ruhu.wint$LONGITUDE,ruhu.wint$LATITUDE))
ruhu.wint <- gIntersection(ruhu.wint,ruhu.wint.range)         

ruhu.brd <- ruhu[which(ruhu$month %in% c(4,5,6) == T),]
ruhu.brd <- SpatialPoints(data.frame(ruhu.brd$LONGITUDE,ruhu.brd$LATITUDE))
ruhu.brd <- gIntersection(ruhu.brd,ruhu.brd.range)

print("Occurrence data loaded. Format: SpatialPoints. Objects: ruhu.wint & ruhu.brd")

#vert <- read.delim("~/Downloads/Srufus_MX_vertNet.txt",header=T)



