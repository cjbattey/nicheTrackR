##niche trackR##
#Building and cross-predicting ENM's for wintering and breeding season hummingbirds with PCA and Maxent 
library(ecodist);library(ggplot2);library(raster);library(dismo);library(rgeos)
setwd("~/Documents/nicheTracker/")
map <- crop(shapefile("~/Documents/worldclim/country_outlines/cntry06/cntry06.shp"),ext)


###################### Multispecies Read-in and processing 

print("reading in range maps")
setwd("~/Documents/nicheTracker/data/ebird/")

files <- list.files()
files <- files[grep(".txt",files)]
range.files <- list.files("~/Documents/Trochilidae_rangeMaps/")
range.files <- range.files[grep(".shp",range.files)]
range.file.index <- c(338,75,57,336,60,334,76,73,58,337,339,335)
range.files <- range.files[range.file.index]
ranges <- lapply(range.files,FUN=function(e) shapefile(paste("~/Documents/Trochilidae_rangeMaps/",e,sep="")))
ranges.winter <- lapply(ranges,function(e) e[which(e@data$SEASONAL == 1 | e@data$SEASONAL == 3),])
ranges.breeding <-  lapply(ranges,function(e) e[which(e@data$SEASONAL == 1 | e@data$SEASONAL == 2),])

brd.months <- list(alhu=c(3,4,5),anhu=c(12,1,2),bchu=c(5,6),bthu=c(5,6,7),buhu=c(6,7),cahu=c(6,7),cohu=c(2:4),
                   luhu=c(6,7),rthu=c(6,7,8),ruhu=c(5,6),schu=c(6,7),vohu=c(6,7))
wnt.months <- list(alhu=c(10,11,12),anhu=c(8,9,10),bchu=c(1,2),bthu=c(11,12,1),buhu=c(12,1),cahu=c(12,1),cohu=c(2:4),
                   luhu=c(12,1),rthu=c(12,1,2),ruhu=c(12,1),schu=c(12,1),vohu=c(12,1))

print("reading in ebird reports")
occ.wnt <- list()
occ.brd <- list()
for (i in 1:12){
  a <- read.delim(files[i],header=T,quote="")
  a$yday <- as.POSIXlt(as.character(a$OBSERVATION.DATE),"%Y-%m-%d",tz="GMT")$yday
  a$month <- as.numeric(substr(a$OBSERVATION.DATE,6,7))
  a.brd <- a[which(a$month %in% brd.months[[i]]==T),]
  a.brd <- SpatialPoints(data.frame(a.brd$LONGITUDE,a.brd$LATITUDE))
  a.brd <- gIntersection(a.brd,ranges.breeding[[i]])
  if(i %in% c(1,10) == T){
    a <- subset(a,a$COUNTRY== "Mexico" | a$STATE_PROVINCE == "California")
  }
  a.wnt <- a[which(a$month %in% wnt.months[[i]]==T),]
  a.wnt <- SpatialPoints(data.frame(a.wnt$LONGITUDE,a.wnt$LATITUDE))
  if (i != 1){
  a.wnt <- gIntersection(a.wnt,ranges.winter[[i]])
  }
  plot(a.brd)+plot(map,add=T)
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

plot(occ.brd$luhu)+plot(map,add=T)
# #Vertnet files (no taxonomy fields... why?)
# Selasphorus <- read.delim("data/vertnet/Selasphorus_map-083d2e517b2b4374aa8a9815051c44d9.tsv",header=T)
# Aheloisa <- read.delim("data/vertnet/Atthis_heloisa.txt")
# Aellioti <- read.delim("data/vertnet/Atthis_ellioti.txt")
# Clucifer <- read.delim("data/vertnet/Calothorax_lucifer.txt")
# Calypte <- read.delim("data/vertnet/Calypte_map-2c2153f96bf549cc8fa884ec0c5e6b2c.tsv")
# Archilochus <- read.delim("data/vertnet/Archilochus_map-da613a5dc40749d4916563afee9eb908.tsv")
# species <- rbind(Selasphorus,Aheloisa,Aellioti,Clucifer,Calypte,Archilochus)
# 

# ### S. rufus only ###
# print("reading in occurrence data for S. rufus")
# #read in rufous hummingbird data
# ruhu <- read.delim("~/Documents/ruhu/sampling/ebd_rufhum_relNov-2014/ebd_rufhum_relNov-2014.txt")
# ruhu$yday <- as.POSIXlt(as.character(ruhu$OBSERVATION.DATE),"%m/%d/%Y",tz="GMT")$yday
# ruhu$date <- as.Date(ruhu$OBSERVATION.DATE,"%m/%d/%Y")
# ruhu$month <- as.numeric(substr(ruhu$date,6,7))
# ruhu$monthName <- months(ruhu$date)
# #subset out birds from expanded wintering range in SE US
# ruhuW <- subset(ruhu,ruhu$LONGITUDE < -96 & ruhu$LATITUDE < 25)
# ruhu.range <- shapefile("~/Documents/Trochilidae_rangeMaps/Selasphorus_rufus_22688296.shp")
# ruhu.wint.range <- ruhu.range[which(ruhu.range@data$SEASONAL == 3),]
# ruhu.brd.range <- ruhu.range[which(ruhu.range@data$SEASONAL == 2),]
# 
# print("subsetting to breeding and wintering ranges")
# #subset to wintering and breeding months and filter by natureserve range maps
# ruhu.wint <- ruhuW[which(ruhuW$month %in% c(11,12,1) == T),]
# ruhu.wint <- SpatialPoints(data.frame(ruhu.wint$LONGITUDE,ruhu.wint$LATITUDE))
# ruhu.wint <- gIntersection(ruhu.wint,ruhu.wint.range)         
# 
# ruhu.brd <- ruhu[which(ruhu$month %in% c(4,5,6) == T),]
# ruhu.brd <- SpatialPoints(data.frame(ruhu.brd$LONGITUDE,ruhu.brd$LATITUDE))
# ruhu.brd <- gIntersection(ruhu.brd,ruhu.brd.range)
# 
# print("Occurrence data loaded. Format: SpatialPoints. Objects: ruhu.wint & ruhu.brd")
