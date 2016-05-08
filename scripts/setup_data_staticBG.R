#gbif data procesing: how many species of hummingbird have enough location reports to use for niche modeling? 
library(data.table);library(dismo);library(ecospat);library(SDMTools);library(pbapply);library(foreach);library(doMC);library(rgeos);library(plyr)
setwd("~/Dropbox/nicheTracker_ec2/")

#set extent
ext <- extent(-155,-30,-65,65)

###############################################################
##### process occurrence data
print("loading occurrence data")
gbif <- fread("./occurrence/gbif_icteridae.csv")
gbif <- gbif[,.(species,infraspecificepithet,countrycode,decimallongitude,decimallatitude,day,month,year,collectioncode,institutioncode)]
gbif$species <- factor(gbif$species)
gbif <- subset(gbif,species != "" & is.na(gbif$decimallatitude)==F & is.na(gbif$decimallongitude)==F)
a <- ddply(gbif,"species",.fun=function(e) nrow(subset(e,month %in% c(12,1,2)==T)))
b <- ddply(gbif,"species",.fun=function(e) nrow(subset(e,month %in% c(6,7,8)==T)))
c <- merge(a,b,by="species")
d <- subset(c,c$V1.x > 5 & c$V1.y > 5)
gbif <- subset(gbif,species %in% d$species)
gbif$species <- factor(gbif$species)

###############################################################
###### load climate data.
print("loading climate data")
setwd("./climate/")
folders <- list.files()

#stack rasters, crop to extent, set layer names to month numbers
climateSummary <- function(e){
  setwd(paste("./",e,sep=""))
  files <- list.files()
  files <- files[grep(".bil",files)]
  crop(stack(files),ext)
  #stack(files)
}
climate <- foreach(i=folders) %dopar% climateSummary(i)
names(climate[[1]]) <- c("January","October","November","December","February","March","April","May","June","July","August","September")
names(climate[[2]]) <- c("January","October","November","December","February","March","April","May","June","July","August","September")
names(climate[[3]]) <- c("January","October","November","December","February","March","April","May","June","July","August","September")
names(climate[[4]]) <- c("January","October","November","December","February","March","April","May","June","July","August","September")
names(climate) <- c("prec","tmax","tmean","tmin")

print("summarizing over winter/summer months")
#mean summer values
clim.sum <- foreach(i=climate) %dopar% mean(i$August,i$June,i$July)
names(clim.sum) <- c("prec","tmax","tmean","tmin")

#mean winter values
clim.wnt <- foreach(i=climate) %dopar% mean(i$February,i$December,i$January)
names(clim.wnt) <- c("prec","tmax","tmean","tmin")

setwd("~/Dropbox/nicheTracker_ec2/")

# load altitude 
alt <- crop(raster("./altitude_2.5min/alt.bil"),ext)

#############################################################
###### load ndvi
#See NDVIDataSetup.R for processing biweekly GIMMS data
print("loading ndvi")
ndvi <- stack("./ndvi_monthly.tif")
ndvi.sum <- resample(mean(ndvi$ndvi_monthly.8, ndvi$ndvi_monthly.6, ndvi$ndvi_monthly.7),alt)
ndvi.wnt <- resample(mean(ndvi$ndvi_monthly.2,ndvi$ndvi_monthly.12,ndvi$ndvi_monthly.1),alt)

#############################################################
###### load range shapefiles
print("loading range maps")
setwd("~/Dropbox/nicheTracker_ec2/ranges/Icteridae_rangeMaps/")
files <- list.files()
files <- files[grep(".shp",files)]
ranges <- foreach(i=1:length(files)) %dopar% shapefile(files[i])
ranges.names <- lapply(ranges, function(e) e@data$SCINAME[1])
names(ranges) <- ranges.names

r <- raster(xmn=-180,xmx=180,ymn=-90,ymx=90,res=0.5,crs=proj4string(alt))
simpleRanger <- function(x) {
  a <- rasterize(x,r,field="SEASONAL",background=NA)
  a <- rasterToPolygons(a,dissolve=T)
}
simple.ranges <- foreach(i=ranges) %dopar% simpleRanger(i)
names(simple.ranges) <- names(ranges)

# #make simple buffered ranges
# print("buffering species ranges to available backgorund data")
# r <- raster(xmn=-180,xmx=180,ymn=-90,ymx=90,res=0.5,crs=proj4string(alt))
# simpleRangeBuffer <- function(e) { #simplify and then buffer rangemaps 
#   a <- rasterize(e,r,field="SEASONAL",background=NA)
#   if(all(is.na(values(a)))){ #very small ranges don't fill any raster cells - convert to point and take buffer for faster processing
#     a <- SpatialPoints(e)
#     buffer(a,width=3e6)
#   } else{
#     a <- buffer(a,width=3e6,doEdge=T)
#     rasterToPolygons(a,dissolve=T)
#   }
# }
# ranges.buffered <- foreach(i=ranges) %dopar% simpleRangeBuffer(i)
# ranges.buffered <- pblapply(ranges.buffered,function(e) crop(e,ext))

#filter for species w/mapped ranges
gbif <- subset(gbif,species %in% names(ranges))
gbif$species <- factor(gbif$species)

setwd("~/Dropbox/nicheTracker_ec2/")

# #############################################################
# ####### load landcover
# print("loading landcover")
# lc <- crop(raster("./landcover/AVHRR_1km_LANDCOVER_1981_1994.GLOBAL.tif"),ext)
# m.forest <- matrix(c(1,5,1,5,15,0),byrow=T,nrow = 2,ncol=3) #set matrices for reclassifying landcover
# m.woodland <- matrix(c(1,5,0,6,7,1,8,15,0),byrow=T,nrow = 3,ncol=3)
# m.shrub <- matrix(c(1,7,0,8,9,1,10,15,0),byrow=T,nrow=3,ncol=3)
# lc.forest <- reclassify(lc,m.forest,right=NA) #reclassify to convert to %cover
# lc.woodland <- reclassify(lc,m.woodland,right=NA)
# lc.shrub <- reclassify(lc,m.shrub,right=NA)
# print("still loading landcover")
# lc.forest <- resample(lc.forest,alt) #resample to 2.5min resolution
# lc.woodland <- resample(lc.woodland,alt)
# lc.shrub <- resample(lc.shrub,alt)
# lc.all <- crop(stack(lc.forest,lc.woodland,lc.shrub),ext)#stack and crop to working extent
# names(lc.all) <- c("pc.forest","pc.woodland","pc.shrub")

###########################################################
########## output seasonal background data as raster stack:
# bg.sum.r <- stack(clim.sum$prec,clim.sum$tmean,clim.sum$tmax,clim.sum$tmin,ndvi.sum,lc.all$pc.forest,lc.all$pc.woodland,lc.all$pc.shrub,alt)
# names(bg.sum.r) <- c("prec","tmax","tmin","ndvi","forest","woodland","shrub","altitude")
# bg.wnt.r <- stack(clim.wnt$prec,clim.sum$tmean,clim.wnt$tmax,clim.wnt$tmin,ndvi.wnt,lc.all$pc.forest,lc.all$pc.woodland,lc.all$pc.shrub,alt)
# names(bg.wnt.r) <- c("prec","tmax","tmin","ndvi","forest","woodland","shrub","altitude")

#run below for no landcover
bg.sum.r <- stack(clim.sum$prec,clim.sum$tmean,clim.sum$tmax,clim.sum$tmin,ndvi.sum,alt)
names(bg.sum.r) <- c("prec","tmean","tmax","tmin","ndvi","altitude")
bg.wnt.r <- stack(clim.wnt$prec,clim.wnt$tmean,clim.wnt$tmax,clim.wnt$tmin,ndvi.wnt,alt)
names(bg.wnt.r) <- c("prec","tmean","tmax","tmin","ndvi","altitude")

bg.sum.df <- na.omit(as.data.frame(bg.sum.r))
bg.wnt.df <- na.omit(as.data.frame(bg.wnt.r))

rm(list=ls()[which(ls() %in% c("alt","bg.sum.r","bg.wnt.r","bg.sum.df","bg.wnt.df","gbif","ranges","ranges.buffered","simple.ranges")==F)])
setwd("~/Dropbox/nicheTracker_ec2/")





