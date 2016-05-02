#gbif data procesing: how many species of hummingbird have enough location reports to use for niche modeling? 
library(data.table);library(dismo);library(ecospat);library(SDMTools);library(pbapply);library(foreach);library(doMC);library(rgeos)

#set extent
ext <- extent(-155,-30,-65,65)

###############################################################
##### process occurrence data
print("loading occurrence data")
gbif <- fread("/Volumes/Seagate Backup Plus Drive/nicheTracker_data/gbif_tyrannidae.csv")
gbif <- gbif[,.(species,infraspecificepithet,countrycode,decimallongitude,decimallatitude,day,month,year,collectioncode,institutioncode)]
gbif$species <- factor(gbif$species)
gbif <- subset(gbif,species != "" & is.na(gbif$decimallatitude)==F & is.na(gbif$decimallongitude)==F)
gbif$species <- factor(gbif$species)

###############################################################
###### load climate data.
print("loading climate data")
setwd("~/Documents/worldclim/worldclim_current_monthly_2.5/")
folders <- list.files()

#stack rasters, crop to extent, set layer names to month numbers
climateSummary <- function(e){
  setwd(paste("~/Documents/worldclim/worldclim_current_monthly_2.5/",e,sep=""))
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
clim.sum <- foreach(i=climate) %dopar% mean(i$May,i$June,i$July)
names(clim.sum) <- c("prec","tmax","tmean","tmin")

#mean winter values
clim.wnt <- foreach(i=climate) %dopar% mean(i$November,i$December,i$January)
names(clim.wnt) <- c("prec","tmax","tmean","tmin")

setwd("/R/nicheTracker/")

# load altitude 
alt <- crop(raster("~/Documents/worldclim/alt_2-5m_bil/alt.bil"),ext)

#############################################################
###### load ndvi
#See NDVIDataSetup.R for processing biweekly GIMMS data
print("loading ndvi")
ndvi <- stack("~/Dropbox/gimms/ndvi_monthly.tif")
ndvi.sum <- resample(mean(ndvi$ndvi_monthly.5, ndvi$ndvi_monthly.6, ndvi$ndvi_monthly.7),alt)
ndvi.wnt <- resample(mean(ndvi$ndvi_monthly.11,ndvi$ndvi_monthly.12,ndvi$ndvi_monthly.1),alt)

#############################################################
###### load range shapefiles
print("loading range maps")
setwd("~/Documents/Tyrannidae_rangeMaps/")
files <- list.files()
files <- files[grep(".shp",files)]
ranges <- foreach(i=1:length(files)) %dopar% shapefile(files[i])
ranges <- pblapply(ranges,function(e) spTransform(e,CRS("+init=epsg:3395")))
ranges.names <- lapply(ranges, function(e) e@data$SCINAME[1])
names(ranges) <- ranges.names
setwd("/R/nicheTracker/")

#############################################################
####### load landcover
print("loading landcover")
lc <- crop(raster("~/Documents/worldclim/LandCover/AVHRR_1km_LANDCOVER_1981_1994.GLOBAL.tif"),ext)
m.forest <- matrix(c(1,5,1,5,15,0),byrow=T,nrow = 2,ncol=3) #set matrices for reclassifying landcover
m.woodland <- matrix(c(1,5,0,6,7,1,8,15,0),byrow=T,nrow = 3,ncol=3)
m.shrub <- matrix(c(1,7,0,8,9,1,10,15,0),byrow=T,nrow=3,ncol=3)
lc.forest <- reclassify(lc,m.forest,right=NA) #reclassify to convert to %cover
lc.woodland <- reclassify(lc,m.woodland,right=NA)
lc.shrub <- reclassify(lc,m.shrub,right=NA)
print("still loading landcover")
lc.forest <- resample(lc.forest,alt) #resample to 2.5min resolution
lc.woodland <- resample(lc.woodland,alt)
lc.shrub <- resample(lc.shrub,alt)
lc.all <- crop(stack(lc.forest,lc.woodland,lc.shrub),ext)#stack and crop to working extent
names(lc.all) <- c("pc.forest","pc.woodland","pc.shrub")

###########################################################
########## output seasonal background data as raster stack:
bg.sum.r <- stack(clim.sum$prec,clim.sum$tmax,clim.sum$tmin,ndvi.sum,lc.all$pc.forest,lc.all$pc.woodland,lc.all$pc.shrub,alt)
names(bg.sum.r) <- c("prec","tmax","tmin","ndvi","forest","woodland","shrub","altitude")
bg.wnt.r <- stack(clim.wnt$prec,clim.wnt$tmax,clim.wnt$tmin,ndvi.wnt,lc.all$pc.forest,lc.all$pc.woodland,lc.all$pc.shrub,alt)
names(bg.wnt.r) <- c("prec","tmax","tmin","ndvi","forest","woodland","shrub","altitude")

bg.sum.df <- na.omit(as.data.frame(bg.sum.r))
bg.wnt.df <- na.omit(as.data.frame(bg.wnt.r))

rm(list=ls()[which(ls() %in% c("alt","bg.sum.r","bg.wnt.r","bg.sum.df","bg.wnt.df","gbif","ranges","ranges.buffered")==F)])
setwd("/R/nicheTracker/")





