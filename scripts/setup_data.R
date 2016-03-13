#gbif data procesing: how many species of hummingbird have enough location reports to use for niche modeling? 
library(data.table);library(dismo);library(ecospat);library(SDMTools);library(pbapply)
setwd("~/R/nicheTracker/")

#set extent
ext <- c(-155,-30,-65,65)

#load species breeding and nonbreeding months
#? tbd

###############################################################
##### process occurrence data
print("loading occurrence data")
gbif <- fread("data/gbif_vireonidae.csv")
gbif <- gbif[,.(species,infraspecificepithet,countrycode,decimallongitude,decimallatitude,day,month,year,collectioncode,institutioncode)]
gbif$species <- factor(gbif$species)
gbif <- subset(gbif,species != "")
s <- c(5,6,7)#set default winter/summer months 
w <- c(11,12,1)
number.reports <- data.frame(species=NA,summer.reports=NA,winter.reports=NA)#checking species w/wo 50 reports in summer/winter
j <- 0
for(i in levels(gbif$species)){
    summer.reports <- nrow(gbif[which(gbif$species == i & gbif$month %in% s == T)])
    winter.reports <- nrow(gbif[which(gbif$species == i & gbif$month %in% w == T)])
    number.reports <- rbind(number.reports,data.frame(species=i,summer.reports=summer.reports,winter.reports=winter.reports))
    j <- j+1
    print(paste(j,"/",nlevels(gbif$species),sep=""))
}
y <- number.reports[which(number.reports$summer.reports > 50 & number.reports$winter.reports > 50),]#222 species with > 50 reports in both seasons
gbif <- gbif[which(gbif$species %in% y$species ==T & is.na(gbif$decimallatitude)==F & is.na(gbif$decimallongitude)==F)]#filter for species with sufficient reports
gbif$species <- factor(gbif$species)

###############################################################
###### load climate data.
print("loading climate data")
setwd("~/Documents/worldclim/worldclim_current_monthly_2.5/")
folders <- list.files()

#stack rasters, crop to extent, set layer names to month numbers
climate <- pblapply(folders,FUN=function(e){
  setwd(paste("~/Documents/worldclim/worldclim_current_monthly_2.5/",e,sep=""))
  files <- list.files()
  files <- files[grep(".bil",files)]
  crop(stack(files),ext)
})
names(climate[[1]]) <- c("January","October","November","December","February","March","April","May","June","July","August","September")
names(climate[[2]]) <- c("January","October","November","December","February","March","April","May","June","July","August","September")
names(climate[[3]]) <- c("January","October","November","December","February","March","April","May","June","July","August","September")
names(climate[[4]]) <- c("January","October","November","December","February","March","April","May","June","July","August","September")
names(climate) <- c("prec","tmax","tmean","tmin")

print("summarizing over winter/summer months")
#mean summer values
clim.sum <- pblapply(climate,FUN=function(e){
  mean(e$May,e$June,e$July)
})
names(clim.sum) <- c("prec","tmax","tmean","tmin")

#mean winter values
clim.wnt <- pblapply(climate,FUN=function(e){
  mean(e$November,e$December,e$January)
})
names(clim.wnt) <- c("prec","tmax","tmean","tmin")

setwd("~/R/nicheTracker/")

# load altitude 
alt <- crop(raster("~/Documents/worldclim/alt_2-5m_bil/alt.bil"),ext)

#############################################################
###### load ndvi
#See NDVIDataSetup.R for processing biweekly GIMMS data
print("loading ndvi")
ndvi <- stack("~/Dropbox/phenology/ndvi_monthly.tif")
ndvi.sum <- resample(mean(ndvi$ndvi_monthly.5, ndvi$ndvi_monthly.6, ndvi$ndvi_monthly.7),alt)
ndvi.wnt <- resample(mean(ndvi$ndvi_monthly.11,ndvi$ndvi_monthly.12,ndvi$ndvi_monthly.1),alt)


#############################################################
###### load range shapefiles
print("loading range maps")
setwd("~/Documents/Trochilidae_rangeMaps/")
files <- list.files()
files <- files[grep(".shp",files)]
names <- strsplit(files,"_")
names <- lapply(names,function(e) e[c(1,2)])
names <- lapply(names,function(e) paste(e[1],e[2]))
names(files) <- names
#files <- files[which(names(files) %in% as.character(gbif$species) ==T)] #209 range shapefiles with names matching gbif species
ranges <- pblapply(files,function(e) shapefile(e))
ranges.names <- lapply(ranges, function(e) e@data$SCINAME[1])
names(ranges) <- ranges.names
setwd("~/Documents/nicheTracker/")

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
#&data frame
bg.sum <- na.omit(as.data.frame(bg.sum.r,xy=T)[3:10])
print("summer background ready")
bg.wnt <- na.omit(as.data.frame(bg.wnt.r,xy=T)[3:10]) 
print("winter background ready")









