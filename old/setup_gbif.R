#gbif data procesing: how many species of hummingbird have enough location reports to use for niche modeling? 
library(data.table);library(dismo);library(ecospat);library(SDMTools);library(pbapply)
setwd("~/Documents/nicheTracker/")

##########################################################
################### load occurrence data #################
##########################################################
#load gbif data
gbif <- fread("data/gbif_trochilidae.csv")
gbif <- gbif[,.(species,infraspecificepithet,countrycode,decimallongitude,decimallatitude,day,month,year,collectioncode,institutioncode)]
gbif$species <- factor(gbif$species)

#set default winter/summer months (replace with breeding & wintering months when possible)
s <- c(5,6,7)
w <- c(11,12,1)

#checking species w/wo 50 reports in summer/winter
number.reports <- data.frame(species=NA,summer.reports=NA,winter.reports=NA)
for(i in levels(gbif$species)){
    summer.reports <- nrow(gbif[which(gbif$species == i & gbif$month %in% s == T)])
    winter.reports <- nrow(gbif[which(gbif$species == i & gbif$month %in% w == T)])
    number.reports <- rbind(number.reports,data.frame(species=i,summer.reports=summer.reports,winter.reports=winter.reports))
}
y <- number.reports[which(number.reports$summer.reports > 50 & number.reports$winter.reports > 50),]#222 species with > 50 reports in both seasons

#filter for species with sufficient reports
gbif <- gbif[which(gbif$species %in% y$species ==T)]
gbif$species <- factor(gbif$species)

print("occurrence data loaded")

#########################################################
######################### climate #######################
#########################################################
setwd("~/Documents/worldclim/worldclim_current_monthly_2.5/")
folders <- list.files()

#read in rasters as stacks
climate <- lapply(folders,FUN=function(e){
  setwd(paste("~/Documents/worldclim/worldclim_current_monthly_2.5/",e,sep=""))
  files <- list.files()
  files <- files[grep(".bil",files)]
  stack(files)
})

print("climate data loaded")
setwd("~/Documents/nicheTracker/")

##########################################################
########################## ndvi ########################## See NDVIDataSetup.R for processing biweekly GIMMS data
##########################################################
ndvi <- stack("~/Dropbox/phenology/ndvi_monthly.tif")
print("ndvi loaded")

##########################################################
######################  range shapefiles #################
##########################################################
setwd("~/Documents/Trochilidae_rangeMaps/")
files <- list.files()
files <- files[grep(".shp",files)]
names <- strsplit(files,"_")
names <- lapply(names,function(e) e[c(1,2)])
names <- lapply(names,function(e) paste(e[1],e[2]))
names(files) <- names
files <- files[which(names(files) %in% as.character(gbif$species) ==T)] #209 range shapefiles with names matching gbif species
ranges <- pblapply(files,function(e) shapefile(e))

print("range maps loaded")

##########################################################
#####################  landcover  ########################
##########################################################
lc <- crop(raster("~/Documents/worldclim/LandCover/AVHRR_1km_LANDCOVER_1981_1994.GLOBAL.tif"),c(-170,-80,15,70))

m.forest <- matrix(c(1,5,1,5,15,0),byrow=T,nrow = 2,ncol=3)
m.woodland <- matrix(c(1,5,0,6,7,1,8,15,0),byrow=T,nrow = 3,ncol=3)
m.shrub <- matrix(c(1,7,0,8,9,1,10,15,0),byrow=T,nrow=3,ncol=3)

lc.forest <- reclassify(lc,m.forest,right=NA)
lc.woodland <- reclassify(lc,m.woodland,right=NA)
lc.shrub <- reclassify(lc,m.shrub,right=NA)

lc.forest <- resample(lc.forest,r)
lc.woodland <- resample(lc.woodland,r)
lc.shrub <- resample(lc.shrub,r)
lc.all <- crop(stack(lc.forest,lc.woodland,lc.shrub),ext)
names(lc.all) <- c("pc.forest","pc.woodland","pc.shrub")



