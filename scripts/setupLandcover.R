##adding landcover data##
lc <- crop(raster("~/Documents/worldclim/LandCover/AVHRR_1km_LANDCOVER_1981_1994.GLOBAL.tif"),c(-170,-80,15,70))

m.forest <- matrix(c(1,5,1,5,15,0),byrow=T,nrow = 2,ncol=3)
m.woodland <- matrix(c(1,5,0,6,7,1,8,15,0),byrow=T,nrow = 3,ncol=3)
m.shrub <- matrix(c(1,7,0,8,9,1,10,15,0),byrow=T,nrow=3,ncol=3)

lc.forest <- reclassify(lc,m.forest,right=NA)
lc.woodland <- reclassify(lc,m.woodland,right=NA)
lc.shrub <- reclassify(lc,m.shrub,right=NA)

lc.forest <- resample(lc.forest,clim.winter)
lc.woodland <- resample(lc.woodland,clim.winter)
lc.shrub <- resample(lc.shrub,clim.winter)
lc.all <- stack(lc.forest,lc.woodland,lc.shrub)
names(lc.all) <- c("pc.forest","pc.woodland","pc.shrub")
