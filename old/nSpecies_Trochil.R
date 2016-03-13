##species co-occurrence data, round 2

#code to generate species diversity raster from birdlife international range maps (~30min run time)
setwd("~/Documents/Trochilidae_rangeMaps/")
files <- list.files()
files <- files[grep(".shp",files)]
r <- raster(xmn=-180,xmx=180,ymn=-90,ymx=90,res=1/6)

##########single file with all ranges
# ranges <- lapply(files,FUN=function(e) rasterize(shapefile(e),r,field=1))
# ranges <- crop(stack(ranges),ext)
# n.species <- stackApply(ranges,indices=rep(1,nlayers(ranges)),fun=sum)
# writeRaster(n.species,"~/Documents/nicheTracker/data/Trochilidae_nSpecies.tif",overwrite=T)
# 
# n.species <- raster("~/Documents/nicheTracker/data/Trochilidae_nSpecies.tif")
# setwd("~/Documents/nicheTracker")

###########splitting resident and migratory species. Warning: slow. (20min?)
#all <- lapply(files,function(e) shapefile(e))
#ranges.brd <- lapply(all,FUN=function(e) e[which(e@data$SEASONAL == 2 |e@data$SEASONAL == 1),])
#ranges.wnt <- lapply(all,FUN=function(e) e[which(e@data$SEASONAL == 3 |e@data$SEASONAL == 1),])
# anhu <- all[[75]]
# anhu <- anhu[which(anhu@data$SEASONAL == 1 | anhu@data$SEASONAL == 2),]
# ranges.wnt[[75]] <- anhu
# ranges.brd.r <- lapply(ranges.brd,function(e) rasterize(e,r,field=1))
# ranges.wnt.r <- lapply(ranges.wnt,function(e) rasterize(e,r,field=1))
# nspecies.brd <- stack(ranges.brd.r)
# nspecies.brd <- stackApply(nspecies.brd,indices=rep(1,length(ranges.brd)),fun=sum)
# nspecies.wnt <- list(nspecies.wnt,ranges.wnt.)
# nspecies.wnt <- stack(ranges.wnt.r)
# nspecies.wnt <- stackApply(nspecies.wnt,indices=rep(1,length(ranges.wnt)),fun=sum)
# 
# writeRaster(nspecies.wnt,"nspecies_Trochilidae_nonbreeding.tif",overwrite=T)
# writeRaster(nspecies.brd,"nspecies_Trochilidae_breeding.tif",overwrite=T)

nspecies.wnt <- crop(raster("nspecies_Trochilidae_nonbreeding.tif"),ext)
nspecies.brd <- crop(raster("nspecies_Trochilidae_breeding.tif"),ext)
print("species diversity data loaded. objects: nspecies.wnt & nspecies.brd. format: raster objects")
setwd("~/Documents/nicheTracker/")
