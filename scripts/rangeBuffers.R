#experimental range buffering: want new set of .shp for each species buffered by X km
#NOTES: added line 15 to project ranges

#read in range shapefiles
print("loading range maps")
setwd("~/Documents/Vireonidae_rangeMaps/Vireonidae/")
files <- list.files()
files <- files[grep(".shp",files)]
names <- strsplit(files,"_")
names <- lapply(names,function(e) e[c(1,2)])
names <- lapply(names,function(e) paste(e[1],e[2]))
names(files) <- names
#files <- files[which(names(files) %in% as.character(gbif$species) ==T)] #209 range shapefiles with names matching gbif species
ranges <- foreach(i=1:length(files)) %dopar% shapefile(files[i])
ranges <- pblapply(ranges,function(e) spTransform(e,CRS("+init=epsg:3395")))
ranges.buffered <- pblapply(ranges,function(e) gSimplify(e,1800))
ranges.buffered <- foreach(i=ranges.buffered) %dopar% buffer(i,width=3e6,dissolve=T)
ranges.names <- lapply(ranges, function(e) e@data$SCINAME[1])
names(ranges) <- ranges.names

r1 <- ranges[[1]]
r2 <- gSimplify(r1,1800)
#testing time for buffering ranges
t1 <- Sys.time()
pblapply(ranges,function(e) buffer(e,width=3e6,dissolve=T))
Sys.time()-t1
t2 <- Sys.time()
foreach(i=ranges) %dopar% buffer(i,width=3e6,dissolve=T)
Sys.time()-t2


