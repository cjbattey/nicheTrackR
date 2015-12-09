#NDVI processing#
library(gimms)
setwd("~/Dropbox/phenology/gimms/")

#This script downloads all available GIMMS NDVI data from NASA servers, rasterizes the data, then averages within and across months. Output is 
#mean monthly NDVI .tif's for 1980-2013. Requires ~2 hrs to run, 18gb hard drive space, and ~12gb free RAM. 

# gimms <- downloadGimms(1980,2013,overwrite=F)
# 
# files <- list.files()
# files <- files[!grepl("VI3g.",files)]
# 
# #read in pair of files, average, store. Need ~14gb free RAM.
# for(i in seq(1,(length(files)/2),by=2)){
#   a <- rasterizeGimms(files[i:(i+1)])
#   b <- mean(a)
#   print(paste(((i/(length(files)/2))*100),"% complete"))
#   if (i == 1){
#     c <- stack(b)
#   }
#   else {
#     c <- addLayer(c,b)
#   }
# }
# names(c) <- substr(files[seq(1,(length(files)/2),2)],4,8)
# months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
# month.indices <- lapply(months,FUN=function(e) grep(e,names(c)))
# index <- c(rep(NA,nlayers(c)))
# 
# for(i in c(1:12)){
#   d <- subset(c,month.indices[[i]])
#   e <- mean(d)
#   if (i ==1){
#     ndvi <- stack(e)
#   }
#   else {
#     ndvi <- addLayer(ndvi,e)
#   }
#   print(paste(i,"/12",sep=""))
# }
# writeRaster(ndvi,"ndvi_monthly.tif")

ndvi <- crop(stack("~/Dropbox/phenology/ndvi_monthly.tif"),ext)
ndvi.winter <- (ndvi[[1]]+ndvi[[2]]+ndvi[[12]])/3
ndvi.breeding <- (ndvi[[4]]+ndvi[[5]]+ndvi[[6]])/3
ndvi.winter <- resample(ndvi.winter,clim.winter)
ndvi.breeding <- resample(ndvi.breeding,clim.breeding)
setwd("~/Documents/nicheTracker/")
print("NDVI loaded. objects: ndvi.winter and ndvi.breeding.")




