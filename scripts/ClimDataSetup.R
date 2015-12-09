##build climate datasets for niche tracking analyses: 
##output format should be x,y,Var1,Var2...VarN, where Var is averaged across either winter (Nov-Jan) or breeding (Apr-Jun) seasons
##Input is CRU2 - separate files for each Var in format: y x month1 month2 month3...month12
setwd("~/Documents/nicheTracker/")
library(raster);library(ecodist);library(plyr)

#Warning: pre and wnd are at slightly different resolutions to other variables - diff number of rows in files. Locations are NOT given 
#in same order among different variables (which: wtf. why?). 
print("reading in data")
pre <- read.table("./data/climate_CRU2/text/grid_10min_pre.dat",header=F,sep="",colClasses=c(rep("numeric",26)))
pre <- pre[,1:14]
frs <- read.table("./data/climate_CRU2/text/grid_10min_frs.dat",header=F,sep="",colClasses=c(rep("numeric",14)))
dtr <- read.table("./data/climate_CRU2/text/grid_10min_dtr.dat",header=F,sep="",colClasses=c(rep("numeric",14)))
#rdo <- read.table("./data/climate_CRU2/text/grid_10min_rd0.dat",header=F,sep="",colClasses=c(rep("numeric",14)))
#reh <-read.table("./data/climate_CRU2/text/grid_10min_reh.dat",header=F,sep="",colClasses=c(rep("numeric",14)))
#sunp <- read.table("./data/climate_CRU2/text/grid_10min_sunp.dat",header=F,sep="",colClasses=c(rep("numeric",14)))
tmp <- read.table("./data/climate_CRU2/text/grid_10min_tmp.dat",header=F,sep="",colClasses=c(rep("numeric",14)))
#wnd <- read.table("./data/climate_CRU2/text/grid_10min_wnd.dat",header=F,sep="",colClasses=c(rep("numeric",14)))
print("aligning points and reformatting")
clim.data <- list(pre,frs,dtr,tmp)
names(clim.data) <- c("pre","frs","dtr","tmp")

clim.data <- lapply(clim.data,FUN=function(e) setNames(e,c("lat","long","Jan","Feb","Mar",
                                                           "Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")))

r <- raster(xmn=-180,xmx=180,ymn=-90,ymx=90,res=1/6)
######prepare climate datasets for PCA's
#rasterize all datasets to harmonize for diff. number/order of points bw precip and all other vars
clim.data.wint <- lapply(clim.data,FUN=function(e) {
  rasterize(
    SpatialPointsDataFrame(data.frame(long=e$long,lat=e$lat),data=(data.frame(var=(e$Nov+e$Dec+e$Jan)/3)))
    ,r,field="var")
}
)
clim.winter <- crop(stack(clim.data.wint),ext)
#transform rasters back to data frame (now w/same points & order)
#clim.winter <- lapply(clim.data.wint,FUN=function(e) {
#  as.data.frame(e)
#})
#cbind into one df
#clim.winter <- do.call(cbind,clim.winter)
#names(clim.winter) <- names(clim.data)
#remove NA points (oceans) and add long,lat as first two columns
#clim.winter <- na.omit(data.frame(as.data.frame(r,xy=T)[,1:2],clim.winter))

##repeat the above for breeding season
clim.data.brd <- lapply(clim.data,FUN=function(e) {
  rasterize(
    SpatialPointsDataFrame(data.frame(long=e$long,lat=e$lat),data=(data.frame(var=(e$Apr+e$May+e$Jun)/3)))
    ,r,field="var")
}
)
clim.breeding <- crop(stack(clim.data.brd),ext)
#clim.breeding <- lapply(clim.data.brd,FUN=function(e) {
#  as.data.frame(e)
#})
#clim.breeding <- do.call(cbind,clim.breeding)
#names(clim.breeding) <- names(clim.data)
#clim.breeding <- na.omit(data.frame(as.data.frame(r,xy=T)[,1:2],clim.breeding))

print("Data laded. Format: raster stack at 10min resolution w/7 layers (see names for var). Objects: clim.winter & clim.breeding")







