##niche trackR##
#Building and cross-predicting ENM's for wintering and breeding season hummingbirds with PCA and Maxent 
library(ecospat);library(ggplot2);library(raster);library(dismo);library(rgeos);library(ade4)
setwd("~/Documents/nicheTracker/")

#read in and reformat occurrence and climate datasets. See screen log for object names and formats.
source("scripts/ClimDataSetup.R")
source("scripts/OccDataSetup.R")
source("scripts/NDVIDataSetup.R")
source("scripts/setupLandcover.R")

#match resolutions, add ndvi and landcover to climate data
ndvi.winter <- resample(ndvi.winter,clim.winter)
ndvi.breeding <- resample(ndvi.breeding,clim.breeding)
clim.winter <- addLayer(clim.winter,ndvi.winter)
clim.breeding <- addLayer(clim.breeding,ndvi.breeding)
names(clim.winter[[4]]) <- "NDVI"
names(clim.breeding[[4]]) <- "NDVI"


#overlay occurrences on climate rasters, extract values. 
ruhu.brd.crop <- SpatialPoints(gridSample(ruhu.brd,clim.winter[[1]])) #subsample breeding data to 1 report per grid cell
occ.breeding <- na.omit(extract(clim.breeding,ruhu.brd.crop))[round(runif(100,min=1,max=1671)),]
occ.winter <- na.omit(extract(clim.winter,ruhu.wint))[round(runif(100,min=1,max=165)),]
occ.res <- na.omit(extract(clim.breeding,ruhu.wint))[round(runif(100,min=1,max=165)),]
occ.all <- rbind(occ.winter,occ.breeding)

#background climate data for PCA
env.winter <- na.omit(extract(clim.winter,as.data.frame(clim.winter,xy=T)[1:2])) #82315 rows per global 10' dataset
env.breeding <- na.omit(extract(clim.breeding,as.data.frame(clim.breeding,xy=T)[1:2]))
env.all <- rbind(env.winter,env.breeding) #164630 total

#combined dataset
df <- rbind(occ.winter,occ.breeding,env.winter,env.breeding)
df.res <- rbind(occ.res,occ.winter,occ.breeding,env.winter,env.breeding)

#row weighting vectors to train PCA on subsets of the data.
weights.env.all <- c(rep(0,nrow(occ.all)),rep(1,nrow(env.all)))
weights.env.brd <- c(rep(0,nrow(occ.all)),rep(1,nrow(env.breeding)),rep(0,nrow(env.winter)))
weights.env.wnt <- c(rep(0,nrow(occ.all)),rep(0,nrow(env.breeding)),rep(1,nrow(env.winter)))
weights.occ.all <- c(rep(1,nrow(occ.all)),rep(0,nrow(env.all)))
weights.occ.brd <- c(rep(0,nrow(occ.breeding)),rep(0,nrow(occ.winter)),rep(1,nrow(env.all)))
weights.occ.wnt <- c(rep(0,nrow(occ.breeding)),rep(1,nrow(occ.winter)),rep(0,nrow(env.all)))
weights.occ.res <- c(rep(1,nrow(occ.res)),rep(0,nrow(occ.breeding)),rep(0,nrow(occ.winter)),rep(0,nrow(env.all)))

#run PCA's
pca.env <- dudi.pca(df,row.w=weights.env.all,nf=2,scannf=F) #train on all available environmental variables 
pca.env.brd <- dudi.pca(df,row.w=weights.env.brd,scannf=F) #train on breeding season climate 
pca.env.wnt <- dudi.pca(df,row.w=weights.env.wnt,nf=2,scannf=F) #train on wintering climate
pca.occ <- dudi.pca(df,row.w=weights.occ.all,nf=2,scannf=F) #train on climate of all occurrences
pca.occ.brd <- dudi.pca(df,row.w=weights.occ.brd,nf=2,scannf=F) #train on breeding occurrence climate
pca.occ.wnt <- dudi.pca(df,row.w=weights.env.wnt,nf=2,scannf=F) #train on wintering ccurrence climate
pca.env.res <- dudi.pca(df.res,row.w=weights.occ.res,nf=2,scannf=F) #include hypothetical "resident" values

allClim <- pca.env.res$li[]

#smooth observations within environmental space
grid.winter <- ecospat.grid.clim.dyn(pca.env$li[201:nrow(env.all),],pca.env$li[201:nrow(env.winter),],sp=pca.env$li[1:100,],R=100)
grid.breeding <- ecospat.grid.clim.dyn(pca.env$li[201:nrow(env.all),],pca.env$li[(nrow(env.winter)+201):(nrow(env.all)+200),],sp=pca.env$li[101:200,],R=100)
grid.resident <- ecospat.grid.clim.dyn(pca.env.res$li[301:nrow(env.all),],pca.env.res$li[(nrow(env.winter)+301):(nrow(env.all)+300),],sp=pca.env.res$li[1:100,],R=100)

ecospat.plot.niche(z=grid.resident)
ecospat.plot.niche(z=grid.winter)
ecospat.plot.niche(z=grid.breeding)
ecospat.plot.niche.dyn(grid.winter,grid.breeding)
ecospat.plot.niche.dyn(grid.resident,grid.breeding,interest=2)
ecospat.niche.overlap(grid.resident,grid.winter,cor=T)
ecospat.niche.overlap(grid.breeding,grid.winter,cor=T)

simTest <- ecospat.niche.similarity.test(grid.winter,grid.breeding,rep=100,one.sided=F)
simTest2 <- ecospat.niche.similarity.test(grid.breeding,grid.winter,rep=100,one.sided=F)
ecospat.plot.overlap.test(simTest,type="D",title="RUHU breeding vs. wintering range")
ecospat.plot.overlap.test(simTest2,type="D",title="RUHU breeding vs. wintering range")

eqTest.wb <- ecospat.niche.equivalency.test(grid.winter,grid.breeding,100)
ecospat.plot.contrib(pca.env.res$c1,pca.env.res$eig)
occ.scores <- pca.env$li[1:300,]
occ.scores$season <- c(rep("resident",nrow(occ.res)),rep("winter",nrow(occ.winter)),rep("breeding",nrow(occ.breeding)))

ggplot(occ.scores,aes(x=Axis1,y=Axis2,col=season))+theme_bw()+geom_point()


#### maxent analyses
m.winter <- maxent(clim.winter,ruhu.wint)
m.breeding <- maxent(clim.breeding,ruhu.brd.crop)
m.resident <- maxent(clim.breeding,ruhu.wint)

p.resident <- predict(m.resident,clim.breeding)
p.winter <- predict(m.winter,clim.winter)
p.wb <- predict(m.winter,clim.breeding)
p.breeding <- predict(m.breeding,clim.breeding)
p.bw <- predict(m.breeding,clim.winter)


