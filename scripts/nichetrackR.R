##niche trackR##
#Building and cross-predicting ENM's for wintering and breeding season hummingbirds with PCA and Maxent 
library(ecospat);library(ggplot2);library(raster);library(dismo);library(rgeos);library(ade4)
setwd("~/Documents/nicheTracker/")

#set extent and resolution (max resolution is 10min (e.g. 1/6 degree))
ext <- c(-155,-60,7,75)
r <- raster(xmn=-180,xmx=180,ymn=-90,ymx=90,res=1/6)

#set breeding and wintering months for each species
brd.months <- list(alhu=c(3,4,5),anhu=c(12,1,2),bchu=c(5,6),bthu=c(5,6,7),buhu=c(6,7),cahu=c(6,7),cohu=c(2:4),
                   luhu=c(6,7),rthu=c(6,7,8),ruhu=c(5,6),schu=c(6,7),vohu=c(6,7))
wnt.months <- list(alhu=c(10,11,12),anhu=c(8,9,10),bchu=c(1,2),bthu=c(11,12,1),buhu=c(12,1),cahu=c(12,1),cohu=c(2:4),
                   luhu=c(12,1),rthu=c(12,1,2),ruhu=c(12,1),schu=c(12,1),vohu=c(12,1))

#read in and format occurrence data (~10 min). Returns two lists of spatialpoints for each species in months given
#by brd.months and wnt.months. occurrence lists are occ.wnt and occ.brd.
source("scripts/OccDataSetup.R")

#read in and format climate and other data (~30min to read in and process all. 2hrs+ from raw data.)
source("scripts/ClimDataSetup.R") #need to adjust to output a list of breeding and a list of wintering climate stacks
source("scripts/NDVIDataSetup.R") #adjust to output a list of breeding and wintering ndvi's
source("scripts/setupLandcover.R") #static, no seasonal changes
source("scripts/nSpecies_Trochil.R") #ready to go - breeding and wintering applicable across species. 

#Analysis:
#1. loop over species
#2. build background data stack of relevant list elements
#3. extract bg and occ data and rbind together as one giant df
#4. build a row index file, including "resident" values (wintering locations, breeding season data)
  # order is: occ.res,occ.winter,occ.breeding,env.winter,env.breeding
#5. run pca
#6. grid and kernel density (ecospat.grid.clim.dyn) on breeding, wintering, resident.
#7. store figures of all of the above
#8. niche similarity test + store output
#9. niche equivalency test + store output

results.df <- data.frame(species=character,wnt.res.I=numeric,wnt.brd.I=numeric,
                         res.brd.I=numeric,sim.wnt.res=numeric,sim.wnt.brd=numeric,sim.res.brd=numeric,
                         eq.wnt.res=numeric,eq.wnt.brd=numeric,eq.res.brd=numeric)
figures.list <- list()

for(i in c(1:12)) { #initiate loop over species index
  print("stacking background data and extracting point values")
  bg.wnt <- stack(clim.winter[[i]],nspecies.wnt,ndvi.winter[[i]],lc.all)  #stack background data
  bg.brd <- stack(clim.breeding[[i]],nspecies.brd,ndvi.breeding[[i]],lc.all)  #stack background data
  names(bg.wnt) <- c("pre","frs","dtr","tmp","n.species","ndvi","pc.forest","pc.woodland","pc.shrub")
  names(bg.brd) <- c("pre","frs","dtr","tmp","n.species","ndvi","pc.forest","pc.woodland","pc.shrub")
  #writeRaster(bg.wnt,paste("envData_wnt_",wnt.months[[i]],".tif",sep="")) #write raster for background values
  sp.res <- na.omit(extract(bg.brd,occ.wnt[[i]]))  #extract hypothetical "resident" bg values
  sp.wnt <- na.omit(extract(bg.wnt,occ.wnt[[i]]))  #extract bg data for nb season
  sp.brd <- na.omit(extract(bg.wnt,occ.brd[[i]]))  #extract bg data for b season
  env.wnt <- na.omit(extract(bg.wnt,as.data.frame(bg.wnt,xy=T)[1:2])) #extract bg data as df
  env.brd <- na.omit(extract(bg.brd,as.data.frame(bg.brd,xy=T)[1:2])) #extract bg data as df
  df <- rbind(sp.res,sp.wnt,sp.brd,env.wnt,env.brd) #rbind full dataset for PCA

  weights <- c(rep(0,nrow(sp.res)),rep(0,nrow(sp.wnt)),rep(0,nrow(sp.brd)),rep(1,nrow(env.wnt)),rep(1,nrow(env.brd)))
  
  print("Conduct PCA")
  pca <- dudi.pca(df,row.w=weights,nf=2,scannf=F,center=T,scale=T) #run PCA
  
  pc.res <- pca$li[1:nrow(sp.res),] #get subset of PC coordinates for "resident" species points
  pc.wnt <- pca$li[nrow(sp.res):(nrow(sp.res)+nrow(sp.wnt)),]  #get subset of PC coordinates for wintering species points 
  pc.brd <- pca$li[(nrow(sp.res)+nrow(sp.wnt)):(nrow(sp.res)+nrow(sp.wnt)+nrow(sp.brd)),]  #get subset of PC coordinates for breeding species points 
  pc.bg.wnt <- pca$li[(nrow(sp.res)+nrow(sp.wnt)+nrow(sp.brd)):(nrow(sp.res)+nrow(sp.wnt)+nrow(sp.brd)+nrow(env.wnt)),]  #get subset of PC coordinates for environment available in winter
  pc.bg.brd <- pca$li[(nrow(sp.res)+nrow(sp.wnt)+nrow(sp.brd)+nrow(env.wnt)):(nrow(sp.res)+nrow(sp.wnt)+nrow(sp.brd)+nrow(env.wnt)+nrow(env.brd)),]  #get subset of PC coordinates for environment available in breeding season
  pc.bg.all <- pca$li[(nrow(sp.res)+nrow(sp.wnt)+nrow(sp.brd)):(nrow(sp.res)+nrow(sp.wnt)+nrow(sp.brd)+nrow(env.wnt)+nrow(env.brd)),]  #get subset of PC coordinates for full annual environment
  
  print("Estimate kernel density in niche space")
  grid.res <- ecospat.grid.clim.dyn(pc.bg.all,pc.bg.brd,pc.res,R=100)
  grid.wnt <- ecospat.grid.clim.dyn(pc.bg.all,pc.bg.wnt,pc.wnt,R=100)
  grid.brd <- ecospat.grid.clim.dyn(pc.bg.all,pc.bg.brd,pc.brd,R=100)
  
  print("Conducting niche similarity and equivalency tests (see Broenniman 2012)")
  sim.wnt.res <- ecospat.niche.similarity.test(grid.wnt,grid.res,rep=100,one.sided=F)
  sim.wnt.brd <- ecospat.niche.similarity.test(grid.wnt,grid.brd,rep=100,one.sided=F)
  sim.res.brd <- ecospat.niche.similarity.test(grid.res,grid.brd,rep=100,one.sided=F)
  eq.wnt.res <- ecospat.niche.equivalency.test(grid.wnt,grid.res,rep=100)
  eq.wnt.brd <- ecospat.niche.equivalency.test(grid.wnt,grid.brd,rep=100)
  eq.res.brd <- ecospat.niche.equivalency.test(grid.res,grid.brd,rep=100)
  
  results.df <- rbind(results.df,c(names(wnt.months[i]),sim.wnt.res$obs$I,sim.wnt.brd$obs$I,sim.res.brd$obs$I,
                                   sim.wnt.res$p.I,sim.wnt.brd$p.I,sim.res.brd$p.I,eq.wnt.res$p.I,eq.wnt.brd$p.I,
                                   eq.res.brd$p.I))
  
  print("Saving figures")
  pdf(file=paste("./figures/",names(wnt.months[i]),"_nicheFigs.pdf",sep=""),width=4,height=4)
  ecospat.plot.niche(grid.res,title=paste(names(wnt.months[i]),"Resident Niche"))
  ecospat.plot.niche(grid.wnt,title=paste(names(wnt.months[i]),"Nonbreeding Niche"))
  ecospat.plot.niche(grid.brd,title=paste(names(wnt.months[i]),"Breeding Niche"))
  ecospat.plot.niche.dyn(grid.wnt,grid.res,title=paste(names(wnt.months[i]),"Nonbreeding vs. Resident\n Niche Overlap"))
  ecospat.plot.niche.dyn(grid.wnt,grid.brd,title=paste(names(wnt.months[i]),"Nonbreeding vs. Breeding\n Niche Overlap"))
  ecospat.plot.niche.dyn(grid.res,grid.brd,title=paste(names(wnt.months[i]),"Resident vs. Breeding\n Niche Overlap"))
  ecospat.plot.overlap.test(sim.wnt.res,type="I",title=paste(names(wnt.months[i]),"Nonbreeding vs. Resident\n Niche Similarity"))
  ecospat.plot.overlap.test(sim.wnt.brd,type="I",title=paste(names(wnt.months[i]),"Nonbreeding vs. Breeding\n Niche Similarity"))
  ecospat.plot.overlap.test(sim.res.brd,type="I",title=paste(names(wnt.months[i]),"Resident vs. Breeding\n Niche Similarity"))
  ecospat.plot.overlap.test(eq.wnt.res,type="I",title=paste(names(wnt.months[i]),"Nonbreeding vs. Resident\n Niche Equivalency"))
  ecospat.plot.overlap.test(eq.wnt.brd,type="I",title=paste(names(wnt.months[i]),"Nonbreeding vs. Breeding\n Niche Equivalency"))
  ecospat.plot.overlap.test(eq.res.brd,type="I",title=paste(names(wnt.months[i]),"Resident vs. Breeding\n Niche Equivalency"))
  ecospat.plot.contrib(pca$c1,pca$eig)
  dev.off()
  
  print(paste(names(wnt.months[i])),"analysis complete")
}



# 
# Old analyses for single species (ruhu) and data exploration
# 
# 
# 
# 
# 
# 
# 
# #match resolutions, add ndvi and landcover to climate data
# env.w.raster <- stack(clim.winter,lc.all,ndvi.winter,n.species)
# env.b.raster <- stack(clim.breeding,lc.all,ndvi.winter,n.species)
# names(all.winter) <- c("pre","frs","dtr","forest","woodland","shrub","ndvi","n.species")
# names(all.breeding) <- c("pre","frs","dtr","forest","woodland","shrub","ndvi","n.species")
# 
# #overlay occurrences on climate rasters, extract values. 
# ruhu.brd.crop <- SpatialPoints(gridSample(ruhu.brd,clim.winter[[1]])) #subsample breeding data to 1 report per grid cell
# occ.breeding <- na.omit(extract(clim.breeding,ruhu.brd.crop))[round(runif(100,min=1,max=1671)),]
# occ.winter <- na.omit(extract(clim.winter,ruhu.wint))[round(runif(100,min=1,max=165)),]
# occ.res <- na.omit(extract(clim.breeding,ruhu.wint))[round(runif(100,min=1,max=165)),]
# occ.all <- rbind(occ.winter,occ.breeding)
# 
# #background climate data for PCA
# env.winter <- na.omit(extract(clim.winter,as.data.frame(clim.winter,xy=T)[1:2])) #82315 rows per global 10' dataset
# env.breeding <- na.omit(extract(clim.breeding,as.data.frame(clim.breeding,xy=T)[1:2]))
# env.all <- rbind(env.winter,env.breeding) #164630 total
# 
# #combined dataset
# df <- rbind(occ.winter,occ.breeding,env.winter,env.breeding)
# df.res <- rbind(occ.res,occ.winter,occ.breeding,env.winter,env.breeding)
# 
# #row weighting vectors to train PCA on subsets of the data.
# weights.env.all <- c(rep(0,nrow(occ.all)),rep(1,nrow(env.all)))
# weights.env.brd <- c(rep(0,nrow(occ.all)),rep(1,nrow(env.breeding)),rep(0,nrow(env.winter)))
# weights.env.wnt <- c(rep(0,nrow(occ.all)),rep(0,nrow(env.breeding)),rep(1,nrow(env.winter)))
# weights.occ.all <- c(rep(1,nrow(occ.all)),rep(0,nrow(env.all)))
# weights.occ.brd <- c(rep(0,nrow(occ.breeding)),rep(0,nrow(occ.winter)),rep(1,nrow(env.all)))
# weights.occ.wnt <- c(rep(0,nrow(occ.breeding)),rep(1,nrow(occ.winter)),rep(0,nrow(env.all)))
# weights.occ.res <- c(rep(1,nrow(occ.res)),rep(0,nrow(occ.breeding)),rep(0,nrow(occ.winter)),rep(0,nrow(env.all)))
# 
# #run PCA's
# pca.env <- dudi.pca(df,row.w=weights.env.all,nf=2,scannf=F) #train on all available environmental variables 
# pca.env.brd <- dudi.pca(df,row.w=weights.env.brd,scannf=F) #train on breeding season climate 
# pca.env.wnt <- dudi.pca(df,row.w=weights.env.wnt,nf=2,scannf=F) #train on wintering climate
# pca.occ <- dudi.pca(df,row.w=weights.occ.all,nf=2,scannf=F) #train on climate of all occurrences
# pca.occ.brd <- dudi.pca(df,row.w=weights.occ.brd,nf=2,scannf=F) #train on breeding occurrence climate
# pca.occ.wnt <- dudi.pca(df,row.w=weights.env.wnt,nf=2,scannf=F) #train on wintering ccurrence climate
# pca.env.res <- dudi.pca(df.res,row.w=weights.occ.res,nf=2,scannf=F) #include hypothetical "resident" values
# 
# allClim <- pca.env.res$li[]
# 
# #smooth observations within environmental space
# grid.winter <- ecospat.grid.clim.dyn(pca.env$li[201:nrow(env.all),],pca.env$li[201:nrow(env.winter),],sp=pca.env$li[1:100,],R=100)
# grid.breeding <- ecospat.grid.clim.dyn(pca.env$li[201:nrow(env.all),],pca.env$li[(nrow(env.winter)+201):(nrow(env.all)+200),],sp=pca.env$li[101:200,],R=100)
# grid.resident <- ecospat.grid.clim.dyn(pca.env.res$li[301:nrow(env.all),],pca.env.res$li[(nrow(env.winter)+301):(nrow(env.all)+300),],sp=pca.env.res$li[1:100,],R=100)
# 
# ecospat.plot.niche(z=grid.resident)
# ecospat.plot.niche(z=grid.winter)
# ecospat.plot.niche(z=grid.breeding)
# ecospat.plot.niche.dyn(grid.winter,grid.breeding)
# ecospat.plot.niche.dyn(grid.resident,grid.breeding,interest=2)
# ecospat.niche.overlap(grid.resident,grid.winter,cor=T)
# ecospat.niche.overlap(grid.breeding,grid.winter,cor=T)
# 
# simTest <- ecospat.niche.similarity.test(grid.winter,grid.breeding,rep=100,one.sided=F)
# simTest2 <- ecospat.niche.similarity.test(grid.breeding,grid.winter,rep=100,one.sided=F)
# ecospat.plot.overlap.test(simTest,type="D",title="RUHU breeding vs. wintering range")
# ecospat.plot.overlap.test(simTest2,type="D",title="RUHU breeding vs. wintering range")
# 
# eqTest.wb <- ecospat.niche.equivalency.test(grid.winter,grid.breeding,100)
# ecospat.plot.contrib(pca.env.res$c1,pca.env.res$eig)
# occ.scores <- pca.env$li[1:300,]
# occ.scores$season <- c(rep("resident",nrow(occ.res)),rep("winter",nrow(occ.winter)),rep("breeding",nrow(occ.breeding)))
# 
# ggplot(occ.scores,aes(x=Axis1,y=Axis2,col=season))+theme_bw()+geom_point()
# 
# 
# #### maxent analyses
# m.winter <- maxent(clim.winter,ruhu.wint)
# m.breeding <- maxent(clim.breeding,ruhu.brd.crop)
# m.resident <- maxent(clim.breeding,ruhu.wint)
# 
# p.resident <- predict(m.resident,clim.breeding)
# p.winter <- predict(m.winter,clim.winter)
# p.wb <- predict(m.winter,clim.breeding)
# p.breeding <- predict(m.breeding,clim.breeding)
# p.bw <- predict(m.breeding,clim.winter)
# 

