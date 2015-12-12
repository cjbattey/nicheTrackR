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
source("scripts/ClimDataSetup.R") 
source("scripts/NDVIDataSetup.R")
#source("scripts/setupLandcover.R") #static, no seasonal changes
source("scripts/nSpecies_Trochil.R")

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

results.df <- data.frame(wnt.res.I=numeric,wnt.brd.I=numeric,
                         res.brd.I=numeric,sim.wnt.res=numeric,sim.wnt.brd=numeric,sim.res.brd=numeric,
                         eq.wnt.res=numeric,eq.wnt.brd=numeric,eq.res.brd=numeric)
figures.list <- list()

for(i in c(1:12)) { #initiate loop over species index
  print("stacking background data and extracting point values")
  bg.wnt <- stack(clim.winter[[i]],nspecies.wnt,ndvi.winter[[i]])  #stack background data
  bg.brd <- stack(clim.breeding[[i]],nspecies.brd,ndvi.breeding[[i]])  #stack background data
  names(bg.wnt) <- c("pre","frs","dtr","tmp","n.species","ndvi")
  names(bg.brd) <- c("pre","frs","dtr","tmp","n.species","ndvi")
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
  pc.brd <- pca$li[(nrow(sp.res)+nrow(sp.wnt)+1):(nrow(sp.res)+nrow(sp.wnt)+nrow(sp.brd)),]  #get subset of PC coordinates for breeding species points 
  pc.bg.wnt <- pca$li[(nrow(sp.res)+nrow(sp.wnt)+nrow(sp.brd)+1):(nrow(sp.res)+nrow(sp.wnt)+nrow(sp.brd)+nrow(env.wnt)),]  #get subset of PC coordinates for environment available in winter
  pc.bg.brd <- pca$li[(nrow(sp.res)+nrow(sp.wnt)+nrow(sp.brd)+nrow(env.wnt)+1):(nrow(sp.res)+nrow(sp.wnt)+nrow(sp.brd)+nrow(env.wnt)+nrow(env.brd)),]  #get subset of PC coordinates for environment available in breeding season
  pc.bg.all <- pca$li[(nrow(sp.res)+nrow(sp.wnt)+nrow(sp.brd)+1):(nrow(sp.res)+nrow(sp.wnt)+nrow(sp.brd)+nrow(env.wnt)+nrow(env.brd)),]  #get subset of PC coordinates for full annual environment
  
  print("Estimate kernel density in niche space")
  grid.res <- ecospat.grid.clim.dyn(pc.bg.all,pc.bg.brd,pc.res,R=100)
  grid.wnt <- ecospat.grid.clim.dyn(pc.bg.all,pc.bg.wnt,pc.wnt,R=100)
  grid.brd <- ecospat.grid.clim.dyn(pc.bg.all,pc.bg.brd,pc.brd,R=100)
  
  print("Conducting niche similarity and equivalency tests (see Broenniman 2012)")
  sim.wnt.res <- ecospat.niche.similarity.test(grid.wnt,grid.res,rep=100,one.sided=F)
  sim.wnt.brd <- ecospat.niche.similarity.test(grid.wnt,grid.brd,rep=100,one.sided=F)
  sim.res.brd <- ecospat.niche.similarity.test(grid.res,grid.brd,rep=100,one.sided=F)
  eq.wnt.res <- ecospat.niche.equivalency.test(grid.wnt,grid.res,rep=10)
  eq.wnt.brd <- ecospat.niche.equivalency.test(grid.wnt,grid.brd,rep=10)
  eq.res.brd <- ecospat.niche.equivalency.test(grid.res,grid.brd,rep=10)
  
  print("running maxent analyses")
  m.res <- maxent(bg.brd,occ.wnt[[i]])
  m.wnt <- maxent(bg.wnt,occ.wnt[[i]])
  m.brd <- maxent(bg.brd,occ.brd[[i]])
  p.res <- predict(m.res,bg.brd)
  p.wnt <- predict(m.wnt,bg.wnt)
  p.brd <- predict(m.brd,bg.brd)
  p.res.wnt <- predict(m.res,bg.wnt)
  p.wnt.brd <- predict(m.wnt,bg.brd)
  p.brd.wnt <- predict(m.brd,bg.wnt)
  
  sim.m.wnt.res <- nicheOverlap(p.wnt,p.res)
  sim.m.wnt.brd <- nicheOverlap(p.wnt,p.brd)
  
  results.df <- rbind(results.df,c(sim.wnt.res$obs$I,sim.wnt.brd$obs$I,sim.res.brd$obs$I,
                                   sim.wnt.res$p.I,sim.wnt.brd$p.I,sim.res.brd$p.I,eq.wnt.res$p.I,eq.wnt.brd$p.I,
                                   eq.res.brd$p.I))
  
  print("Saving figures")
  pdf(file=paste("./figures/",names(wnt.months[i]),"_nicheFigs.pdf",sep=""),width=4,height=4)
  ecospat.plot.niche(grid.res,title=paste(names(wnt.months[i]),"Resident Niche"))
  ecospat.plot.niche(grid.wnt,title=paste(names(wnt.months[i]),"Nonbreeding Niche"))
  ecospat.plot.niche(grid.brd,title=paste(names(wnt.months[i]),"Breeding Niche"))
  ecospat.plot.niche.dyn(grid.wnt,grid.res,quant=0.02,title=paste(names(wnt.months[i]),"Nonbreeding vs. Resident\n Niche Overlap"))
  ecospat.plot.niche.dyn(grid.wnt,grid.brd,quant=0.02,title=paste(names(wnt.months[i]),"Nonbreeding vs. Breeding\n Niche Overlap"))
  ecospat.plot.niche.dyn(grid.res,grid.brd,quant=0.02,title=paste(names(wnt.months[i]),"Resident vs. Breeding\n Niche Overlap"))
  ecospat.plot.overlap.test(sim.wnt.res,type="I",title=paste(names(wnt.months[i]),"Nonbreeding vs. Resident\n Niche Similarity"))
  ecospat.plot.overlap.test(sim.wnt.brd,type="I",title=paste(names(wnt.months[i]),"Nonbreeding vs. Breeding\n Niche Similarity"))
  ecospat.plot.overlap.test(sim.res.brd,type="I",title=paste(names(wnt.months[i]),"Resident vs. Breeding\n Niche Similarity"))
  ecospat.plot.overlap.test(eq.wnt.res,type="I",title=paste(names(wnt.months[i]),"Nonbreeding vs. Resident\n Niche Equivalency"))
  ecospat.plot.overlap.test(eq.wnt.brd,type="I",title=paste(names(wnt.months[i]),"Nonbreeding vs. Breeding\n Niche Equivalency"))
  ecospat.plot.overlap.test(eq.res.brd,type="I",title=paste(names(wnt.months[i]),"Resident vs. Breeding\n Niche Equivalency"))
  ecospat.plot.contrib(pca$c1,pca$eig)
  plot(p.res,main="resident")
  plot(p.wnt)
  plot(p.brd)
  plot(p.res.wnt)
  plot(p.wnt.brd)
  plot(p.brd.wnt)
  ggplot()+theme_bw()+
    geom_point(data=pc.brd,col="red",aes(x=Axis1,y=Axis2))+
    geom_point(data=pc.res,col="green",aes(x=Axis1,y=Axis2))+
    geom_point(data=pc.wnt,col="blue",aes(x=Axis1,y=Axis2))
  
  dev.off()
  
  print(paste(names(wnt.months[i]),"analysis complete"))
}

colnames(results.df) <- c("wnt.res.I","wnt.brd.I",
                           "res.brd.I","sim.wnt.res","sim.wnt.brd","sim.res.brd",
                           "eq.wnt.res","eq.wnt.brd","eq.res.brd")
rownames(results.df) <- names(wnt.months)
results.df$species <- names(wnt.months)

ggplot(data=results.df,aes(x=wnt.res.I,y=wnt.brd.I,col=species),size=4)+
  xlim(0,1)+ylim(0,1)+geom_point()+theme_bw()+geom_abline(intercept=0,slope=1,linetype="dotted")
  
  




