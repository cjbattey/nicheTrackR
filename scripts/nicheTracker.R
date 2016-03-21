nicheTracker <- function(i) {
  require(data.table);require(dismo);require(ecospat);require(pbapply);require(rgeos);require(ggplot2)  
  
  loc <- subset(gbif,species == i) 
  loc.sum <- subset(loc,month%in%c(5,6,7)==T)
  loc.wnt <- subset(loc,month%in%c(11,12,1)==T)
  loc.sum <- SpatialPoints(data.frame(loc.sum[,.(decimallongitude,decimallatitude)]), proj4string=CRS(proj4string(alt)))
  loc.wnt <- SpatialPoints(data.frame(loc.wnt[,.(decimallongitude,decimallatitude)]), proj4string=CRS(proj4string(alt)))
  
  #find centroid of seasonal occurrences
  centroid.sum <- gCentroid(loc.sum)
  centroid.wnt <- gCentroid(loc.wnt)
  centroid.distance <- pointDistance(centroid.sum,centroid.wnt,lonlat=T)
  
  #extract background data
  sp.sum <- na.omit(extract(bg.sum.r,loc.sum))
  sp.wnt <- na.omit(extract(bg.wnt.r,loc.wnt))
  sp.res.s <- na.omit(extract(bg.wnt.r,loc.sum)) #resident on summer territory
  sp.res.w <- na.omit(extract(bg.sum.r,loc.wnt)) #resident on winter territory
  
  #continue if > N reports
  if(nrow(sp.sum) > 15 & nrow(sp.wnt) > 15) {
    
  #combo dataset for pca
  df <- rbind(sp.sum,sp.wnt,sp.res.s,sp.res.w,data.frame(bg.sum.dt),data.frame(bg.wnt.dt)) 
  
  #set weighting to train PCA on full background data
  weights <- c(rep(0,nrow(sp.sum)),rep(0,nrow(sp.wnt)),rep(0,nrow(sp.res.s)),
               rep(0,nrow(sp.res.w)),rep(1,nrow(bg.sum.dt)),rep(1,nrow(bg.wnt.dt))) 
  
  
  #run pca
  pca <- dudi.pca(df,row.w=weights,nf=2,scannf=F,center=T,scale=T) 
  
  #pull rows for species & background pc coords
  pc.sum <- pca$li[1 : nrow(sp.sum),] 
  pc.wnt <- pca$li[(1+nrow(sp.sum)):(nrow(sp.sum)+nrow(sp.wnt)),]
  pc.res.s <- pca$li[(1+nrow(sp.sum)+nrow(sp.wnt)) : (1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)),]
  pc.res.w <- pca$li[(1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)) : (1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)+nrow(sp.res.w)),]
  pc.bg.sum <- pca$li[(1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)+nrow(sp.res.w)) : (1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)+nrow(sp.res.w)+nrow(bg.sum.dt)),]
  pc.bg.wnt <- pca$li[(1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)+nrow(sp.res.w)+nrow(bg.sum.dt)) : nrow(pca$li),]
  pc.bg.all <- pca$li[(1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)+nrow(sp.res.w)) : nrow(pca$li),]
  
  #grid with kernal density estimator (see Broenniman et al. 2012, J. Biogeography)
  grid.sum <- ecospat.grid.clim.dyn(glob=pc.bg.all,glob1=pc.bg.sum,sp=pc.sum,R=100) 
  grid.wnt <- ecospat.grid.clim.dyn(glob=pc.bg.all,glob1=pc.bg.wnt,sp=pc.wnt,R=100)
  grid.res.s <- ecospat.grid.clim.dyn(glob=pc.bg.all,glob1=pc.bg.wnt,sp=pc.res.s,R=100)
  grid.res.w <- ecospat.grid.clim.dyn(glob=pc.bg.all,glob1=pc.bg.sum,sp=pc.res.w,R=100)
  
  #run similarity test
  sim.s.w <- ecospat.niche.similarity.test.noplot(z1=grid.sum,z2=grid.wnt,rep=100,one.sided=F)
  overlap.res.s <- ecospat.niche.overlap(grid.sum,grid.res.s,cor=T)
  overlap.res.w <- ecospat.niche.overlap(grid.wnt,grid.res.w,cor=T)
  
  data.frame(species=paste(i),centroid.distance=centroid.distance,I.obs=sim.s.w$obs$I,I.res.s=overlap.res.s$I,I.res.w=overlap.res.w$I,p.similar=sim.s.w$p.I)
  }
}
