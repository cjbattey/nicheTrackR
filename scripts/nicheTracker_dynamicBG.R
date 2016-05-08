#niche tracking analysis function using PCA & kernel density estimator via ecospat
#improvements: should limit background to a 5000km(?) buffer around the full range for each species. 

nicheTracker <- function(i) {
  require(data.table);require(dismo);require(ecospat);require(pbapply);require(rgeos);require(ggplot2)  
  
  #subset species/seasons; remove occurrences outside known range
  loc <- subset(gbif,species == i) 
  loc.sum <- subset(loc,month%in%c(6,7,8)==T)
  loc.wnt <- subset(loc,month%in%c(12,1,2)==T)
  loc.sum <- SpatialPoints(data.frame(loc.sum[,.(decimallongitude,decimallatitude)]), proj4string=CRS(proj4string(ranges[[i]])))
  loc.wnt <- SpatialPoints(data.frame(loc.wnt[,.(decimallongitude,decimallatitude)]), proj4string=CRS(proj4string(ranges[[i]])))
  #loc.sum <- gIntersection(loc.sum,ranges[[i]][ranges[[i]]@data$SEASONAL == 2 | ranges[[i]]@data$SEASONAL == 1,])
  #loc.wnt <- gIntersection(loc.wnt,ranges[[i]][ranges[[i]]@data$SEASONAL == 3 | ranges[[i]]@data$SEASONAL == 1,]) #problem: austral migrants!
  
  #note: vireos fails on task 13 bc there are no wintering localities w/in (tiny) mapped range...
  #options          - add an if statement and run intersect only if total area is > n
  #**this for now   - add if statement after filtering and don't analyze species w/o n reports within suitable range areas 
  #                 ---note this will cause load balancing issues with doMC so won't run at 100% efficiency :(
  #do this later?   - Could extract from the full range for sp. w/o occurrence data. 
  #                 ---can't do "resident" niche that way tho
  if(length(loc.wnt) > 15 & length(loc.sum) > 15){
    
  #find centroid of seasonal occurrences
  centroid.sum <- gCentroid(loc.sum)
  centroid.wnt <- gCentroid(loc.wnt)
  centroid.distance <- pointDistance(centroid.sum,centroid.wnt,lonlat=T)
  
  #extract occurrence pt background data
  sp.sum <- na.omit(extract(bg.sum.r,loc.sum))
  sp.wnt <- na.omit(extract(bg.wnt.r,loc.wnt))
  sp.res.s <- na.omit(extract(bg.wnt.r,loc.sum)) #resident on summer territory
  sp.res.w <- na.omit(extract(bg.sum.r,loc.wnt)) #resident on winter territory
  
  #extract available background data
  bg.sum.df <- na.omit(data.frame(extract(bg.sum.r,ranges.buffered[[i]])))
  bg.wnt.df <- na.omit(data.frame(extract(bg.wnt.r,ranges.buffered[[i]])))
  
  #combo dataset for pca
  df <- rbind(sp.sum,sp.wnt,sp.res.s,sp.res.w,bg.sum.df,bg.wnt.df)
  
  #set weighting to train PCA on full background data
  weights <- c(rep(0,nrow(sp.sum)),rep(0,nrow(sp.wnt)),rep(0,nrow(sp.res.s)),
               rep(0,nrow(sp.res.w)),rep(1,nrow(bg.sum.df)),rep(1,nrow(bg.wnt.df))) 
  
  #run pca
  pca <- dudi.pca(df,row.w=weights,nf=2,scannf=F,center=T,scale=T) 
  
  #pull rows for species & background pc coords
  pc.sum <- pca$li[1 : nrow(sp.sum),] 
  pc.wnt <- pca$li[(1+nrow(sp.sum)):(nrow(sp.sum)+nrow(sp.wnt)),]
  pc.res.s <- pca$li[(1+nrow(sp.sum)+nrow(sp.wnt)) : (1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)),]
  pc.res.w <- pca$li[(1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)) : (1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)+nrow(sp.res.w)),]
  pc.bg.sum <- pca$li[(1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)+nrow(sp.res.w)) : (1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)+nrow(sp.res.w)+nrow(bg.sum.df)),]
  pc.bg.wnt <- pca$li[(1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)+nrow(sp.res.w)+nrow(bg.sum.df)) : nrow(pca$li),]
  pc.bg.all <- pca$li[(1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)+nrow(sp.res.w)) : nrow(pca$li),]
  rm(list=c("pca"))
  
  #grid with kernal density estimator (see Broenniman et al. 2012, J. Biogeography)
  grid.sum <- ecospat.grid.clim.dyn(glob=pc.bg.all,glob1=pc.bg.sum,sp=pc.sum,R=100) 
  grid.wnt <- ecospat.grid.clim.dyn(glob=pc.bg.all,glob1=pc.bg.wnt,sp=pc.wnt,R=100)
  grid.res.s <- ecospat.grid.clim.dyn(glob=pc.bg.all,glob1=pc.bg.wnt,sp=pc.res.s,R=100)
  grid.res.w <- ecospat.grid.clim.dyn(glob=pc.bg.all,glob1=pc.bg.sum,sp=pc.res.w,R=100)
  
  #run similarity test
  sim.s.w <- ecospat.niche.similarity.test.noplot(z1=grid.sum,z2=grid.wnt,rep=1000,one.sided=F)
  overlap.res.s <- ecospat.niche.overlap(grid.sum,grid.res.s,cor=T)
  overlap.res.w <- ecospat.niche.overlap(grid.wnt,grid.res.w,cor=T)
  
  #experimental maxent runs
  #bg.sum.clim <- dropLayer(bg.sum.r,c(5,6,7))
  #bg.wnt.clim <- dropLayer(bg.wnt.r,c(5,6,7))
  
  #max.sum <- maxent(bg.sum.clim,loc.sum)
  #pred.sum <- predict(max.sum,bg.sum.clim)
  #plot(pred.sum)
  
  data.frame(species=paste(i),centroid.distance=centroid.distance,I.obs=sim.s.w$obs$I,I.res.s=overlap.res.s$I,I.res.w=overlap.res.w$I,p.similar=sim.s.w$p.I)
  }
}

