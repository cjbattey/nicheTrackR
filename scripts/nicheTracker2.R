####Hummingbird seasonal movement and niche-tracking analysis
library(data.table);library(dismo);library(ecospat);library(pbapply);library(rgeos);library(ggplot2)
setwd("~/R/nicheTracker/scripts/")

#load data
source("setup_data.R")

#analysis loop
out <- data.frame(species=character(),centroid.distance=numeric(),I.obs=numeric(),I.res.s=numeric(),
                  I.res.w=numeric(),p.similar=numeric(),stringsAsFactors = F)

for(i in levels(gbif$species)){
#for(i in c("Calothorax lucifer","Calypte anna","Calypte costae","Selasphorus sasin")){  
  t1 <- Sys.time()
  #load species occurrence points, split by season
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
  
  #combo dataset for pca
  df <- rbind(sp.sum,sp.wnt,sp.res.s,sp.res.w,bg.sum,bg.wnt) 
  
  #set weighting to run PCA on full background data
  weights <- c(rep(0,nrow(sp.sum)),rep(0,nrow(sp.wnt)),rep(0,nrow(sp.res.s)),
               rep(0,nrow(sp.res.w)),rep(1,nrow(bg.sum)),rep(1,nrow(bg.wnt))) 
  
  #run pca
  pca <- dudi.pca(df,row.w=weights,nf=2,scannf=F,center=T,scale=T) 
  
  #pull rows for species & background pc coords
  pc.sum <- pca$li[1 : nrow(sp.sum),] 
  pc.wnt <- pca$li[(1+nrow(sp.sum)):(nrow(sp.sum)+nrow(sp.wnt)),]
  pc.res.s <- pca$li[(1+nrow(sp.sum)+nrow(sp.wnt)) : (1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)),]
  pc.res.w <- pca$li[(1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)) : (1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)+nrow(sp.res.w)),]
  pc.bg.sum <- pca$li[(1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)+nrow(sp.res.w)) : (1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)+nrow(sp.res.w)+nrow(bg.sum)),]
  pc.bg.wnt <- pca$li[(1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)+nrow(sp.res.w)+nrow(bg.sum)) : nrow(pca$li),]
  pc.bg.all <- pca$li[(1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)+nrow(sp.res.w)) : nrow(pca$li),]
  
  #grid with kernal density estimator (see Broenniman et al. 2012, J. Biogeography)
  grid.sum <- ecospat.grid.clim.dyn(glob=pc.bg.all,glob1=pc.bg.sum,sp=pc.sum,R=100) 
  grid.wnt <- ecospat.grid.clim.dyn(glob=pc.bg.all,glob1=pc.bg.wnt,sp=pc.wnt,R=100)
  grid.res.s <- ecospat.grid.clim.dyn(glob=pc.bg.all,glob1=pc.bg.wnt,sp=pc.res.s,R=100)
  grid.res.w <- ecospat.grid.clim.dyn(glob=pc.bg.all,glob1=pc.bg.sum,sp=pc.res.w,R=100)
  
  #run similarity test
  sim.s.w <- ecospat.niche.similarity.test(z1=grid.sum,z2=grid.wnt,rep=100,one.sided=F)
  overlap.obs <- ecospat.niche.overlap(grid.sum,grid.wnt,cor=T)
  overlap.res.s <- ecospat.niche.overlap(grid.sum,grid.res.s,cor=T)
  overlap.res.w <- ecospat.niche.overlap(grid.wnt,grid.res.w,cor=T)
  #sim.s.rs <- ecospat.niche.similarity.test(z1=grid.sum,z2=grid.res.s,rep=100,one.sided=F)
  #sim.w.rw <- ecospat.niche.similarity.test(z1=grid.wnt,z2=grid.res.w,rep=100,one.sided=F)

  t2 <- Sys.time()
  print(paste(i,t2-t1)) 
  
  x <- c(species=paste(i),centroid.distance=centroid.distance,I.obs=sim.s.w$obs$I,I.res.s=overlap.res.s$I,
         I.res.w=overlap.res.w$I,p.similar=sim.s.w$p.I) #collect output
  out[nrow(out)+1,] <- x #store in a df
}

###### summary tables for niche overlap & similarity
#overlap <- c(lapply(out,function(e) e$obs$D))
#similarity <- c(lapply(out,function(e) e$p.D))
#results <- data.frame(species=species,overlap=unlist(overlap),similarity=unlist(similarity))

#classify species migratory status by range map layers (there is probably a better way to do this...)
#seasonal status numbers are: 1=resident, 2=breeding, 3=nonbreeding 
r1 <- lapply(ranges,function(e) unlist(factor(e@data$SEASONAL))) #list of species and range map classifications
migStatus <- unlist(lapply(r1,function(e) {
  if(1 %in% levels(e)==T & 2 %in% levels(e)==T & 3 %in% levels(e)==F){ #if there's a resident and breeding range, partial migrant
    2 #2=partial migrant
  } else if(1 %in% levels(e)==T & 2 %in% levels(e)==F & 3 %in% levels(e)==T){ #resident and nonbreeding=partial migrant
    2
  } else if(1 %in% levels(e)==T & 2 %in% levels(e)==T & 3 %in% levels(e)==T){ #all three=partial migrant
    2
  } else if(1 %in% levels(e)==T & 2 %in% levels(e)==F & 3 %in% levels(e)==F) { #resident only=resident
    1 #1=resident
  } else if(1 %in% levels(e)==F & 2 %in% levels(e)==T & 3 %in% levels(e)==T) { #breeding & nonbreeding but NO resident = migratory
    3 #3=obligate migrant
  }
}))

migStatus <- data.frame(species=names(migStatus),migStatus)

results <- na.omit(merge(hummers,migStatus,by="species",all.x=T,all.y=T))
results[,2:6] <- as.numeric(as.character(unlist(results[,2:6])))
results$migStatus <- factor(results$migStatus)

ggplot(data=results,aes(x=centroid.distance,y=I.obs))+theme_bw()+
  geom_point()

ggplot(data=results,aes(x=I.obs,y=I.res.w,col=migStatus))+theme_bw()+xlim(0.25,1)+ylim(0.25,1)+
  geom_point()+geom_abline(slope=1)


hummers <- data.frame(fread("~/R/nicheTracker/nicheTracker2_hummingbirds.txt",drop=1))
names(hummers) <- c("species","centroid.distance", "I.obs", "I.res.s", "I.res.w", "p.similar")






# ### analysis as a function
# nicheTracker <- function(k) {
# 
#   #load species occurrence points, split by season
#   loc <- subset(gbif,species == k)
#   loc.sum <- subset(loc,month%in%c(5,6,7)==T)
#   loc.wnt <- subset(loc,month%in%c(11,12,1)==T)
#   loc.sum <- SpatialPoints(data.frame(loc.sum[,.(decimallongitude,decimallatitude)]), proj4string=CRS(proj4string(alt)))
#   loc.wnt <- SpatialPoints(data.frame(loc.wnt[,.(decimallongitude,decimallatitude)]), proj4string=CRS(proj4string(alt)))
# 
#   #find centroid of seasonal occurrences
#   centroid.sum <- gCentroid(loc.sum)
#   centroid.wnt <- gCentroid(loc.wnt)
#   centroid.distance <- pointDistance(centroid.sum,centroid.wnt,lonlat=T)
# 
#   #extract background data
#   sp.sum <- na.omit(extract(bg.sum.r,loc.sum))
#   sp.wnt <- na.omit(extract(bg.wnt.r,loc.wnt))
#   sp.res.s <- na.omit(extract(bg.wnt.r,loc.sum)) #resident on summer territory
#   sp.res.w <- na.omit(extract(bg.sum.r,loc.wnt)) #resident on winter territory
# 
#   #combo dataset for pca
#   df <- rbind(sp.sum,sp.wnt,sp.res.s,sp.res.w,bg.sum,bg.wnt)
# 
#   #set weighting to run PCA on full background data
#   weights <- c(rep(0,nrow(sp.sum)),rep(0,nrow(sp.wnt)),rep(0,nrow(sp.res.s)),
#                rep(0,nrow(sp.res.w)),rep(1,nrow(bg.sum)),rep(1,nrow(bg.wnt)))
# 
#   #run pca
#   pca <- dudi.pca(df,row.w=weights,nf=2,scannf=F,center=T,scale=T)
# 
#   #pull rows for species & background pc coords
#   pc.sum <- pca$li[1 : nrow(sp.sum),]
#   pc.wnt <- pca$li[(1+nrow(sp.sum)):(nrow(sp.sum)+nrow(sp.wnt)),]
#   pc.res.s <- pca$li[(1+nrow(sp.sum)+nrow(sp.wnt)) : (1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)),]
#   pc.res.w <- pca$li[(1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)) : (1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)+nrow(sp.res.w)),]
#   pc.bg.sum <- pca$li[(1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)+nrow(sp.res.w)) : (1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)+nrow(sp.res.w)+nrow(bg.sum)),]
#   pc.bg.wnt <- pca$li[(1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)+nrow(sp.res.w)+nrow(bg.sum)) : nrow(pca$li),]
#   pc.bg.all <- pca$li[(1+nrow(sp.sum)+nrow(sp.wnt)+nrow(sp.res.s)+nrow(sp.res.w)) : nrow(pca$li),]
# 
#   #grid with kernal density estimator (see Broenniman et al. 2012, J. Biogeography)
#   grid.sum <- ecospat.grid.clim.dyn(glob=pc.bg.all,glob1=pc.bg.sum,sp=pc.sum,R=100)
#   grid.wnt <- ecospat.grid.clim.dyn(glob=pc.bg.all,glob1=pc.bg.wnt,sp=pc.wnt,R=100)
#   grid.res.s <- ecospat.grid.clim.dyn(glob=pc.bg.all,glob1=pc.bg.wnt,sp=pc.res.s,R=100)
#   grid.res.w <- ecospat.grid.clim.dyn(glob=pc.bg.all,glob1=pc.bg.sum,sp=pc.res.w,R=100)
# 
#   #run similarity test
#   sim.s.w <- ecospat.niche.similarity.test(z1=grid.sum,z2=grid.wnt,rep=100,one.sided=F)
#   overlap.res.s <- ecospat.niche.overlap(grid.sum,grid.res.s,cor=T)
#   overlap.res.w <- ecospat.niche.overlap(grid.wnt,grid.res.w,cor=T)
#   #sim.s.rs <- ecospat.niche.similarity.test(z1=grid.sum,z2=grid.res.s,rep=100,one.sided=F)
#   #sim.w.rw <- ecospat.niche.similarity.test(z1=grid.wnt,z2=grid.res.w,rep=100,one.sided=F)
# 
#   #t2 <- Sys.time()
#   #print(paste(i,t2-t1))
# 
#   c(species=paste(i),centroid.distance=centroid.distance,I.obs=sim.s.w$obs$I,I.res.s=overlap.res.s$I,
#     I.res.w=overlap.res.w$I,p.similar=sim.s.w$p.I) #output
# }
# 
# out <- pblapply(c("Selasphorus rufus","Calypte anna"),function(e) nicheTracker(e))




