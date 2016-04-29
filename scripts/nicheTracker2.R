#### NicheTracker2: Seasonal movement and niche-tracking analysis 
#this script implements and summarizes results - see setup_data and nicheTracker for data prep and analysis, respectively
library(data.table);library(dismo);library(ecospat);library(pbapply);library(rgeos);library(ggplot2);library(doMC);library(foreach)
registerDoMC(cores=2) 

#load data
source("scripts/setup_data.R")
source("scripts/ecospat.niche.similarity.test.noplot.R")
source("scripts/nicheTracker.R")

#run analyses in parallel (~6gb memory per core; swap ok with solid-state drive on mbp)
out <- foreach(i=levels(gbif$species), .combine=rbind) %dopar% nicheTracker(i)

#add migratory status from range maps
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

results <- na.omit(merge(out,migStatus,by="species",all.x=T,all.y=T))
results[,2:6] <- as.numeric(as.character(unlist(results[,2:6])))
results$migStatus <- factor(results$migStatus)

#plot centroid distance vs niche overlap
ggplot(data=m,aes(x=centroid.distance,y=I.obs,col=migStatus))+theme_bw()+
  geom_point()

#plot observed versus winter-resident niche overlap, colored by migratory status
ggplot(data=results,aes(x=I.obs,y=I.res.w,col=migStatus))+theme_bw()+xlim(0.25,1)+ylim(0.25,1)+
  geom_point()+geom_abline(slope=1)

m <- subset(results,migStatus==3 | migStatus == 2)
nm <- subset(results,migStatus==1)

nrow(subset(m,m$p.similar < 0.05))/nrow(m)
nrow(subset(nm,nm$p.similar < 0.05))/nrow(nm)

write.csv(results,"/R/nicheTracker/out/vireonidae_29Apr16.csv",row.names = F)
