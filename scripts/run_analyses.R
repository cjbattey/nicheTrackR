#### NicheTracker2: Seasonal movement and niche-tracking analysis 
#this script implements and summarizes results - see setup_data and nicheTracker for data prep and analysis, respectively
library(data.table);library(dismo);library(ecospat);library(pbapply);library(rgeos);library(ggplot2);library(doMC);library(foreach);library(plyr)
registerDoMC(cores=16) 
setwd("~/Dropbox/nicheTracker_ec2/")

#load data (warning: setup_data.R clears all objects from the workspace - see last line.)
source("scripts/setup_data_staticBG.R")
source("scripts/ecospat.niche.similarity.test.noplot.R")
source("scripts/nicheTracker_staticBG.R")

#staticBG runs same BG data for all species (all areas are "available" in similarity test)
out <- foreach(i=levels(gbif$species), .combine=rbind) %dopar% nicheTracker_staticBG(i)

#########################################
####### data summary and plots ##########
#########################################

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
results$I.diff <- results$I.obs - results$I.res.w

#plot centroid distance vs niche overlap
ggplot(data=results,aes(x=centroid.distance,y=I.obs))+theme_bw()+
  geom_point()+geom_smooth(method = "lm")

#plot observed versus winter-resident niche overlap, colored by migratory status
ggplot(data=results,aes(x=I.obs,y=I.res.w,col=migStatus))+theme_bw()+xlim(0.25,1)+ylim(0.25,1)+
  geom_point()+geom_abline(slope=1)

results.sum <- ddply(results,"migStatus",function(e) c(percent.similar=nrow(subset(e,e$p.similar < 0.05))/nrow(e),
                                                       number.species=nrow(e),
                                                       centroid.disance=mean(e$centroid.distance),
                                                       I.obs=mean(e$I.obs)))

results$migStatusbin <- gsub("3","2",results$migStatus)
ggplot(data=results,aes(x=migStatusbin,y=I.diff))+theme_bw()+
  ggtitle("Observed - Resident Niche Overlap\nParulidae")+
  geom_boxplot(outlier.color=NA)+
  geom_point(alpha=0.3,col="red",size=1)

write.csv(results,"./out/Icteridae_staticBG_1occPerCell.csv",row.names = F)
