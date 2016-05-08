
minnlocs <- function(i) {
  loc <- subset(gbif,species == i) 
  loc.sum <- subset(loc,month%in%c(6,7,8)==T)
  loc.wnt <- subset(loc,month%in%c(12,1,2)==T)
  loc.sum <- SpatialPoints(data.frame(loc.sum[,.(decimallongitude,decimallatitude)]), proj4string=CRS(proj4string(ranges[[i]])))
  loc.wnt <- SpatialPoints(data.frame(loc.wnt[,.(decimallongitude,decimallatitude)]), proj4string=CRS(proj4string(ranges[[i]])))

  #crop to range if range map is present and species is not an austral migrant
  if( !is.null(simple.ranges[[i]])){  
    if(ymin(gCentroid(simple.ranges[[i]][simple.ranges[[i]]@data$layer %in% c(1,2),])) >= ymin(gCentroid(simple.ranges[[i]][simple.ranges[[i]]@data$layer %in% c(1,3),])) ){ 
      loc.sum <- gIntersection(loc.sum,simple.ranges[[i]][simple.ranges[[i]]@data$layer %in% c(1,2),])
      loc.wnt <- gIntersection(loc.wnt,simple.ranges[[i]][simple.ranges[[i]]@data$layer %in% c(1,3),])
    }
  }
  min(length(loc.sum),length(loc.wnt))
}