## vertnet datasets ##
setwd("./data/vertnet")
files <- list.files()
species <- gsub(".tsv","",files)
species <- gsub(".txt","",species)
vert <- lapply(files,FUN=function(e) read.delim(e,quote="")) #read in 
fields <- c("institutioncode","catalognumber","sex","lifestage","reproductivecondition","oganismname","year","month","day", 
            "country","stateprovince","county","locality","decimallatitude","decimallongitude","coordinateuncertaintyinmeters",
            "scientificname","genus","specificepithet","infraspecificepiphet")  #fields I want
vert <- lapply(vert,FUN=function(e) e[,which(names(e) %in% fields == T)] )  #subset fields for each speces
b <- list()

for(i in 1:length(vert)){  #dumb loop to reclassify "months" and pull out breeding/wintering months
  a <- vert[[i]]
  a$month <- as.numeric(as.character(a$month))
  occ.brd[[i]] <- subset(a,month %in% brd.months[[i]] ==T)
  
  #attempting to georeference specimens w/o coordinates via google. Will need to repeat with wintering reports. 
  Coords <- occ.brd[[i]][which(is.na(occ.brd[[i]]$decimallongitude) == T),]
  noCoords <- occ.brd[[i]][which(is.na(occ.brd[[i]]$decimallongitude) == T),]
  noCoords$georef.loc <- paste(noCoords$locality,", ",noCoords$county," county",", ",noCoords$stateprovince,", ",noCoords$country, sep="")
  latlong <- geocode(noCoords$georef.loc,oneRecord = T)
  
  
  occ.wnt[[i]] <- subset(a,month %in% wnt.months[[i]] ==T)
}

n.brd.occs <- lapply(occ.brd,function(e) nrow(e)) #check number of reports in each season
n.wnt.occs <- lapply(occ.wnt,function(e) nrow(e))

