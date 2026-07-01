###get midpoint of huc#####

getMidLat<-function(huc, wbdDir){
  library(sf)
  source("functions/fix8.R")
  #my wbdDir="/Users/cearatalbot/HUC04Example/myWBD/"
  
  fc <- sf::st_read(paste(wbdDir, "WBD_National_GDB.gdb", sep=""),layer="WBDHU8")
  fc<-fc[fc$huc8==huc,]
  fc_geom<-st_geometry(fc)
  fcNP<-st_transform(fc_geom,crs="EPSG:4326") #correct 3857 when needing projection
  
  myHUC = st_as_sf(data.frame(HUC8=as.numeric(fc$huc8), geometry=fcNP, crs ="+proj=longlat +datum=WGS84 +no_defs"))
  myHUC$HUC8<-fix8(myHUC$HUC8)
  
  midPoint = st_centroid(myHUC) #get centroid of HUC
  coords<-data.frame(matrix(ncol=3, nrow=1))
  colnames(coords)<-c("huc", "x", "y")
  coords$huc[1]<-huc
  coords$x[1]<-st_coordinates(midPoint)[1] #extract the x and y
  coords$y[1]<-st_coordinates(midPoint)[2] 
  
  return(coords)
}
