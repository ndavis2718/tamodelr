####Calculate ET radiation and convert to PAR given latitude and DOY####
library(sirad) 

calcPAR<-function(df, lat){ #df should have a DOY column
  
  pi<-3.14159265 #pi
  lat<-as.numeric(lat) #latitude
  lat_radian<-lat * pi/180   #decimal deg to radians
  
  n<-extrat(df$DOY, lat_radian)
  df$PARest<-n$ExtraTerrestrialSolarRadiationDaily*4.6*0.45 #units are MJm^-2 to to Em^-2 and converted to PAR with *0.45 after Britton & Dodd 1993
  df$ETrad_W<-n$ExtraTerrestrialSolarRadiationDaily/0.0864 #convert to W m^-2 for calculating aquatic evaporation

  return(df)
}

