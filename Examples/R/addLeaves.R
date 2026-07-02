####Leaf module setup, based on CABLE phenology Haverd et al. 2018
#need to make a function that generates a timeseries output that can be used
#fro phenology (similar to previous TAM; Talbot et al. 2022)

addLeaves<-function(df, a, k, b){
  df$GDD<-0 #growing degree days
  df$DOY<-0 #day of year
  df$GDDc<-0
  df$GDD0<-0
  df$grow<-0
  df$greenup<-0
  df$sen<-0 #10-day cumulative degree days (for leaf off)
  for(i in 1:length(unique(df$Year))){
    sub<-df[df$Year==unique(df$Year)[i], ]
    df$GDD[df$Year==unique(df$Year)[i]]<-GddCalc(sub$TA_F, ths=0, chill=FALSE)
    df$DOY[df$Year==unique(df$Year)[i]]<-1:nrow(sub)
    df$GDDc[df$Year==unique(df$Year)[i]]<-GddCalc(sub$TA_F, ths=0, chill=TRUE)
    ##GDD0 following Sykes et al 1996: GDD0=a+be-kC; C = length of chilling period, ab&k are species specific 
    df$GDD0[df$Year==unique(df$Year)[i]]<-a+b*exp(-k*(length(which(df$GDDc[df$Year==unique(df$Year)[i]] == 1 & df$DOY[df$Year==unique(df$Year)[i]] < 200)))) #create relationship between this and chilling days
    df$greenup[df$Year==unique(df$Year)[i]]<-ifelse(df$GDD[df$Year==unique(df$Year)[i]]-df$GDD0[df$Year==unique(df$Year)[i]] < 200 & df$GDD[df$Year==unique(df$Year)[i]] > df$GDD0[df$Year==unique(df$Year)[i]], 1, 0) #when gdd > gdd0 leaf on starts, adds 1s until GDD-GDD0 > 200 
    df$sen[df$Year==unique(df$Year)[i]]<-c(rep(0, (which(df$greenup[df$Year==unique(df$Year)[i]] > 0)[1]+200)), rep(1, nrow(df[df$Year==unique(df$Year)[i], ])-(which(df$greenup[df$Year==unique(df$Year)[i]] > 0)[1]+200))) #200 days after the onset of leaf growth, leaves begin to senesce 
    df$grow[df$Year==unique(df$Year)[i]]<-ifelse(df$GDD[df$Year==unique(df$Year)[i]] > df$GDD0[df$Year==unique(df$Year)[i]] & df$greenup[df$Year==unique(df$Year)[i]] == 0 & df$sen[df$Year==unique(df$Year)[i]] == 0, 1, 0) #when gdd > gdd0 leaf on starts, adds 1s until GDD-GDD0 > 200 
    
     #need to turn these into timeseries for NPP allocation
    }
  
  return(df)
}
