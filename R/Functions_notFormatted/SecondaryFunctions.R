####Functions required to run TAM####
#By: Ceara J. Talbot

library(lubridate)
library(parallel)
library(reshape2)
library(deSolve)
library(zoo)

###formatting neon forcings
summarize_data<-function(data_in, neon_var, val_col){
  
  if(neon_var=="Triple Aspirated Air Temperature" | neon_var=="Relative Humidity" | neon_var=="Soil Temperature" | neon_var=="Barometric Pressure" | neon_var=="Photosynthetically Active Radiation" | neon_var=="Net Radiation"){
    data_in$Date<-substr(data_in$startDateTime, 1, 10)
    data_out<-aggregate(data_in[,val_col]~data_in$Date+data_in$siteID, FUN=mean)
    colnames(data_out)[1:3]<-c("date", "siteID", "value")
    
  } else if(neon_var=="Precipitation"){
    data_in$Date<-substr(data_in$startDateTime, 1, 10)
    data_out<-aggregate(data_in[,val_col]~data_in$Date+data_in$siteID, FUN=sum)
    colnames(data_out)[1:3]<-c("date", "siteID", "value")
    
  }else if(neon_var=="SW Chemistry"){
    data_in$Date<-substr(data_in$startDate, 1, 10)
    data_out<-data_in[which(data_in$analyte=="DOC" | data_in$analyte=="DIC" | data_in$analyte=="TOC"),]
    data_out<-aggregate(data_out, analyteConcentration~Date+siteID+analyte, FUN=mean)
    #data_out_sd<-aggregate(data_out$analyteConcentration~data_out$Date+data_out$siteID+data_out$analyte, FUN=sd)
    #data_out_count<-aggregate(data_out$analyteConcentration~data_out$Date+data_out$siteID+data_out$analyte, FUN=length)
    #data_out<-data.frame(cbind(data_out_me, data_out$analyteUnits, data_out$belowDetectionQF))
    #colnames(data_out)[4:6]<-c("analyte", "units", "BDL")
  }else if(neon_var=="Continuous Q"){
    data_in$Date<-as.Date(data_in$endDate, format="%Y-%mm-%dd")
    data_out<-aggregate(data_in, continuousDischarge~Date+siteID, FUN=mean)
  }else if(neon_var=="Field Q"){
    data_in$Date<-as.Date(data_in$collectDate, format="%Y-%mm-%dd")
  }
  
  return(data_out)
}

###fill in missing data in forcings
fill_forcings<-function(data_in){
  
  data_in$DOY<-yday(data_in$date)
  for(i in 3:ncol(data_in)){
    for(r in 1:nrow(data_in)){
      if(is.na(data_in[r,i])){
        data_in[r,i]<-mean(data_in[which(data_in$DOY==data_in$DOY[r]),i], na.rm=T)
      }
    }
  }
  return(data_in)
}

###Leaf phenology functions
#Growing degree day calculations
  #input is temps vector for a single year
  GddCalc<-function(temps, ths, chill=FALSE){
    outTemps<-c(ifelse(temps[1] > ths, temps[1], 0), rep(0, (length(temps)-1))) #setup the initial temps
    if(chill==TRUE){
      for(i in 1:length(temps)){
        if(temps[i] < ths & temps[i] > 0){
          outTemps[i]<-1 #yes a chilling day
        } else{
          outTemps[i]<-0 #not a chilling day
        }
      }
    }else{
      for(i in 2:length(temps)){
        if(temps[i] > ths){
          outTemps[i]<-outTemps[(i-1)]+temps[i]
        } else{
          outTemps[i]<-outTemps[(i-1)]
        }
      }
    }
    return(outTemps)
  } #end function

#Leaf on and leaf off date
  ####Leaf module setup, based on CABLE phenology Haverd et al. 2018
  #need to make a function that generates a timeseries output that can be used
  #from phenology (similar to previous TAM; Talbot et al. 2022)
  
  addLeaves<-function(df, a, k, b, glength){
    df$GDD<-0 #growing degree days
    df$DOY<-0 #day of year
    df$GDDc<-0
    df$GDD0<-0
    df$grow<-0
    df$greenup<-0
    df$sen<-0 #10-day cumulative degree days (for leaf off)
    for(i in 1:length(unique(df$year))){
      sub<-df[df$year==unique(df$year)[i], ]
      df$GDD[df$year==unique(df$year)[i]]<-GddCalc(temps=sub$TA_F, ths=0, chill=FALSE)
      df$DOY[df$year==unique(df$year)[i]]<-1:nrow(sub)
      df$GDDc[df$year==unique(df$year)[i]]<-GddCalc(temps=sub$TA_F, ths=0, chill=TRUE)
      ##GDD0 following Sykes et al 1996: GDD0=a+be-kC; C = length of chilling period, ab&k are species specific 
      df$GDD0[df$year==unique(df$year)[i]]<-a+b*exp(-k*(length(which(df$GDDc[df$year==unique(df$year)[i]] == 1 & df$DOY[df$year==unique(df$year)[i]] < glength)))) #create relationship between this and chilling days
      df$greenup[df$year==unique(df$year)[i]]<-ifelse(df$GDD[df$year==unique(df$year)[i]]-df$GDD0[df$year==unique(df$year)[i]] < glength & df$GDD[df$year==unique(df$year)[i]] > df$GDD0[df$year==unique(df$year)[i]], 1, 0) #when gdd > gdd0 leaf on starts, adds 1s until GDD-GDD0 > 200 
      df$sen[df$year==unique(df$year)[i]]<-c(rep(0, (which(df$greenup[df$year==unique(df$year)[i]] > 0)[1]+glength)), rep(1, nrow(df[df$year==unique(df$year)[i], ])-(which(df$greenup[df$year==unique(df$year)[i]] > 0)[1]+glength))) #200 days after the onset of leaf growth, leaves begin to senesce 
      df$grow[df$year==unique(df$year)[i]]<-ifelse(df$GDD[df$year==unique(df$year)[i]] > df$GDD0[df$year==unique(df$year)[i]] & df$greenup[df$year==unique(df$year)[i]] == 0 & df$sen[df$year==unique(df$year)[i]] == 0, 1, 0) #when gdd > gdd0 leaf on starts, adds 1s until GDD-GDD0 > 200 
    }
    
    return(df)
  }
  
###Filling in missing climate drivers
  genWeather2<-function(climateData, soilk=5){
    ###soil temperature (surface)###
    #generate daily soil T w/trailing moving average
    #default=5 days
    climateData$tsoil<-rollmean(climateData$TA_F, k=soilk, fill=NA, align="right") 
    for(i in 1:(soilk-1)){
      if(i==1){
        climateData$tsoil[i]<-climateData$TA_F[i]
      } else{
        climateData$tsoil[i]<-mean(c(climateData$TA_F[1:i],climateData$TA_F[soilk]))
      }
    }
    ######PAR######
    #default coeff from Britton and Dodd 1976
    #1 W/m2 ≈ 4.6 μmole.m2/s, 1 μmole.m2/s = 1 μEinstein/m2/s
    #mstmip is J/m2.. 1 J = 1 W, 1 W/m2 ≈ 4.6 μmole.m2/s,
    #climateData$PARest<-climateData$swdown*4.6/1000000*parC #W/m2 to E/m2
    #####VPD####
    #helpful: https://cran.r-project.org/web/packages/humidity/vignettes/humidity-measures.html
    #hPa to kPa=0.1; Pa to kPa = 0.001
    climateData$svp<-(6.11*exp((2.5*10^6/461.52)*((1/273.15)-(1/(climateData$TA_F+273.15)))))*0.1  #kPa; CHECKED
    #climateData$vp<-climateData$RH_30min*climateData$svp*0.001 #kPa 
    #climateData$vpd<-climateData$svp-climateData$vp #kPa
    #loop evap calculation over each DOY
    #evap in mm day^-1
    climateData$aqEvap<-0
    climateData$netrad<-climateData$NETRAD #W m-2; CHECKED
    for(b in 1:nrow(climateData)){
      n<-evap_calc(tmean=(climateData$TA_F[b]),svp=(climateData$svp[b]*1000), #convert KPa to Pa
                   vpd=(climateData$VPD_F[b]*100), netrad=climateData$NETRAD[b]) #convert hPa to Pa for vpd; *100
      climateData$aqEvap[b]<-n
    }
    return(climateData)
  }
  
###Aquatic surface evaporation
  #Based on code provided by SEJ
  evap_calc<-function(tmean, svp, vpd, netrad){
    # REQUIRED INPUTS FOR THE LOCATION
    # elevation [elev; m]
    # air temperature [tair; C]
    # max air temp [tmax; K]
    # min air temp [tmin; k]
    # svp [kPa]
    # vp [kPa]
    # vpd [kPa]
    # daily extraterrestrial radiation [ETrad; MJ m^-2 day^-1]
    # incoming shortwave radiation [shortwave; ????]
    # incoming longwave radiation [longwave]
    # specific humidity [SH; Kg Kg-1]
    # OUTPUT
    # evaporation [evap; mm day-1]
    # PARAMETERS
    LAPSE_PM = -0.006  #environmental lapse rate [C m-1]
    PS_PM = 101300  #sea level air pressure [Pa]
    CP_PM = 1013  #specific heat of moist air at constant pressure [J kg-1 C-1]
    Z0_Lower = 0.001  #roughness
    d_Lower = 0.0054  #displacement
    von_K = 0.4  #Von Karman constant for evapotranspiration
    K2 = von_K*von_K
    C = 2.16679 #constant for specific humidity to vapor pressure conversion [g K J-1]
    A_SVP = 0.61078  #A term in saturated vapor pressure calculation
    B_SVP = 17.269  #B term in saturated vapor pressure calculation
    C_SVP = 237.3  #C term in saturated vapor pressure calculation
    SEC_PER_DAY = 86400  #seconds per day
    H2O_SURF_ALBEDO = 0.08  #albedo of water surface
    STEFAN_B = 5.6696e-8 #5.6696e-8, stefan-boltzmann constant [W/m^2/K^4]
    elev = 1
    # INTERMEDIATE EQUATIONS
    h=287/9.81*((tmean)+0.5*elev*LAPSE_PM)  #scale height in the atmosphere [m?]
    pz=PS_PM*exp(-elev/h)  #surface air pressure [Pa]
    #pz=pressure  #surface air pressure [Pa]
    lv=2501000-2361*tmean  #latent heat of vaporization [J Kg-1]
    gamma=1628.6*pz/lv  #psychrometric constant [Pa C-1]
    r_air=0.003486*pz/(275+tmean)  #air density [Kg m-3]
    rs=0  #minimal stomatal resistance [s m-1]
    rc=0  #
    ra=log((2+(1/0.63-1)*d_Lower)/Z0_Lower)*log((2+(1/0.63-1)*d_Lower)/(0.1*Z0_Lower))/K2 #aerodynamic resistance [s m-1]
    rarc=0  #architectural resistance [s m-1]
    #svp=svp#saturated vapor pressure [Pa]
    #vp=vp#vapor pressure [Pa]
    vpd=vpd  #vapor pressure deficit [Pa]
    slope=((B_SVP*C_SVP)/(C_SVP+tmean)*(C_SVP+tmean))*svp  #slope of saturated vapor pressure curve [Pa K-1]
    #shortwave=(0.75 + 2*10e-5 * elev)*ETrad #clear day solar radiation W m^2 day^-1 ... eqn 37 http://www.fao.org/3/X0490E/x0490e07.htm#air%20temperature
    #net_short = (1-H2O_SURF_ALBEDO)*shortwave  # [W m-2]
    # From VIC func_surf_energy_bal.c
    #Tmp = Ts + KELVIN; Ts is soil temperature or surface temperature and KELVIN=273.15
    #LongBareOut = STEFAN_B * Tmp * Tmp * Tmp * Tmp;
    #a lake study (Binyamin et al. 2006; Int. J. Climatol. 26: 2261-2273) used E*sigma*T^4 as outgoing and E*longwave as incoming, where
    # E is emissivity of the water surface (0.97), sigma is the Stefan-Boltzman constant, and T is the
    # surface water temperature
    #net_long= STEFAN_B *(tmean)*0.34-0.14*sqrt(vp)*(1.35*1-0.35) #W m-2
    rad= netrad #net_short+net_long  #[W m-2]
    evap=((slope*rad+r_air*CP_PM*vpd/ra)/(lv*(slope+gamma*(1+(rc+rarc)/ra)))*SEC_PER_DAY)/10 #evaporation [cm day-1]
    
    return(evap)
  }
  
  ##function to summarize outputs for multiple PFTs into a single daily val per site
  pft_summaries<-function(df, pfts_site){
    
    summ_df<-data.frame(matrix(nrow=length(unique(df$Date)), ncol=ncol(df)))
    colnames(summ_df)<-colnames(df)
    summ_df$Date<-unique(df$Date)
    
    for(i in 1:nrow(summ_df)){
      summ_df[i,1:85]<-colMeans(df[which(df$Date==summ_df$Date[i]), 1:85])
    }
    return(summ_df)
  }
  
  