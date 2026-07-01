#function for calculating lake evaporation

evap_calc<-function(tmean, svp, vp, vpd, shortwave, longwave){
  #evap_calc<-function(elev,tair,shortwave,longwave,SH){
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
  svp=svp#saturated vapor pressure [Pa]
  vp=vp#vapor pressure [Pa]
  vpd=svp-vp  #vapor pressure deficit [Pa]
  slope=((B_SVP*C_SVP)/(C_SVP+tmean)*(C_SVP+tmean))*svp  #slope of saturated vapor pressure curve [Pa K-1]
  #shortwave=(0.75 + 2*10e-5 * elev)*ETrad #clear day solar radiation W m^2 day^-1 ... eqn 37 http://www.fao.org/3/X0490E/x0490e07.htm#air%20temperature
  net_short = (1-H2O_SURF_ALBEDO)*shortwave  # [W m-2]
  # From VIC func_surf_energy_bal.c
  #Tmp = Ts + KELVIN; Ts is soil temperature or surface temperature and KELVIN=273.15
  #LongBareOut = STEFAN_B * Tmp * Tmp * Tmp * Tmp;
  #a lake study (Binyamin et al. 2006; Int. J. Climatol. 26: 2261-2273) used E*sigma*T^4 as outgoing and E*longwave as incoming, where
  # E is emissivity of the water surface (0.97), sigma is the Stefan-Boltzman constant, and T is the
  # surface water temperature
  net_long= STEFAN_B *(tmean)*0.34-0.14*sqrt(vp)*(1.35*1-0.35) #W m-2
  rad=net_short+net_long  #[W m-2]
  evap=((slope*rad+r_air*CP_PM*vpd/ra)/(lv*(slope+gamma*(1+(rc+rarc)/ra)))*SEC_PER_DAY)/10 #evaporation [cm day-1]
  
  return(evap)
}
