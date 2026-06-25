
### Eavaporation from surface water
# Penman equation
# arguments: mean air temp [deg C], vpd [Pa], net radiation [W m-2], elevation [m]
evap_calc<-function(tair, vpd, netrad, elev){

  #parameters
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

  # INTERMEDIATE EQUATIONS
  h = 287/9.81*((tair)+0.5*elev*LAPSE_PM)  #scale height in the atmosphere [m?]
  pz = PS_PM*exp(-elev/h)  #surface air pressure [Pa]
  lv = 2501000-2361*tair  #latent heat of vaporization [J Kg-1]
  gamma = 1628.6*pz/lv  #psychrometric constant [Pa C-1]
  r_air = 0.003486*pz/(275+tair)  #air density [Kg m-3]
  rs = 0  #minimal stomatal resistance [s m-1]
  rc = 0  #
  ra = log((2+(1/0.63-1)*d_Lower)/Z0_Lower)*log((2+(1/0.63-1)*d_Lower)/(0.1*Z0_Lower))/K2 #aerodynamic resistance [s m-1]
  rarc = 0  #architectural resistance [s m-1]
  svp = (6.11*exp((2.5*10^6/461.52)*((1/273.15)-(1/(tair+273.15)))))*100  # saturated vapor pressure, [Pa]
  slope =((B_SVP*C_SVP)/(C_SVP+tair)*(C_SVP+tair))*svp  #slope of saturated vapor pressure curve [Pa K-1]

  evap = ((slope*netrad+r_air*CP_PM*vpd/ra)/(lv*(slope+gamma*(1+(rc+rarc)/ra)))*SEC_PER_DAY)/10 #evaporation [cm day-1]

  return(evap)
}
