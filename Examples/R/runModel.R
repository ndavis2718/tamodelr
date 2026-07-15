### example code for model simulations - using OSBS site from NEON w/ SUGG lake site

setwd("/Users/nicodavis/Desktop/R/TERC/tamodelr") # set working directory to tamodelr folder

library(deSolve) # load differential equation solver from library

# load functions
source("Examples/R/GddCalc.R") # calculate growing degree days
source("Examples/R/TAMmodel.R") # model code
source("Examples/R/addLeaves.R") # calculate leaf on/ leaf off

forcings<-read.csv("Examples/Data/OSBS_neonForcings_ex.csv") # read in forcing data
  #there are a bunch of variables from this file that we don't need right now; let's keep only what we need
  site_forcings<-data.frame(cbind(TIMESTAMP=forcings$TIMESTAMP, TA_F=forcings$TA_F, forcings[,106:113]))
  PFTtable<-read.csv("Examples/Data/paramTable.csv")# read in parameter table
 #vector of plant functional types at OSBS neon site
PFTtable$"EGNE1" <- PFTtable$"EGNE"#add new EGNE- Same PFT diff lake properties
PFTtable$"EGNE2" <- PFTtable$"EGNE"
pfts_site<-c("EGNE","EGNE1","EGNE2")

#### we need to find and assign the lake-specific parameters for the SUGG neon site
# Aa = surface area (m^2), Ac = catchment area (m^2), zbar = mean depth (m)
zbar = 2.06
Ac = 3.96 * 10^7
Aa = 730000
#the lake qualities are encoded in the plant types; every PFT has the same lake depth, SA, etc
PFTtable[which(PFTtable$pName=="zbar"), 2:8]<-zbar
PFTtable[which(PFTtable$pName=="Ac"), 2:8]<-Ac
PFTtable[which(PFTtable$pName=="Aa"), 2:8]<-Aa
#where the row of col pNmae = "zbar" and where colnames == EGNE1)
PFTtable[which(PFTtable$pName == "zbar"), which(colnames(PFTtable) == "EGNE1")] <- zbar*2
PFTtable[which(PFTtable$pName == "Ac"), which(colnames(PFTtable) == "EGNE1")] <-  Ac*2
PFTtable[which(PFTtable$pName == "Aa"), which(colnames(PFTtable) == "EGNE1")] <- Aa*2
#where the row of col pNmae = "zbar" and where colnames == EGNE2)
PFTtable[which(PFTtable$pName == "zbar"), which(colnames(PFTtable) == "EGNE2")] <- zbar*1/2
PFTtable[which(PFTtable$pName == "Ac"), which(colnames(PFTtable) == "EGNE2")] <-  Ac *1/2
PFTtable[which(PFTtable$pName == "Aa"), which(colnames(PFTtable) == "EGNE2")] <- Aa*1/2


### STOP: check forcing data to make sure that there are no negative values (negative temperatures are ok but everything else should be positive)

#1. Model spinup OR input own inital conditions if data exist for site.
##run equilibrium to get initial conditions | equilibrium means a linear regression is about 0 for each pool (things arent changing)
#create a data frame that holds initial conditions for each pft
colnames <-c("Cw","Cl", "Cs1", "Cs2", "Cs3", "Cs4", "Cdoc1", "Cdoc2", "W1", "W2", "Ca", "Cr", "Ccwd", "Cdic1", "Cdic2", "W3", "Cdic3", "Cdoc3", "Calat", "PFT")
initDF<-data.frame(matrix(ncol=length(colnames), nrow=length(pfts_site))) #initial conditions setup. 1 col for each item we track (types of carbon)
colnames(initDF)<-colnames
initDF$PFT<-pfts_site

#data.frame to store the initial conditions for dynamic simulations
inits<-data.frame(matrix(ncol=length(colnames), nrow=length(pfts_site)))
colnames(inits)<-colnames
inits$PFT<-pfts_site

#Store linear regressions that will be used to confirm equilibrium at the end of spinup
regsDF<-data.frame(matrix(ncol=length(colnames), nrow=length(pfts_site)))
colnames(regsDF)<-colnames
regsDF$PFT<-pfts_site

#Fill in the initial conditions data.frame with generic pft-specific values prior to model spinup start
#in PP version ive mad calat init guess 10 for each PFT
for(y in 1:nrow(initDF)){
  params<-as.numeric(PFTtable[,which(colnames(PFTtable)==initDF$PFT[y])])  #params for each pft. goes to pft[y] in the pft table, extracts table as pure list of numbers
  names(params)<-PFTtable$pName #name params
  if(initDF$PFT[y]=="EGBR"){
    initDF[y,1:(ncol(initDF)-1)]<-c(13000, as.numeric(params[[55]]*params[[20]]*params[[21]]), 60, 60, 300, 2000, 50, 30, as.numeric(params[[19]]*.90), as.numeric(params[[41]]*.90), 600000, 3000, 50, 0.01, 0.01, as.numeric(params[[41]]*.90), 0.0001, 10, 200000000) #fill initial conditions
  } else if(initDF$PFT[y]=="EGNE"){
    initDF[y,1:(ncol(initDF)-1)]<-c(13000, as.numeric(params[[55]]*params[[20]]*params[[21]]), 60,60, 300, 2000, 50, 30, as.numeric(params[[19]]*.90), as.numeric(params[[41]]*.90), 600000, 3000, 50, 0.01, 0.01,  as.numeric(params[[41]]*.90), 0.0001, 10, 200000000) #fill initial conditions
  } else if(initDF$PFT[y]=="DEBR"){
    initDF[y,1:(ncol(initDF)-1)]<-c(13000, as.numeric(params[[55]]*params[[20]]*params[[21]]), 60, 60, 300, 2000, 50, 30, as.numeric(params[[19]]*.90), as.numeric(params[[41]]*.90), 600000, 3000, 50, 0.01, 0.01,  as.numeric(params[[41]]*.90), 0.01, 10, 200000000) #fill initial conditions
  } else if(initDF$PFT[y]=="DENE"){
    initDF[y,1:(ncol(initDF)-1)]<-c(13000, as.numeric(params[[55]]*params[[20]]*params[[21]]), 60,60, 300, 2000, 50, 30, as.numeric(params[[19]]*.90), as.numeric(params[[41]]*.90), 600000, 3000, 50, 0.01, 0.01,  as.numeric(params[[41]]*.90), 0.01, 10, 200000000) #fill initial conditions
  } else if(initDF$PFT[y]=="SH"){
    initDF[y,1:(ncol(initDF)-1)]<-c(10000, as.numeric(params[[55]]*params[[20]]*params[[21]]), 60,60, 300, 2000, 50, 30, as.numeric(params[[19]]*.90), as.numeric(params[[41]]*.90), 600000, 2500, 50, 0.01, 0.01,  as.numeric(params[[41]]*.90), 0.01, 10, 200000000) #fill initial conditions
  } else if(initDF$PFT[y]=="GR" | initDF$PFT[y]=="CR"){
    initDF[y,1:(ncol(initDF)-1)]<-c(0, as.numeric(params[[55]]*params[[20]]*params[[21]]), 60,60, 300, 2000, 50, 30, as.numeric(params[[19]]*.90), as.numeric(params[[41]]*.90), 600000, 500, 50, 0.01, 0.01,  as.numeric(params[[41]]*.90), 0.01, 10, 200000000) #fill initial conditions
  }
}

#4a. Loop over each plant functional type
for(i in 1:nrow(initDF)){
  nruns=0 #keep track of how many times spinup simulations are done for each pft
  repeat{
    params<-PFTtable[,which(colnames(PFTtable)==initDF$PFT[i])]  #params for each pft
    names(params)<-c(PFTtable$pName) #name params

    site_forcings$year<-substr(site_forcings$TIMESTAMP, 1,4) # add a year column
    #site forcings makes unuique forcings every day of the year for x years

    # needs to start on a january first for generating leaf on / leaf off days properly
    site_forcings$mon_day<-paste(substr(site_forcings$TIMESTAMP, 5,6), substr(site_forcings$TIMESTAMP, 7,8), sep="-")
    if(site_forcings$mon_day[1] != "01-01"){
      site_forcings<-site_forcings[which(site_forcings$mon_day== "01-01")[1]:nrow(site_forcings), ]#ensuring we span from J1 to the end. J1 must be first
    }

    # create the leaf phenology from temperature (growing degree days determine when leaves pop up)
    site_forcings<-addLeaves(df=site_forcings, a=params[[57]] , k=params[[56]], b=params[[58]])
    site_forcings$runDay<-1:nrow(site_forcings) # time steps for the model runs

    spin_times<-1:nrow(site_forcings) #list of time-steps in days; should match forcing data inputs

    length(which(is.na(site_forcings))) #make sure that there are no missing data (as NAs)

     #Set initial conditions from the data.frame created in step #4
    S0<<-c(Cw=initDF[i,1],Cl=initDF[i,2],Cs1=initDF[i,3],Cs2=initDF[i,4], Cs3=initDF[i,5],
           Cs4=initDF[i,6], Cdoc1=initDF[i,7], Cdoc2=initDF[i,8], W1=initDF[i,9],
           W2=initDF[i,10], Ca=initDF[i,11], Cr=initDF[i,12], Ccwd=initDF[i,13],
           Cdic1=initDF[i,14], Cdic2=initDF[i,15], W3=initDF[i,16], Cdic3=initDF[i,17],
           Cdoc3=initDF[i,18], Ci = initDF[i,19])

    #define forcing approx functions #interpolate smoothly through time for ode
    PARapprox<<-approxfun(x=as.numeric(site_forcings$runDay), y = as.numeric(site_forcings$PAR_e))
    Papprox<<-approxfun(x =as.numeric(site_forcings$runDay), y = as.numeric(site_forcings$Precip))
    VPDapprox<<-approxfun(x=as.numeric(site_forcings$runDay), y = as.numeric(site_forcings$VPD_kPa))
    Tair_approx<<-approxfun(x=as.numeric(site_forcings$runDay), y = as.numeric(site_forcings$TA_F))
    Tsoil_approx<<-approxfun(x=as.numeric(site_forcings$runDay), y = as.numeric(site_forcings$tsoil))
    Evap_approx<<-approxfun(x=as.numeric(site_forcings$runDay), y= as.numeric(site_forcings$aqEvap))
    Snowapprox<<-approxfun(x=as.numeric(site_forcings$runDay), y=as.numeric(rep(0, times=nrow(site_forcings))))
    CO2_approx<<-approxfun(x=as.numeric(site_forcings$runDay), y=as.numeric(rep(params[[53]], times=nrow(site_forcings))))
    gdd_approx<<-approxfun(x=as.numeric(site_forcings$runDay), y=as.numeric(site_forcings$GDD))
    Doy_approx<<-approxfun(x=as.numeric(site_forcings$runDay), y=as.numeric(site_forcings$DOY))
    sen_approx<<-approxfun(x=as.numeric(site_forcings$runDay), y=as.numeric(site_forcings$sen))
    grow_approx<<-approxfun(x=as.numeric(site_forcings$runDay), y=as.numeric(site_forcings$grow))
    green_approx<<-approxfun(x=as.numeric(site_forcings$runDay), y=as.numeric(site_forcings$greenup))

    spinup=ode(y=S0,times=spin_times,func=PPStep, parms=params, method="euler") #run the model!
    spinup=data.frame(spinup) 

    nruns<-nruns+1 #continue keeping track of the number of runs until equilibrium is achieved..

    #update initDF with mean pool size from last year, so that each run is closer to equilibrium
    for(z in 1:(ncol(initDF)-2)){
      initDF[i,z]<-round(spinup[nrow(spinup),which(colnames(spinup)==colnames(regsDF)[z])], digits=3)
    }
    initDF[i,9:10]<-c(as.numeric(params[[19]]*.90), as.numeric(params[[41]]*.90))
    #check for equilibrium w/linear regression
    for(x in 1:(ncol(regsDF)-2)){
      regMod<-lm(spinup[,which(colnames(spinup)==colnames(regsDF)[x])]~spinup$time)
      regsDF[i,x]<-regMod$coefficients[[2]]
    }

    tol = .05
    if((((length(which(regsDF[i,1:(ncol(regsDF)-2)] > -1*tol))==(ncol(regsDF)-2))) & (length(which(regsDF[i,1:(ncol(regsDF)-2)]< tol))==(ncol(regsDF)-2))) | nruns==20){
      break
    } #if the slope of linear regressions is close enough to 0, stop repeating the loop
  }

  spinup$pft<-initDF$PFT[i] #label PFT

  #update the initial conditions for dynamic runs with the model outputs from the last day of spin-up
  #we can ask: how do these final initial conditions change when we adjust values about the lake? 
  for(z in 1:(ncol(inits)-1)){ #i = PFT =row
    inits[i,z]<-spinup[nrow(spinup),which(colnames(spinup)==colnames(inits)[z])]
  }
}#end model loop



