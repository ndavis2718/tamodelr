####setup TAM DA for each HUC#####

setupTAM<-function(huc){
  
      #paramTable
      pTable<-read.csv("paramTable.csv") #at some point need to change to vic params avg'd for the watershed
      pTable[37,2:8]<-c(0, 0.0027, 0.0027, 0,0.0027,0,0.0027)
      pTable[35,2:8]<-c(0.1, 0.1, 0.1, 0.1, 0, 0.1, 0)
      pTable[12,5]<-2
      pTable[15,2:8]<-c(rep(0.07, times= 7))  #old val is .004
      pTable[18,2:8]<-c(10.9, 10.9, 12.5, 12.5, 12,8,15) #from Zhou et al 2016... NF modified to account for maritime -- 37% higher than temperate (so, scaled up by 0.5*.37=0.185)
      pTable[41, 2:8]<-20
      write.csv(pTable, paste("AssimSites/HUC_", huc, "/pTable.csv", sep=""), row.names=F)
      
      #read in PFT split from MsTMIP grid
      lulc_1<-read.csv("DataProcessing/MsTMIP/lulcc_allDA.csv", stringsAsFactors = F)
      lulc_1$HUC8<-fix8(lulc_1$HUC8)
      lulc_1<-lulc_1[lulc_1$HUC8==huc, ]
      lulc_1$Yr<-as.numeric(lulc_1$yr)
      lulc_2<-lulc_1[which(lulc_1$yr == 2009), ]
      lulc_2<-lulc_2[which(lulc_2$PercentCov > 0),]
      
      lulc_2$cov_each_pft<-NA
      for(i in 1:nrow(lulc_2)){
        if(sum(lulc_2[i,6:12]) > 0){
          frac_each<-1/sum(lulc_2[i,6:12]) #get the fraction of the param cov for each pft
          lulc_2$cov_each_pft[i]<-(lulc_2$PercentCov[i]*frac_each)*lulc_2$cell_cov_area[i] #calculate the area of the cell within the huc covered by each pft
        }
      }
      
      #calculate area of each PFT in huc
      PFTs<-names(which(colSums(lulc_2[,6:12]) > 0))
      pft_cov<-data.frame(matrix(ncol=3, nrow=length(PFTs)))
      colnames(pft_cov)<-c("year", "pft", "cov_m2")
      pft_cov$year<-2009
      for(i in 1:length(PFTs)){
        pft_cov$cov_m2[i]<-sum(lulc_2[which(lulc_2$PercentCov > 0), which(names(lulc_2)==PFTs[i])]*lulc_2$cov_each_pft[which(lulc_2$PercentCov > 0)], na.rm=T)
        pft_cov$pft[i]<-PFTs[i]
      }
      pft_cov$pct<-pft_cov$cov_m2/sum(pft_cov$cov_m2)*100
      pft_cov<-pft_cov[which(pft_cov$pct > 5), ]
      pftList<-pft_cov$pft
      
      
      write.csv(pft_cov, paste("AssimSites/HUC_", huc, "/pft_cov.csv", sep=""), row.names=F)
      ##check: these should be roughly the same, but pft_cov prob slightly lower
      #sum(pft_cov$cov_m2)
      #sum(unique(lulc_2$cell_cov_area))
      
      #read in gpp data
      GPP_8day<-read.csv(paste("DataProcessing/MODIS/", "MODIS_GPP_HUC_", huc,".csv", sep=""))#pull in LAI
      GPP_8day$year<-substr(GPP_8day$date, 1, 4)
      
      #get the start date by seeing when we start to have LAI-- most limiting timeseries
      minStart<-"2003-01-01" #lai_4day$date[1]
      maxEnd<-"2021-12-31"
      
      #subset GPP
      GPP_8day<-GPP_8day[which(GPP_8day$date >= minStart), ]
      
      #create environmental forcings
      ppt<-read.csv("DataProcessing/PRISM/PRISM_pptallHUC.csv") #pull in PRISM data
        ppt$huc<-fix8(ppt$huc)
        ppt<-ppt[which(ppt$huc==huc), ]
        ppt$value<-ppt$value/10 #mm to cm
      vpdmin<-read.csv("DataProcessing/PRISM/PRISM_vpdminallHUC.csv") #pull in PRISM data
        vpdmin$huc<-fix8(vpdmin$huc)
        vpdmin<-vpdmin[which(vpdmin$huc==huc), ]
      vpdmax<-read.csv("DataProcessing/PRISM/PRISM_vpdmaxallHUC.csv") #pull in PRISM data
        vpdmax$huc<-fix8(vpdmax$huc)
        vpdmax<-vpdmax[which(vpdmax$huc==huc), ]
      tmean<-read.csv("DataProcessing/PRISM/PRISM_tmeanallHUC.csv") #pull in PRISM data
        tmean$huc<-fix8(tmean$huc)
        tmean<-tmean[which(tmean$huc==huc), ]
      forc_data<-cbind(ppt, vpdmin$value , vpdmax$value , tmean$value)
      colnames(forc_data)[3:6]<-c("P", "vpdmin", "vpdmax", "tmean")
      forc_data$Date<-NA
      
      for(r in 1:nrow(forc_data)){
        datStr<-strsplit(forc_data$file[r], split="_" )[[1]][5]
        forc_data$Date[r]<-paste(substr(datStr, 1,4), substr(datStr, 5,6), substr(datStr, 7,8), sep="-")
      }
      
      forc_data<-forc_data[which(forc_data$Date >= minStart), ]
      forc_data$DOY<-yday(forc_data$Date)#add DOY column
      
      outCoords<-getMidLat(huc=huc, wbdDir="/Users/cearatalbot/HUC04Example/myWBD/") #get the lat of middle of watershed
      
      forc_data<-calcPAR(df=forc_data, lat=(as.numeric(outCoords$x)*-1)) #calculate PAR from ET rad using library(sirad)
      forc_data$vpdmean<-((forc_data$vpdmin+forc_data$vpdmax)/2)/10 #hPa to kPa
      #aqEvap() #calculate aquatic evaporation
      forc_data$SnowMelt_cm<-0
      
      GPP_8day$diffs<-NA #calculate the difference between days
      GPP_8day$diffs[1]<-8 #if the df is already subsetted so that the first days include the full 8 days
      for(i in 2:nrow(GPP_8day)){
        GPP_8day$diffs[i]<-round(as.numeric(difftime(GPP_8day$date[i], GPP_8day$date[(i-1)], units="days")), digits=0)
      }
      
      write.csv(GPP_8day, paste("AssimSites/HUC_", huc, "/MODISData/gpp_data.csv", sep=""), row.names=F)
      write.csv(forc_data, paste("AssimSites/HUC_", huc, "/ForcingData/clim.csv", sep=""), row.names=F)

      
      #setup initial conditions
      ##run equilibrium to get initial conditions
      initDF<-data.frame(matrix(ncol=16, nrow=length(pftList))) #initial conditions setup
      colnames(initDF)<-c("Cw","Cl", "Cs1", "Cs2", "Cs3", "Cs4", "Cdoc1", "Cdoc2", "W1", "W2", "Ca", "Cr", "Ccwd", "Cdic1", "Cdic2", "PFT")
      initDF$PFT<-pftList
      
      for(i in 1:nrow(initDF)){
        params<-pTable[,which(colnames(pTable)==initDF$PFT[i])]  #parms for each pft
        names(params)<-pTable$pName #name params
        if(initDF$PFT[i]=="EGBR"){
          initDF[i,1:(ncol(initDF)-1)]<-c(13000, 30, 60, 60, 300, 1000, 50, 30, as.numeric(params[[19]]*.1), as.numeric(params[[41]]*.1), 600000, 3000, 50, 0.01, 0.01) #fill initial conditions
        } else if(initDF$PFT[i]=="EGNE"){
          initDF[i,1:(ncol(initDF)-1)]<-c(13000, 30, 60,60, 300, 1000, 50, 30, as.numeric(params[[19]]*.1), as.numeric(params[[41]]*.1), 600000, 3000, 50, 0.01, 0.01) #fill initial conditions
        } else if(initDF$PFT[i]=="DEBR"){
          initDF[i,1:(ncol(initDF)-1)]<-c(13000, 0, 60,60, 300, 1000, 50, 30, as.numeric(params[[19]]*.1), as.numeric(params[[41]]*.1), 600000, 3000, 50, 0.01, 0.01) #fill initial conditions
        } else if(initDF$PFT[i]=="DENE"){
          initDF[i,1:(ncol(initDF)-1)]<-c(13000, 0, 60,60, 300, 1000, 50, 30, as.numeric(params[[19]]*.1), as.numeric(params[[41]]*.1), 600000, 3000, 50, 0.01, 0.01) #fill initial conditions
        } else if(initDF$PFT[i]=="SH"){
          initDF[i,1:(ncol(initDF)-1)]<-c(10000, 0, 60,60, 300, 1000, 50, 30, as.numeric(params[[19]]*.1), as.numeric(params[[41]]*.1), 600000, 2500, 50, 0.01, 0.01) #fill initial conditions
        } else if(initDF$PFT[i]=="GR" | initDF$PFT[i]=="CR"){
          initDF[i,1:(ncol(initDF)-1)]<-c(0, 0, 60,60, 300, 1000, 50, 30, as.numeric(params[[19]]*.1), as.numeric(params[[41]]*.1), 600000, 500, 50, 0.01, 0.01) #fill initial conditions
        }
      }
      
      write.csv(initDF, paste("AssimSites/HUC_", huc, "/inits.csv", sep=""), row.names=F)
  
}