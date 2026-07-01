
# generic initial conditions based on pft

gen_initials<-function(pftList, PFTtable){
  
        initDF<-data.frame(matrix(ncol=19, nrow=length(pftList))) #initial conditions setup
        colnames(initDF)<-c("Cw","Cl", "Cs1", "Cs2", "Cs3", "Cs4", "Cdoc1", "Cdoc2", "W1", "W2", "Ca", "Cr", "Ccwd", "Cdic1", "Cdic2", "W3", "Cdic3", "Cdoc3", "PFT")
        initDF$PFT<-pftList
        
        for(i in 1:nrow(initDF)){
          params<-as.numeric(PFTtable[,which(colnames(PFTtable)==initDF$PFT[i])])  #parms for each pft
          names(params)<-PFTtable$pName #name params
          if(initDF$PFT[i]=="EGBR"){
            initDF[i,1:(ncol(initDF)-1)]<-c(13000, as.numeric(params[[55]]*params[[20]]*params[[21]]), 60, 60, 300, 2000, 50, 30, as.numeric(params[[19]]*.90), as.numeric(params[[41]]*.90), 600000, 3000, 50, 0.01, 0.01, as.numeric(params[[41]]*.90), 0.0001, 10) #fill initial conditions
          } else if(initDF$PFT[i]=="EGNE"){
            initDF[i,1:(ncol(initDF)-1)]<-c(13000, as.numeric(params[[55]]*params[[20]]*params[[21]]), 60,60, 300, 2000, 50, 30, as.numeric(params[[19]]*.90), as.numeric(params[[41]]*.90), 600000, 3000, 50, 0.01, 0.01,  as.numeric(params[[41]]*.90), 0.0001, 10) #fill initial conditions
          } else if(initDF$PFT[i]=="DEBR"){
            initDF[i,1:(ncol(initDF)-1)]<-c(13000, as.numeric(params[[55]]*params[[20]]*params[[21]]), 60,60, 300, 2000, 50, 30, as.numeric(params[[19]]*.90), as.numeric(params[[41]]*.90), 600000, 3000, 50, 0.01, 0.01,  as.numeric(params[[41]]*.90), 0.01, 10) #fill initial conditions
          } else if(initDF$PFT[i]=="DENE"){
            initDF[i,1:(ncol(initDF)-1)]<-c(13000, as.numeric(params[[55]]*params[[20]]*params[[21]]), 60,60, 300, 2000, 50, 30, as.numeric(params[[19]]*.90), as.numeric(params[[41]]*.90), 600000, 3000, 50, 0.01, 0.01,  as.numeric(params[[41]]*.90), 0.01, 10) #fill initial conditions
          } else if(initDF$PFT[i]=="SH"){
            initDF[i,1:(ncol(initDF)-1)]<-c(10000, as.numeric(params[[55]]*params[[20]]*params[[21]]), 60,60, 300, 2000, 50, 30, as.numeric(params[[19]]*.90), as.numeric(params[[41]]*.90), 600000, 2500, 50, 0.01, 0.01,  as.numeric(params[[41]]*.90), 0.01, 10) #fill initial conditions
          } else if(initDF$PFT[i]=="GR" | initDF$PFT[i]=="CR"){
            initDF[i,1:(ncol(initDF)-1)]<-c(0, as.numeric(params[[55]]*params[[20]]*params[[21]]), 60,60, 300, 2000, 50, 30, as.numeric(params[[19]]*.90), as.numeric(params[[41]]*.90), 600000, 500, 50, 0.01, 0.01,  as.numeric(params[[41]]*.90), 0.01, 10) #fill initial conditions
          }
        }
     
return(initDF)       
}
