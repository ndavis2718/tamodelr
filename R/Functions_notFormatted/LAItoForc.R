
LAItoForc<-function(lai_df, pctAbv, SLW_param, Cfrac_param){

  thrs<-min(lai_df$LAI[which(lai_df$LAI >= 0)])+min(lai_df$LAI[which(lai_df$LAI >= 0)])*.50
  lai_df$LAI[which(lai_df$LAI==-999)]<-0

  
  #lai_df$LAI[which(lai_df$LAI < thrs)]<-0
  #these were hard to figure out. L and Ll are the flux of leaves on and off. So they are 0 unless leaves are actively coming in (L) or actively falling (Ll)          lai_df$L<-0 #ifelse(lai_df$LAI >= thrs, lai_df$LAI*params[[20]]*params[[21]], 0)             
  lai_df$L<-0
  for(y in 2:nrow(lai_df)){
    #if(lai_df$LAI[y-1] >= thrs){
      lai_df$L[y]<-ifelse((lai_df$LAI[y]-lai_df$LAI[y-1])>= 0, ((lai_df$LAI[y]-lai_df$LAI[y-1])*SLW_param*Cfrac_param), 0)
    #} 
  }
  lai_df$L<-round(lai_df$L, digits=2)
  lai_df$Ll<-0
  for(y in 2:nrow(lai_df)){
   # if(lai_df$LAI[y-1] >= thrs){
      lai_df$Ll[y]<-ifelse((lai_df$LAI[y-1]-lai_df$LAI[y])>= 0, ((lai_df$LAI[y-1]-lai_df$LAI[y])*SLW_param*Cfrac_param), 0)
   # }
  }
  lai_df$Ll<-round(lai_df$Ll, digits=2)
  for(y in 1:length(unique(lai_df$Year))){
    phenolYr<-lai_df[which(lai_df$Year==unique(lai_df$Year)[y]), ]
    if((sum(phenolYr$Ll)-sum(phenolYr$L)) > 0 | (sum(phenolYr$Ll)-sum(phenolYr$L))< 0)
      dif=sum(phenolYr$L) - sum(phenolYr$Ll)
    if(dif < 0){
      phenolYr$L<-phenolYr$L+(dif/nrow(phenolYr)*-1)
    } else{
      phenolYr$Ll[which(phenolYr$Ll > 0)][length(phenolYr$Ll[which(phenolYr$Ll > 0)])]<-(phenolYr$Ll[which(phenolYr$Ll > 0)][length(phenolYr$Ll[which(phenolYr$Ll > 0)])])+dif#last val of Ll
    }
    lai_df$Ll[which(lai_df$Year==unique(lai_df$Year)[y])]<-phenolYr$Ll
    lai_df$L[which(lai_df$Year==unique(lai_df$Year)[y])]<-phenolYr$L
  }
  

return(lai_df)
} #end fun
