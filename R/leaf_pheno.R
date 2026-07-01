###Functions to calculate leaf phenology
#' @title leaf_pheno
#'
#' @description Calculates leaf phenology forcings (leaf on and leaf off)
#'
#' @param temps Vector of mean daily temperature (deg C)
#' @param ths Growing degree day threshold (deg C)
#' @param chill Calculate chilling days (set to TRUE)
#'
#' @return A data frame object that contains calculated growing degree days
#'              \code{summarize_tam}
#' @examples
#' data(sing_watershed)
#' output_table <- overview_tab(t = toydata, S = ccode, p = year)
#' @export


#Growing degree day calculations
#input is temps vector for a single year, ths is the temperature threshold
gddCalc<-function(temps, ths, chill=FALSE){

  outTemps<-c(ifelse(temps[1] > ths, temps[1], 0), rep(0, (length(temps)-1))) #setup the initial temps

  if(chill==TRUE){ # to count chilling days
    for(i in 1:length(temps)){
      if(temps[i] < ths & temps[i] > 0){
        outTemps[i]<-1 #yes, a chilling day
      } else{
        outTemps[i]<-0 #no, not a chilling day
      }
    }
  }else{ #growing degree days
    for(i in 2:length(temps)){
      if(temps[i] > ths){
        outTemps[i]<-outTemps[(i-1)]+temps[i]
      } else{
        outTemps[i]<-outTemps[(i-1)]
      }
    }
  }

  return(outTemps)
}

# Calculate leaf on and leaf off date, formatted for forcing input into proc_model()
####Leaf module setup, based on CABLE phenology from Haverd et al. 2018; https://doi.org/10.5194/gmd-11-2995-2018
# required arguments include data frame with "tair" and "year" columns and pft-specific phenology parameters a, b, k
addLeaves<-function(df, a, k, b, glength){
  df$GDD<-0 #g rowing degree days
  df$DOY<-0 # day of year
  df$GDDc<-0
  df$GDD0<-0
  df$grow<-0
  df$greenup<-0
  df$sen<-0 #10-day cumulative degree days (for leaf off)
    for(i in 1:length(unique(df$year))){
      sub<-df[df$year==unique(df$year)[i], ]
      df$GDD[df$year==unique(df$year)[i]]<-gddCalc(temps=sub$tair, ths=0, chill=FALSE)
      df$DOY[df$year==unique(df$year)[i]]<-1:nrow(sub)
      df$GDDc[df$year==unique(df$year)[i]]<-gddCalc(temps=sub$tair, ths=0, chill=TRUE)
      ##GDD0 following Sykes et al 1996: GDD0=a+be-kC; C = length of chilling period, ab&k are species specific
      df$GDD0[df$year==unique(df$year)[i]]<-a+b*exp(-k*(length(which(df$GDDc[df$year==unique(df$year)[i]] == 1 & df$DOY[df$year==unique(df$year)[i]] < glength)))) #create relationship between this and chilling days
      df$greenup[df$year==unique(df$year)[i]]<-ifelse(df$GDD[df$year==unique(df$year)[i]]-df$GDD0[df$year==unique(df$year)[i]] < glength & df$GDD[df$year==unique(df$year)[i]] > df$GDD0[df$year==unique(df$year)[i]], 1, 0) #when gdd > gdd0 leaf on starts, adds 1s until GDD-GDD0 > 200
      df$sen[df$year==unique(df$year)[i]]<-c(rep(0, (which(df$greenup[df$year==unique(df$year)[i]] > 0)[1]+glength)), rep(1, nrow(df[df$year==unique(df$year)[i], ])-(which(df$greenup[df$year==unique(df$year)[i]] > 0)[1]+glength))) #200 days after the onset of leaf growth, leaves begin to senesce
      df$grow[df$year==unique(df$year)[i]]<-ifelse(df$GDD[df$year==unique(df$year)[i]] > df$GDD0[df$year==unique(df$year)[i]] & df$greenup[df$year==unique(df$year)[i]] == 0 & df$sen[df$year==unique(df$year)[i]] == 0, 1, 0) #when gdd > gdd0 leaf on starts, adds 1s until GDD-GDD0 > 200
    }

  return(df)
}
