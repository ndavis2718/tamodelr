#####Growing degree days function
#input is temps vector for a single year
GddCalc<-function(temps, ths, chill=FALSE){
                outTemps<-c(ifelse(temps[1] > 0, temps[1], 0), rep(0, (length(temps)-1))) #setup the initial temps
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
