#add time for each PFT. ave applies seq_along to each pft group
library(ggplot2)
library(tidyr)
library(dplyr)
spinupALL$time <- unlist(tapply(seq_len(nrow(spinupALL)), spinupALL$pft, seq_along))#deconstruct df along time


plot_spinup <- function(data) {
  ggplot(data, aes(x = time, y = Ca, color = pft)) +
    geom_line() +
    labs(x = "time", y = "Ca", color = "PFT")
}

plot_spinup(spinupALL)

plot_data <- spinupALL %>%
  filter(pft=="EGNE")
plot_data <- plot_data[-(1:5),]
  

ggplot(plot_data, aes(x=time))+
         geom_line(aes(y=Ci), color = "lightblue4")+
         geom_line(aes(y=Ca), color = "blue")+
         geom_line(aes(y=Alg),color= 'olivedrab3')+
         labs(x = "time", y = "Carbon")

ggplot(plot_data,aes(x=time,y=Ca/Ci))+geom_line()
ggplot(plot_data,aes(x=time,y=Alg))+geom_line()

spinupALL[which(spinupALL$pft=='EGNE'),]