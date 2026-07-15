#add time for each PFT. ave applies seq_along to each pft group
library(ggplot2)
spinupAll$time <- unlist(tapply(seq_len(nrow(spinupAll)), spinupAll$pft, seq_along))


plot_spinup <- function(data) {
  ggplot(data, aes(x = time, y = Ca, color = pft)) +
    geom_line() +
    labs(x = "time", y = "Ca", color = "PFT")
}

plot_spinup(spinupAll)

plot_spinup1 <- function(data) {
  ggplot(data, aes(x = time, y = Calat, color = pft)) +
    geom_line() +
    labs(x = "time", y = "Calat", color = "PFT")
}


#plot_spinup1(spinupAll)