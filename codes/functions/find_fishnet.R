library(sp)
library(magrittr)

find_fishnet <- function(fishnet_data,coord){
  #load("data/BE_raster_data/fishnet_data.RData")
  fishnet_data %>% 
    split(.$fishnet) %>%
    sapply(function(pol) point.in.polygon(coord[,1], coord[,2], pol$x, pol$y) > 0) %>%
    apply(1,which)
}

