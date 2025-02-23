library(raster)
library(spatstat)
library(RColorBrewer)

plot_design <- function(tr_design){
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- rep(col_vector,50)
  
  load("data/BE_raster_data/BE_raster.RData")
  
  plot(covP3,  xlim=c(610400,614200), xlab="Easting", ylab="Northing",legend.only=FALSE,
       legend.args=list(text='Depth (m)',side=4,font=2, line=2.5, cex=0.8))
  abline(v=xd)
  abline(h=yd)
  text(inter_points)
  
  nt <- length(lnsp[[1]]$transects)
  
  load("data/transect_data/data_for_design_fshsize_predang100_n50_dep50_len500.RData")
  trd1 <- sort(tr_design)
  print(trd1)
  point <- ceiling(trd1/nt)
  print(point)
  rem <- trd1 - (point-1)*nt
    
  for(i in 1:length(tr_design)){
    plot(owin(poly=lnsp[[point[i]]]$transects[[rem[i]]]),col=col_vector[i],add=TRUE)
  }
  text(inter_points)
  
}


plot_design_AIMS <- function(tr_design,dy){
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- rep(col_vector,50)
  
  load("data/BE_raster_data/BE_raster.RData")
  
  plot(covP3,  xlim=c(610400,614200), xlab="Easting", ylab="Northing",legend.only=FALSE,
       legend.args=list(text='Depth (m)',side=4,font=2, line=2.5, cex=0.8))
  abline(v=xd)
  abline(h=yd)
  
  if(dy == 2010){
    load(paste0("data/transect_data/data_for_design_2010.RData"))
    for(i in tr_design){
      plot(owin(poly=lnsp[[i]]$transects),col=col_vector[i],add=TRUE)
    }
  }else if(dy==2011){
    load(paste0("data/transect_data/data_for_design_2011_13.RData"))
    for(i in tr_design){
      plot(owin(poly=lnsp[[i]]$transects),col=col_vector[i],add=TRUE)
    }
  }else{
    
  }
}



