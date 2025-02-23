## ---------------------------------------------------------------------------------------------
## Authors: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  Generate the following required for the design algorithm:
#              i) Transects  
#             ii) Coordinates for data collection locations inside the 
#                 generated transects (Easting, Northing)
#             iii) Covaraite values and fishnet number at these locations
## ---------------------------------------------------------------------------------------------

################################ CODE ###########################################################

library(raster)
library(sf)
library(RColorBrewer)
library(spatstat)
library(DSpat)

source("codes/functions/find_fishnet.R")
source("codes/functions/basic_functions.R")

load("data/BE_raster_data/BE_raster.RData") 

pdf(file = paste0("plots/generate_data/plot_generate_transect_data.pdf"),   
    width = 8, # The width of the plot in inches
    height = 6,
)

indexL <- 1              ## index of transect length
nt <-100                 ## number of transects around one point
n <- 50*indexL           ## number of data per transect
max_depth <- -50         ## maximum depth for the design
min_depth <- -19.14501         ## maximum depth for the design
size <- seq(50,250,50)   ## different fishnet sizes
  
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- rep(col_vector,50)

length <- c(500,1000,1500)
length <- length[indexL]

study.area=owin(xrange=c(range(xd)[1]-500,range(xd)[2]+500),
                yrange=c(range(yd)[1]-100,range(yd)[2]+100))

## generate 100 transects around each of 22 starting points, for ACE algorithm 
plot(covP3,  xlim=c(610000,614200), xlab="Easting", ylab="Northing",legend.only=FALSE,      
     legend.args=list(text='Depth (m)',side=4,font=2, line=2.5, cex=0.8))
abline(v=xd)
abline(h=yd)

zp <- seq(0,2,by=2/nt)[-(nt+1)] ## 100 different angles to predict from Kringing model
lnsp <- list()
for(j in 1:nrow(inter_points)){
  x1 <- vector()
  y1 <- vector()
  for(i in 1: length(zp)){
    x1[i] <- inter_points[j,1] + length*cospi(zp[i])
    y1[i] <- inter_points[j,2] + length*sinpi(zp[i])
  }
  
  lnsp[[j]] <- lines_to_strips(lines = data.frame(x0=inter_points[j,1],
                                                  x1=x1,
                                                  y0=inter_points[j,2],
                                                  y1=y1),
                               study.area=study.area
                               ,width=50)
  print(j)
  plot(owin(poly=lnsp[[j]]$transects),col=col_vector[j],add=TRUE)
}
text(inter_points)

z <- seq(0,2,by=2/20)[-(20+1)]  ## 20 different angles to fit the Kringing model    
indz <- match(round(z,2),round(zp,2))
plot(covP3,  xlim=c(610000,614200), xlab="Easting", ylab="Northing",legend.only=FALSE,
     legend.args=list(text='Depth (m)',side=4,font=2, line=2.5, cex=0.8))
abline(v=xd)
abline(h=yd)

study.area=owin(xrange=c(range(xd)[1]-500,range(xd)[2]+500),
                yrange=c(range(yd)[1]-100,range(yd)[2]+100))
for(j in 1:nrow(inter_points)){
  print(j)
  plot(owin(poly=lnsp[[j]]$transects[indz]),col=col_vector[j],add=TRUE)
}
text(inter_points)

## coordinates inside generated transects
set.seed(2)
coordp <- list()
index <- 1
for(i in 1:length(lnsp)){
  for(j in 1:length(zp)){
    coordp[[index]] <- coords(runifpoint(n, win=owin(poly=lnsp[[i]]$transects[j]),giveup = 10000000))
    index <- index +1
  }
}

plot(covP3,  xlim=c(610000,614200), xlab="Easting", ylab="Northing",legend.only=FALSE,      
     legend.args=list(text='Depth (m)',side=4,font=2, line=2.5, cex=0.8))
abline(v=xd)
abline(h=yd)

for(i in 1:length(coordp)){
  points(coordp[[i]])
}

depthp <- list();ll <- vector();j <- 1
for(i in 1:length(coordp)){
  print(i)
  depthp[[i]] <- raster::extract(covP1, coordp[[i]])
  rm <- which((depthp[[i]] < max_depth) | (depthp[[i]] > min_depth))
  print(rm)
  
  if(length(rm)>=1){
    coordp[[i]] <- NA
    depthp[[i]] <- NA
  }else{
    ll[j] <- (length(depthp[[i]]))
    j <- j+1
  }
}
all(ll==50)

plot(covP3,  xlim=c(610000,614200), xlab="Easting", ylab="Northing",legend.only=FALSE,      
     legend.args=list(text='Depth (m)',side=4,font=2, line=2.5, cex=0.8))
abline(v=xd)
abline(h=yd)

for(i in 1:length(coordp)){
  points(coordp[[i]])
}
text(inter_points,col="red")

## find fishnet number for each location
load("data/BE_raster_data/fishnets_diff_sizes.RData")
fishnet_all <- list()

for(ind_fish_size in 1:5){
  fishnetp <- list()
  for (i in 1:length(coordp)) {
    if(length(coordp[[i]])==1){
      fishnetp[[i]] <- NA
      print(i)
    }else{
      fishnetp[[i]] <- find_fishnet(fishnet_data[[ind_fish_size]],coordp[[i]])
      print("***")
      print(i)
    }
  }
  fishnet_all[[ind_fish_size]] <- fishnetp
}

fll <-list()
fl <- vector()
for (ind_fish_size in 1:5) {
  for(i in 1:length(fishnetp)){
    fl[i] <- (length(unique(fishnet_all[[ind_fish_size]][[i]])))
  } 
  fll[[ind_fish_size]] <- fl 
  print(table(fl))
}

save(list = c("inter_points","lnsp","coordp","depthp","fishnet_all"),
     file = paste0("data/transect_data/data_for_design_fshsize_predang",nt,"_n",n,"_dep",
                   (-max_depth),"_len",length,".RData"))
dev.off()
