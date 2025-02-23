## ---------------------------------------------------------------------------------------------
## Authors: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  1. Save BE raster for depth range -50 to -19.14501
##           2. Generate fishnet polygon data to use in find_fishnet function
## ---------------------------------------------------------------------------------------------

################################ CODE ###########################################################

library(NCmisc)
## check the packages used in the R script
list.functions.in.file("codes/generate_data/generate_fishnet_data.R", alphabetic = TRUE) 

library(raster)
library(sf)

pdf(file = paste0("plots/generate_data/plot_generate_fishnet_data.pdf"))

# Download Barracouta East raster data from https://doi.org/10.25845/p1sr-s997.
covP <- raster("data/BE_raster_data/depth.img")
out <- match(covP,table = -9999)
covP1 <- covP                          # assigning NA's for outer values
covP1[out] <- NA
rp <- rasterToPoints(covP1)            # converting the depth values of shoal values to a data frame
max_depth <- -50
# min_depth <-  -19.14501  
rb <- subset(rp,rp[,3] < max_depth)    # depth greater than 30
mac <- (unique(rb[,3]))
out1 <- match(covP1,table = mac)       # removing them from the raster
covP2<- covP1
covP2[out1] <- NA
lin <- rasterToContour(is.na(covP2))
rp <- rasterToPoints(covP2)
range(rp[,3])
# rb <- subset(rp,rp[,3]> min_depth)     # depth less than 18.202
# mac <- (unique(rb[,3]))
# out1 <- match(covP2,table = mac)       # removing them from the raster
covP3<- covP2
covP3[out1] <- NA
rp <- rasterToPoints(covP3)
range(rp[,3])

## a grid that covers the shoal
re <-  c( 610812.6, 613862.6) 
rn <-  c(8611923, 8613923) 
xd <- seq(from=re[1],to =re[2],by= 500)
yd <- seq(from=rn[1],to =rn[2],by=500) 

raster::plot(covP3,  xlim=c(610000,614200), xlab="Easting", ylab="Northing",legend.only=FALSE,
     legend.args=list(text='Depth (m)',side=4,font=2, line=2.5, cex=0.8))
abline(v=xd)
abline(h=yd)

## starting points of the transects
inter_points <- NULL              
for(i in 1:length(xd)){
  for (j in 1:length(yd)) {
    inter_points <- rbind(inter_points, c(x0=xd[i],y0=yd[j]))
  }
}

## choose starting points that are inside the shoal
dep <- raster::extract(covP, inter_points)
print(dep)
pp <- which(dep<(max_depth))
inter_points <- inter_points[-pp,]
text(inter_points)
dep <- raster::extract(covP, inter_points)
print(dep)

save(list = c("covP1","covP3","xd","yd","inter_points"),file = "data/BE_raster_data/BE_raster.RData") 

## fishnet raster as fishnet polygon data
fishnet_BE <- st_read("data/BE_raster_data/BarracoutaEast_Fishnet.shp")
dim(fishnet_BE)
f <- (fishnet_BE$geometry)
plot(f,add=TRUE)

par(mfrow=c(3,2))
size <- seq(50,250,50)
fishnet_data <- list()

for(i in 1:length(size)){
  e <- as(raster::extent(610812.6-500, 613862.6+500, 8611923-250, 8613923+250), "SpatialPolygons") %>% 
    st_as_sf()
  f <- st_make_grid(e, cellsize = c(size[i], size[i]))
  plot(covP3,  xlim=c(610000,614200), xlab="Easting", ylab="Northing",legend.only=FALSE,
       legend.args=list(text='Depth (m)',side=4,font=2, line=2.5, cex=0.8))
  plot(f, add=TRUE)
  abline(v=xd,col="red")
  abline(h=yd,col="red")
  
  co <- raster::geom(as(f, "Spatial"))
  dim(co)
  length(unique(co[,1]))
  fishnet_data[[i]] <- data.frame(fishnet=co[,1],co[,5:6])
  print(length(unique(fishnet_data[[i]]$fishnet)))
  
}

save(list = c("fishnet_data"),file = "data/BE_raster_data/fishnets_diff_sizes.RData")

dev.off()
