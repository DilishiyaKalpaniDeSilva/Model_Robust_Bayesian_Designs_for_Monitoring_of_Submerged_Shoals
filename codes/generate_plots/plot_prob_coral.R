## ---------------------------------------------------------------------------------------------
## Authors: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  1. Predict mean probability of hard coral for whole shoal
##           2. Plot optimal GAMM designs with mean probability of hard coral
## ---------------------------------------------------------------------------------------------

##################################### CODE #####################################################

library(DSpat)
library(raster)
library(sf)
library(rgdal)
library(RColorBrewer)
library(spatstat)
library(ggplot2)
library(ggpubr)
library(rasterVis)

source("codes/functions/cal_z.R")
source("codes/functions/basic_functions.R")

load(paste0("data/transect_data/data_for_design_fshsize_predang100_n50_dep50_len500.RData"))
load("data/BE_raster_data/BE_raster.RData")

cbbPalette <- viridis::viridis(6)
rp <- rasterToPoints(covP3)
data_s <- data.frame(rp)

## 1. Predict mean probability of hard coral for whole shoal
year_all = c(2010,2011,2013, "yre")
data <- NULL
for(indexY in 1:4){
  
  year = year_all[indexY]
  load(paste0("outputs/priors/GAMM/GAMM_prior_es_FishNet50_",year,".RData")) ## loading priors
  bu <- c(out_Lap[[1]]$par[1:2],out_Lap[[2]]$par[1:25])
  knot_loc <- out_Lap[[3]]
  num.knots <- k <- length(knot_loc)
  intKnots <- knot_loc[-c(1, length(knot_loc))]
  
  print(year)
  num.knots <- 25
  n <- nrow(data)
  x <- as.matrix(norm_01(data_s$depth))
  zz <- cal_z(numBasis = num.knots,XsplPreds=x,knot_type="es" ,intKnots,
              boundKnots=range(x), basis ="OS")
  z.spline <- zz$Z
  C <- cbind(1,x,z.spline)
  
  data_s$prob <- expit(C%*%bu)
  if(indexY==4){
    data_s$Prior <- "yre"
  }else{
    data_s$Prior <- year    
  }
  
  data <- rbind(data,data_s)
}

## 2. Plot optimal GAMM designs with mean probability of hard coral

year <- c("2010","2011","2013","yre")
nt <- length(lnsp[[1]]$transects)

data_t <- NULL
for(j in 1:4){
  load(paste0("outputs/designs/GAMM/",year[j],"/Optimal.RData"))
  tr <- sort(tr_design)
  print(tr)
  
  for(i in 1:18){
    point <- ceiling(tr[i]/nt)
    print(i)
    t_no <- tr[i]-((point-1)*nt)
    data_t1 <- data.frame((lnsp[[point]]$transects)[[t_no]],Tr=paste0("Tr",i), Prior=year[j])
    data_t <- rbind(data_t, data_t1)
  }
}

cbbPalette <- viridis::viridis(6)
rp <- rasterToPoints(covP3)
data_s <- data.frame(rp)

ggplot() +
  geom_tile(data=data,aes(x=x,y=y,fill = prob)) +
  geom_polygon(data = data_t,aes(x=x,y=y,group=Tr)) +
  xlab("Easting") + ylab("Northing") +
  scale_fill_gradient(na.value="white", name="Probability",
                      low = "yellow",
                      high = "red") +
  scale_x_continuous(breaks = (xd), labels=(xd),limits=range(xd)+c(-500,250))+
  scale_y_continuous(breaks = yd, labels=yd,limits=range(yd)+c(-250,250))+
  coord_equal() + theme_minimal() +
  geom_vline(xintercept = xd)+
  geom_hline(yintercept = yd) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle = 270),
        panel.background = element_blank(),
        strip.background = element_rect(fill="lightgrey",color="black"),
        strip.text = element_text(colour = "black",size = 12,face="bold"),
        title = element_text(size=12,face="bold"),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(color = "black",size=14,face="bold"),
        axis.text.x = element_text(color = "black", size=8,angle = 90),
        axis.text.y = element_text(color = "black", size=8,angle = 0),
        legend.text = element_text(color = "black", size=12),
        legend.position = "right", 
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.margin=margin(0,0,0,0),
        plot.margin=grid::unit(c(5,-5,3,-5), "mm")) +
  facet_wrap(~ Prior, ncol = 2)

ggsave(paste0("plots/designs/Ex2_optimal_designs_prob.jpeg"),width = 3*3,height=2*3, bg="white")

