## ---------------------------------------------------------------------------------------------
## Authors: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  Plot depth vs predicted for GAMM models
## ---------------------------------------------------------------------------------------------

################################# CODE #########################################################

library(NCmisc)
list.functions.in.file("codes/generate_plots/plot_depth_vs_pred_GAMM.R", alphabetic = TRUE)

library(plyr)
library(ggplot2)
library(readxl)
library(tidyr)
library(ggpubr)
library(grid)

source("codes/functions/cal_z.R")
source("codes/functions/basic_functions.R")

BE_data <- read_excel("data/BE_data/Barracouta_East_Rdata_FS.xlsx")
year <- c(2010,2011,2013,"yre")

sum_data <- function(pname, yvec, kp){
  data_all <- NULL
  for(indexY in yvec){
    if(year[indexY]=="yre"){
      data <- BE_data
    }else{
      data <- subset(BE_data,BE_data$Year==year[indexY])
    }
    range(data$depth)
    unique(data$Year)

    suppressMessages(attach(data))
    dim(data)
    print(range(data$depth))
    print(unique(data$Year))
    
    load(paste0("outputs/priors/GAMM/GAMM_prior_",pname,"FishNet50_",year[indexY],".RData")) ## loading priors
    knot_loc <- out_Lap[[3]]
    num.knots <- k <- length(knot_loc)
    intKnots <- knot_loc[-c(1, length(knot_loc))]
    
    # predictions for GAMM model
    n <- nrow(data)
    y <- data$HardCoral
    x <- as.matrix(norm_01(data$depth))
    zz <- cal_z(numBasis = num.knots,XsplPreds=x,knot_type=kp ,intKnots,
                boundKnots=range(x), basis ="OS")
    z.spline <- zz$Z
    C <- cbind(1,x,z.spline)
    m <- length(unique(data$FishNet50))

    
    print(out_Lap[[1]]$convergence)
    print(out_Lap[[2]]$convergence)
    bu <- c(out_Lap[[1]]$par[1:2],out_Lap[[2]]$par[1:25])
    U <- out_Lap[[2]]$par[26:(25+m)]
    random <- factor(data$FishNet50)
    random <- mapvalues(random, from =levels(random), to=c(1:m))
    Year <- mapvalues(data$Year, from =levels(as.factor(data$Year)), to=c(1:length(unique(data$Year))))

    if(year[indexY]=="yre"){
      Y <- tail(out_Lap[[2]]$par,length(unique(Year))-1)
      Y1 <- c(0,Y)

      datapp <- data.frame(Depth= data$depth,
                           Linear_Predictor_without_RE=C%*%bu+ + Y1[Year],
                           Linear_Predictor_with_RE=C%*%bu+ U[random] + Y1[Year],
                           p_with_RE=expit(C%*%bu+ U[random]+ Y1[Year]),
                           p_without_RE=expit(C%*%bu + Y1[Year]),
                           year=data$Year, Method="with yre",degree=4)
      datap <- datapp
      print(dim(datap))
    }else{
      datap <- data.frame(Depth= data$depth,
                          Linear_Predictor_without_RE=C%*%bu, Linear_Predictor_with_RE=C%*%bu+ U[random],
                          p_without_RE=expit(C%*%bu),p_with_RE=expit(C%*%bu+ U[random]),
                          year=year[indexY], Method="without yre",degree=4)
    }
    data_all <- rbind(data_all, datap)
  }  
  ## plot depth vs predictions for GAMM model
  pred <- names(data_all)[2:5]
  data_all_GAMM1 <- data_all %>%
    pivot_longer(cols = pred , names_to = "pred", values_to = "val")
  data_all_GAMM1$pred <- factor(data_all_GAMM1$pred,
                              levels = c("Linear_Predictor_without_RE","p_without_RE",
                                         "Linear_Predictor_with_RE", "p_with_RE"),
                              labels = c(bquote(bold(X*beta + Zu)),
                                         bquote(bold(logit^-1 * (X*beta + Zu))),
                                         bquote(bold(X*beta + Zu + Z^s*u^s)),
                                         bquote(bold(logit^-1 * (X*beta + Zu + Z^s*u^s)))))

  
  p <- ggplot(data = data_all_GAMM1,aes(x = Depth, y = val, col=Method, shape=Method, group=Method)) +
    geom_point(size=0.5) + theme_bw() +
    xlab("Depth") + ylab("Predictor") +
    scale_color_viridis_d(direction = -1) +
    theme(#panel.grid = element_line(color = "grey",size = 0.3,linetype = 1),
          #strip.background = element_rect(fill="grey",color="black"),
          strip.text = element_text(colour = "black",size = 14,face="bold"),
          strip.text.y = element_text(colour = "black",size = 11,face="bold",angle = 270),
          axis.title = element_text(color = "black",size=16,face="bold"),
          axis.text = element_text(color = "black", size=14),
          axis.text.x = element_text(color = "black", size=14,angle = 90),
          legend.title = element_text(size=16,face="bold"),
          legend.text = element_text(color = "black", size=16),
          legend.position = "bottom",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-5,-5,-5,-5),
          plot.margin = margin(t = 30, r = 30, b = 30, l = 30, unit = "pt")) + 
    guides(colour = guide_legend(override.aes = list(size=1.5))) +
    facet_grid(pred~year,scales = "free_y",labeller = labeller(
      pred = label_parsed
    ))
  return(p)

}

yvec <- c(1:4)
sum_data(pname="es_",yvec,  kp="es") + labs(title = "")
ggsave(paste0("plots/models/Ex2_depth_vs_predictor_GAMM.jpeg"),width = 2.5*3,height=2.3*4,bg="white")

