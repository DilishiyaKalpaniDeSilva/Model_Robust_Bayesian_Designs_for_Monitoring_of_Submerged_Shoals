## ---------------------------------------------------------------------------------------------
## Authors: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  Plot relative efficiency  
## ---------------------------------------------------------------------------------------------

##################################### CODE #####################################################

library(ggplot2)
library(raster)
library(rasterVis)
library(rgdal)
library(readxl)

True_Model_all <- c("Linear", "Quadratic", "Cubic")
year_all <- c(2010, 2011, 2013, "yre")
Assumed_Model_all <-  c("Linear"  , "Quadratic" ,"Cubic" ,"GAMM")

data_all <- NULL
for (year in year_all) {
  for (model in True_Model_all) {
    path <- paste0("outputs/rel_utility/rel_util_", year, "_pol_", model, ".RData")
    if (file.exists(path)) {
      load(path)
      for(i in 1:length(Assumed_Model_all)){
        data <- data.frame(Prior=year, True_Model=model, Assumed_Model=Assumed_Model_all[i], 
                           Relative_Efficiency=rel_util[i])  
        data_all <- rbind(data_all, data)
      }
    }else {
      for(i in 1:length(Assumed_Model_all)){
        data <- data.frame(Prior=year, True_Model=model, Assumed_Model=Assumed_Model_all[i], 
                           Relative_Efficiency=NA) 
        print("not found")
      }
    }
  }
}

data_all$Assumed_Model<-factor(data_all$Assumed_Model,
                               levels = c("Linear"  , "Quadratic" ,"Cubic" ,"GAMM"),
                               labels=c("Linear"  , "Quadratic" ,"Cubic" ,"GAMM"))

data_all$True_Model<-factor(data_all$True_Model,
                            levels = c("Linear"  , "Quadratic" ,"Cubic" ),
                            labels=c("Linear"  , "Quadratic" ,"Cubic" ))

ggplot(data_all,aes(x=True_Model, y=Relative_Efficiency,
                          colour=Assumed_Model,group=Assumed_Model,
                          shape=Assumed_Model)) +
  theme_bw() + 
  geom_point(size=2)+ 
  geom_line(size=1)+
  scale_color_viridis_d(direction = 1)+
  xlab("True Model") + ylab("Relative Efficiency")+
  facet_wrap(~ Prior, ncol=4) + 
  theme(strip.text.y = element_text(angle = 270),
        #strip.background = element_rect(fill="white",color="black"),
        strip.text = element_text(colour = "black",size = 12,face="bold"),
        title = element_text(size=12,face="bold"),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(color = "black",size=12,face="bold"),
        axis.text.x = element_text(color = "black", size=10,angle = 0),
        axis.text.y = element_text(color = "black", size=12,angle = 0),
        legend.text = element_text(color = "black", size=10),
        legend.position = "bottom",
        legend.margin=margin(0,0,0,0)) + 
  labs(color="Assumed Model",shape="Assumed Model")

ggsave(paste0("plots/designs/Ex2_relative_efficiency.jpeg"),width = 2*4,height=3,bg="white")
