## ---------------------------------------------------------------------------------------------
## Authors: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  Plot run-times of approximated expected utility under different numbers of design 
##          points: Supplement A, Figure S4
## ---------------------------------------------------------------------------------------------

################################ CODE ##########################################################

library(ggplot2)
library(dplyr)
library(tidyr)

B_all <- c(500,250)
B1_all <- c(500,250)
i <- 1
Design_no <- B_vec <- B1_vec <- N_all <- timeR_vec <- vector()
year <- 2010

for(indexD in 1:50){
  print(indexD)
  for (indexB in 1:2) {
    for(indexB1 in 1:2){
      for(indexT in 1:9){
        print("*******")
        B_vec[i] <- B_all[indexB]
        B1_vec[i] <- B1_all[indexB1]
        Design_no[i] <- indexD
        N_all[i] <- indexT*2*50
        
        path1 <- paste0("outputs/comparisons/run_time/",year,"/Design","_D",indexD,"_B_",B_vec[i],
                        "_B1_",B_vec[i],"_N",N_all[i],".RData")
        if(file.exists(path1)){
          load(path1)
          timeR_vec[i] <- time
          print(time)
        }else{
          timeR_vec[i] <- NA
        }
        i <- i + 1  
      }
    }    
  }
}

time_and_util_data <- data.frame(Design_no , L = paste0("L = ",B_vec) , E = paste0("E = ",B1_vec) , 
                                 N = N_all , timeR =  timeR_vec/60 )
time_and_util_data <- na.omit(time_and_util_data)
# save(list = c("time_and_util_data"), file = "summary/time_and_util_data.RData")

mean_data <- time_and_util_data %>%
  group_by(L,E,N) %>%
  dplyr::summarise(mean_timeR = mean(timeR, na.rm=TRUE))

mean_data %>% 
  tidyr::pivot_longer(c("mean_timeR"),names_to="method",values_to ="mean_time") %>%
  ggplot(.,aes(x=N,y=mean_time)) +
  xlab("Number of design points") + ylab("Average run time per design (mins)")+
  geom_line()+
  geom_point() +
  theme_bw()+
  scale_color_viridis_d() +
  facet_wrap(L~E, scales = "free_y")+
  scale_x_continuous(breaks = seq(from=100,to=900, by=100)) +
  theme(strip.text.y = element_text(angle = 270),
        #strip.background = element_rect(fill="white",color="black"),
        strip.text = element_text(colour = "black",size = 12,face="bold"),
        title = element_text(size=12,face="bold"),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(color = "black",size=12,face="bold"),
        axis.text.x = element_text(color = "black", size=12,angle = 90),
        axis.text.y = element_text(color = "black", size=12,angle = 0),
        legend.text = element_text(color = "black", size=12),
        legend.position = "bottom",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-3,-5,-5))


time_and_util_data %>% 
  ggplot(aes(x=N, y=timeR)) +
  xlab("Number of design points (n)") + ylab("Run-time (mins)") +
  geom_line(data = mean_data, aes(x = N, y = mean_timeR, group = interaction(L, E), color = "Mean"), size=1) +
  geom_point(aes(color="Each Design")) +
  theme_bw() +
  scale_color_manual(values = c("Mean" = "#440154FF", "Each Design" = "#FDE725FF"),
                     name = "Run-time",
                     labels = c("Each Design", "Mean")) +
  facet_grid(L~E, scales = "free_y") +
  scale_x_continuous(breaks = seq(from=100, to=900, by=100)) +
  theme(strip.text.y = element_text(angle = 270),
        strip.text = element_text(colour = "black", size = 12, face = "bold"),
        title = element_text(size=12, face="bold"),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(color = "black", size=12, face="bold"),
        axis.text.x = element_text(color = "black", size=12, angle = 90),
        axis.text.y = element_text(color = "black", size=12, angle = 0),
        legend.text = element_text(color = "black", size=12),
        legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-5, -3, -5, -5))

ggsave(paste0("plots/comparisons/Ex2_run_time.jpeg"),width = 3*3,height=2*3)
