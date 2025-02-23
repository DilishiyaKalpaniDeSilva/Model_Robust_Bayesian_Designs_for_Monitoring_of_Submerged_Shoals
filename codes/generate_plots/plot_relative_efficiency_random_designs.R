## ---------------------------------------------------------------------------------------------
## Authors: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  Plot relative utility of 100 random designs compared to optimal design: 
##          Supplement A, Figure S3
## ---------------------------------------------------------------------------------------------

##################################### CODE #####################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(matrixcalc)

dtype <- c("shallow","deep","both","clustered")
fshall <- paste0("FishNet",seq(50,250,by=50))

time1  <- vector()
mu1 <- sig1 <- list()
kld1 <- kld_betam <- kld_beta10 <- kld_beta11 <- kld_lprec_r1 <- kld_lprec_s1 <- vector()

Design_no <-deg <- Year <- fsize <- type <- vector()
k <- 1

dname="es"

for (year in c(2010)) {
  for (indexDC in 1:4) {
    for(indexD in 1:25){
      
      Design_no[k] <- (indexDC-1)*25 + indexD
      type[k] <- dtype[indexDC]
      Year[k] <- year
      
      mu <- sig <- kld <-  kld_beta <- kld_beta0 <- kld_beta1 <- kld_lprec_r <- kld_lprec_s <- vector()
      
      iname <- paste0("outputs/comparisons/exp_utility/GAMM/",year,"/Design_",dname,"_",dtype[indexDC],"_D",indexD,".RData")
      
      if(file.exists(iname)){
        load(iname)
        #print(time/60)
        print(k)
        
        for (i in 1:length(util)) {
          kld[i] <- util[[i]][[3]]
          kld_beta0[i] <- util[[i]][[4]][1]
          kld_beta1[i] <- util[[i]][[4]][2]
          kld_beta[i] <- util[[i]][[5]]
          kld_lprec_s[i] <- util[[i]][[6]]
          kld_lprec_r[i] <- util[[i]][[7]]
        }
        kld1[k] <- mean(unlist(kld))
        kld_beta10[k] <- mean(unlist(kld_beta0))
        kld_beta11[k] <- mean(unlist(kld_beta1))
        kld_betam[k] <-  mean(unlist(kld_beta))
        kld_lprec_s1[k] <- mean(unlist(kld_lprec_s))
        kld_lprec_r1[k] <- mean(unlist(kld_lprec_r))
        time1[k] <- time
      }else{
        kld1[k] <- NA
        kld_beta10[k] <- NA
        kld_beta11[k] <- NA
        kld_betam[k] <-  NA
        kld_lprec_s1[k] <- NA
        kld_lprec_r1[k] <- NA
      }
      k <- k+1
    }
  }    
}

util_data <- data.frame(Design_no , kld1 , kld_betam,kld_beta10,kld_beta11,kld_lprec_s1,kld_lprec_r1, type, Year)
util_data$Design_no <- 1:100
#util_data1 <- na.omit(util_data)

util_data11 <- util_data %>% 
  pivot_longer(cols = c("kld1", "kld_betam","kld_beta10","kld_beta11","kld_lprec_s1","kld_lprec_r1" ),
               names_to = "Util_Cat",values_to = "Exp_util")

rel_data <- util_data %>% 
  pivot_longer(cols = c("kld1" ),
               names_to = "Util_Cat",values_to = "Exp_util")

opt_util <- 134.11
rel_data$Rel_util <- rel_data$Exp_util/opt_util

ggplot(rel_data,aes(x=factor(Design_no),y=Rel_util,colour=type,shape=type, group=Util_Cat)) +
  xlab("Design Number") + ylab("Relative Efficiency")+
  ylim(0,1) +
  theme_bw() +
  geom_point(size=1.8) +
  scale_x_discrete(breaks=seq(0,100,5)) +
  labs(shape="Design category",color="Design category") +
  scale_color_viridis_d(direction = -1) +
  theme(legend.position = "bottom",
        axis.title = element_text(color = "black",size=12,face="bold"),
        legend.text = element_text(color = "black", size=12),
        legend.title = element_text(color = "black",size=12,face="bold"),
        strip.text = element_text(color = "black",size=12,face="bold")) +
  facet_grid(~Year)

ggsave(paste0("plots/comparisons/Ex2_rel_util_rand_designs.jpeg"),width = 3*3,height=2*3)

