##########################################  HPC-INDEX  ##########################################

indexL <- 1
indexY <- as.integer(Sys.getenv("indexY"))
iter <- as.integer(Sys.getenv("iter"))
indexTr <- as.integer(Sys.getenv("indexTr"))
B <- as.integer(Sys.getenv("B"))            # no of MC draws for the utility approximation
B1 <- as.integer(Sys.getenv("B1"))          # no of MC draws in Laplace approximation 
ncl <- as.integer(Sys.getenv("ncl"))
indexI <- as.integer(Sys.getenv("indexI"))
P <- 22

##########################################  CODE  ################################################

library(foreach)
library(doParallel)

set.seed(2008)

year_all <- c("yre") 
year <- year_all[indexY]

tlength_all <- c(500,1000,1500)
tlength <- tlength_all[indexL]

fol<-paste0("designs/GAMM/",year,"/",indexI,"/")

if(indexTr==1){ ## initial design at the start of each iteration
  iname <- paste0("outputs/",fol,"Iteration_",iter-1,".RData")
  load(iname)
}else{ ## rest of the transect replacements
  load(paste0("outputs/",fol,"Iteration",iter,"_transect_opt_",indexTr-1,".RData"))
}

u_maxO <- curr_u

u <- tr <- vector()
for(indexP in 1:P){
  fnamep <- paste0(paste0("outputs/",fol,"/Iteration_",iter,"/Transect_",indexTr,"/point_",indexP),".RData")
  if(file.exists(fnamep)){
    load(fnamep)
  }else{
    print(indexP)
    source("codes/designs/GAMM/with_yre/6_optimize_angle_ACE_rest.R")
    load(fnamep)
  }
  u[indexP] <- outACE$curr_u
  tr[indexP] <- outACE$curr_t
}

ind_c <- which(u==max(u,na.rm = TRUE))#index of maximum utility
u_maxN <- u[ind_c[1]] #x value which maximizes the utility, for each x value
if(u_maxN > u_maxO){
  tr_design[indexTr] <- tr[ind_c[1]]
  curr_u <- u[ind_c[1]]
}

pname <- paste0("outputs/",fol,"Iteration",iter,"_transect_opt_",indexTr,".RData")
save(list = c("tr_design","curr_u"), file = pname)


