##########################################  HPC-INDEX  ##########################################

indexL <- 1
indexY <- as.integer(Sys.getenv("indexY"))
indexO <- as.integer(Sys.getenv("indexO"))
iter <- as.integer(Sys.getenv("iter"))
indexTr <- as.integer(Sys.getenv("indexTr"))
indexP <- as.integer(Sys.getenv("indexP"))
indexA <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
B <- as.integer(Sys.getenv("B"))            # no of MC draws for the utility approximation
B1 <- as.integer(Sys.getenv("B1"))          # no of MC draws in Laplace approximation 
ncl <- as.integer(Sys.getenv("ncl"))
indexI <- as.integer(Sys.getenv("indexI"))

##########################################  CODE  ################################################

library(gtools) # to normalize data
library(HRW) #O'Sullivan spline design matrices
library(mvtnorm) # to generate spatial random effect
library(LaplacesDemon)
library(foreach)
library(doParallel)
library(dplyr)
library(Matrix)
library(plyr)
library(readxl)
library(DiceKriging)

source("codes/functions/cal_z.R")
source("codes/functions/basic_functions.R")
source("codes/functions/Rcpp_functions.R")
source("codes/designs/polynomial/with_yre/approximate_utility_pol.R")

set.seed(2008)

year_all <- "yre"
year <- year_all[indexY]
print(year)

tlength_all <- c(500,1000,1500)
tlength <- tlength_all[indexL]
print(tlength)

ntrans <- 18
print(ntrans)

print(B)# no of MC for the utility approximation
print(B1) # no of MC in Laplace approximation 

fol <- paste0("designs/polynomial/",year,"/",indexO,"/",indexI,"/")
fol1 <- paste0("outputs/",fol,"Iteration_",iter,"/Transect_",indexTr,"/")
if (!dir.exists(fol1)) {
  dir.create(fol1, recursive = TRUE)
}

cl <- makePSOCKcluster(ncl,setup_timeout=2400)
clusterCall(cl, RcppFunctions)
registerDoParallel(cl)

load(paste0("outputs/priors/polynomial/poly_prior_lg_FishNet50_",year,"_deg",indexO,".RData")) ## loading priors
load(paste0("data/transect_data/data_for_design_fshsize_predang100_n50_dep50_len",tlength,".RData"))
fishnetp <- fishnet_all[[1]]

data <- read_excel(paste0("data/BE_data/BE_data_",year,".xlsx"))
range(data$depth)
FishNet <- unname(unlist(as.vector(data[,"FishNet50" ])))
data$FishNet <- FishNet

## to normalize the data according to original depth range from different years
depth_range <- range(data$depth)
X_all <- list()
for(i in 1:length(depthp)){
  X_all[[i]] <- scaleP(depthp[[i]],depth_range[1],depth_range[2],-1,1)
  # print(range(X_all[[i]], na.rm = TRUE))
}
print(range(X_all, na.rm = TRUE))

data$depth <- scaleP(data$depth,depth_range[1],depth_range[2],-1,1)
range(data$depth)
Xd_hist <- cbind( 1,outer(X = data$depth,Y = 1:indexO,FUN = "^"))

# mu and covariate matrix for all both parameters and u
fitOptim_prior <- out_Lap[[1]]
fitOptim_prior_u <- out_Lap[[2]]

# posterior from the model as prior for the design for beta and log(1/(sigma_r^2))
mu_prior <- fitOptim_prior$par
sigma_prior <- solve(-fitOptim_prior$hessian)

para_path <- paste0("codes/designs/polynomial/para_all_",year,"_",indexO,".RData")
load(para_path)

if(indexTr==1){ ## initial design at the start of each iteration
  if(iter==1){ ## initializing design when iteration=1
    iname <- paste0("outputs/",fol,"Iteration_",0,".RData")
    if(!file.exists(iname)){
      curr_u <- mean(exp_crit(d=tr_design,B),na.rm=TRUE)
      save(list = c("tr_design","curr_u"), file = iname)
    }else{
      load(iname)
    }
  }else{ ## initializing design when iteration>1, saved in compare utility iteration
    load(paste0("outputs/",fol,"Iteration_",iter-1,".RData"))
  }
}else{ ## rest of the transect replacements
  load(paste0("outputs/",fol,"Iteration",iter,"_transect_opt_",indexTr-1,".RData")) 
}

# calculating utility when each of 20 transects are replaced with the t transect of the design
Q1 <- 20
tr_no <- (5*(Q1*(indexP-1) + indexA))-4 
tr_design[indexTr] <- tr_no
print(tr_design)
uz_allA <- exp_crit(d=tr_design,B)
uzA <- mean(uz_allA,na.rm=TRUE)

fnameA <- paste0(fol1,"point",indexP,"_angle",indexA,".RData")
save(list = c("tr_design","uz_allA","uzA"), file = fnameA)

stopCluster(cl)

