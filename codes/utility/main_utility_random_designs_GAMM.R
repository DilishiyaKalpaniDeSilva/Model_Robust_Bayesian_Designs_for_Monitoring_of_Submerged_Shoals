## ---------------------------------------------------------------------------------------------
## Authors: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  Evaluate expected utility of 100 random designs 
## ---------------------------------------------------------------------------------------------

################################ HPC-INDEX #####################################################
indexD <- as.numeric(Sys.getenv("indexD"))          ## design number
B <- as.numeric(Sys.getenv("B"))                    ## Monte Carlo draws for utility
B1 <- as.numeric(Sys.getenv("B1"))                  ## Monte Carlo draws for Laplace
indexY <- as.numeric(Sys.getenv("indexY"))          ## prior
indexDC <- as.numeric(Sys.getenv("indexDC"))        ## depth category
indexF <- 1                                         ## fishnet size
indexL <- 1
ncl <- 12

################################ CODE ##########################################################

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

source("codes/functions/cal_z.R")
source("codes/functions/basic_functions.R")
source("codes/functions/Rcpp_functions.R")
source("codes/utility/approximate_utility_GAMM.R")

set.seed(2008)

year_all <- c(2010,2011,2013,"yre") 
year <- year_all[indexY]
print(year)

tlength_all <- c(500,1000,1500)
tlength <- tlength_all[indexL]
print(tlength)

ntrans <- 18
print(ntrans)

print(B)# no of MC for the utility approximation
print(B1) # no of MC in Laplace approximation 

Sys.sleep((indexD+indexDC)*2)

cl <- makePSOCKcluster(ncl)
clusterCall(cl, RcppFunctions)
registerDoParallel(cl)

load(paste0("outputs/priors/GAMM/GAMM_prior_es_FishNet50_",year,".RData")) ## loading priors
load(paste0("data/transect_data/random_designs_fs.RData"))
load(paste0("data/transect_data/data_for_design_fshsize_predang100_n50_dep50_len",tlength,".RData"))
fishnetp <- fishnet_all[[1]]

data <- read_excel(paste0("data/BE_data/BE_data_",year,".xlsx"))
data <- subset(data, data$depth >= (-52))
range(data$depth)
FishNet <- unname(unlist(as.vector(data[,"FishNet50" ])))
data$FishNet <- FishNet

## to normalize the data according to original depth range from different years
depth_range <- range(data$depth)
X_all <- list()
range(depthp, na.rm = TRUE)
for(i in 1:length(depthp)){
  X_all[[i]] <- norm_01_mm(depthp[[i]],depth_range[1], depth_range[2])
}

data$depth <- norm_01(data$depth)
range(data$depth)
Xd_hist <- cbind(1,data$depth)

for (i in 1:100) {
  d <- design[[i]]
  Xd <- unlist(X_all[d]) 
  print(range(Xd))
  if(range(Xd)[2]>range(data$depth)[2]){
    print(i)
  }
  
}

# mu and covariate matrix for all both parameters and u
fitOptim_prior <- out_Lap[[1]]
fitOptim_prior_u <- out_Lap[[2]]
knot_loc <- out_Lap[[3]]
num.knots <- k <- length(knot_loc)
intKnots <- knot_loc[-c(1, length(knot_loc))]

zz <- cal_z(numBasis = num.knots,XsplPreds=as.matrix(data$depth),knot_type="es" ,intKnots,
            boundKnots=range(data$depth), basis ="OS")
z.spline <- zz$Z
C <- cbind(1,data$depth,z.spline)
bu <- c(out_Lap[[1]]$par[1:2],out_Lap[[2]]$par[1:num.knots])
plot(data$depth, (C%*% bu))

# posterior from the model as prior for the design for beta and log(1/(sigma_r^2))
mu_prior <- fitOptim_prior$par
sigma_prior <- solve(-fitOptim_prior$hessian)

# posterior from the model as prior for the design for wiggliness parameter u
mu_prior_u <- fitOptim_prior_u$par[1:k]
sigma_prior_u <- solve(-fitOptim_prior_u$hessian)[1:k,1:k]

## remove empty transects which are outside the shoal
j <- 1; ind <- vector()
for(i in 1:length(fishnetp)){
  if(length(fishnetp[[i]])==1){
    ind[j] <- i
    j <- j+1
  }
}

nt <- (1:length(fishnetp))[-ind]
num_fish_all1 <- range(unique(unlist(fishnetp[nt]))) ## total number of different fishnets for designs
length(unique(unlist(fishnetp[nt])))

num_fish_all2 <- range(unique(data$FishNet)) ## total number of different fishnets AIMS data
length(unique(data$FishNet))
length(fitOptim_prior_u$par)-k
nnf <- length(unique(data$FishNet))
num_fish_all <- max(c(num_fish_all1,num_fish_all2))

# as the fishnet random effect values are unknown for new location
mu_prior_U1 <- 0 
sigma_prior_U1 <- unname(1/sqrt(exp(fitOptim_prior$par["log_prec_r"])))
#sigma_prior_U1 <- unname(exp(fitOptim_prior$par["log_sigma_r"]))

# prior from running model for the previous data
mu_prior_U2 <- fitOptim_prior_u$par[(k+1):(k+nnf)]
sigma_prior_U2 <- solve(-fitOptim_prior_u$hessian)[(k+1):(k+nnf),(k+1):(k+nnf)]

mu_prior_Uall <- rep(mu_prior_U1,num_fish_all)
sigma_prior_Uall <- diag(rep(sigma_prior_U1^2,num_fish_all))
# sigma_prior_Uall <- as(sigma_prior_Uall, "dgCMatrix")

fsh_d <- unique(data$FishNet)
mu_prior_Uall[fsh_d] <- mu_prior_U2
sigma_prior_Uall[fsh_d,fsh_d] <- sigma_prior_U2

## to integrate out in the full data likelihood
rand_U <- matrix(runif(B1*num_fish_all),num_fish_all,B1) # generating probabilities
rand_u <- matrix(runif(B1*num.knots),num.knots,B1)

## set up for data generation in utility approximation (should be posterior from historical data)
# model parameters 
theta_prior <- rmvnorm(B, mu_prior, sigma_prior)  #### Should be posterior from historical data
u_all <- rmvnorm(B, mu_prior_u, sigma_prior_u)
u_all <- as(u_all, "dgCMatrix")

U_path <- "codes/utility/U_all.RData"
if(file.exists(U_path)){
  load(U_path)
}else{
  U_all <- rmvnorm(B, mu_prior_Uall, sigma_prior_Uall)
  U_all <- as(U_all, "dgCMatrix")
  save(list=c("U_all"), file=U_path)
}

n_all <- 50*ntrans #+ nrow(data)
rand_y <- matrix(runif(n_all*B),nrow = B , ncol = n_all, byrow = TRUE)

# design
dtype <- c("shallow","deep","both","clustered")
dc <- 25*(indexDC-1) + indexD
tr_design <- design[[dc]]

t1 <- Sys.time()
util <- exp_crit(d=tr_design,B)
t2 <- Sys.time()
time <- difftime(t2, t1, units = "secs")

if (!dir.exists(paste0("outputs/comparisons/exp_utility/GAMM/",year,"/")) ) {
  dir.create(paste0("outputs/comparisons/exp_utility/GAMM/",year,"/"), recursive = TRUE)
}

iname <- paste0("outputs/comparisons/exp_utility/GAMM/",year,"/Design_","es","_",dtype[indexDC],"_D",indexD,".RData")
save(list = c("tr_design","util","time"), file = iname)

stopCluster(cl)

