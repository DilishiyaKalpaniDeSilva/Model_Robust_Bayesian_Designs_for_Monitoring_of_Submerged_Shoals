##########################################  HPC-INDEX  ##########################################

indexL <- 1
indexY <- as.integer(Sys.getenv("indexY"))
B <- as.integer(Sys.getenv("B"))            
B1 <- as.integer(Sys.getenv("B1"))          # no of MC draws in Laplace approximation 
ncl <- as.integer(Sys.getenv("ncl"))
indexI <- as.integer(Sys.getenv("indexI"))

##########################################  CODE  ################################################

library(gtools)      # to normalize data
library(HRW)         # O'Sullivan spline design matrices
library(mvtnorm)     # to generate spatial random effect
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
source("codes/designs/GAMM/with_yre/approximate_utility_GAMM.R")

set.seed(2008)

year_all <- c("yre") 
year <- year_all[indexY]
print(year)

tlength_all <- c(500,1000,1500)
tlength <- tlength_all[indexL]
print(tlength)

ntrans <- 18
print(ntrans)

print(B)        # no of MC draws for the utility approximation
print(B1)       # no of MC draws in Laplace approximation 

cl <- makePSOCKcluster(ncl,setup_timeout=2400)
clusterCall(cl, RcppFunctions)
registerDoParallel(cl)

load(paste0("outputs/priors/GAMM/GAMM_prior_es_FishNet50_",year,".RData")) ## loading priors
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
  X_all[[i]] <- norm_01_mm(depthp[[i]],depth_range[1], depth_range[2])
}

data$depth <- norm_01(data$depth)
range(data$depth)
Xd_hist <- cbind(1,data$depth)

# mu and covariate matrix 
fitOptim_prior <- out_Lap[[1]]
fitOptim_prior_u <- out_Lap[[2]]
knot_loc <- out_Lap[[3]]
num.knots <- k <- length(knot_loc)
intKnots <- knot_loc[-c(1, length(knot_loc))]

zz <- cal_z(numBasis = num.knots,XsplPreds=as.matrix(data$depth),knot_type="es" ,intKnots,
            boundKnots=range(data$depth), basis ="OS")
z.spline <- zz$Z

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

# prior from running model for the previous data
mu_prior_U2 <- fitOptim_prior_u$par[(k+1):(k+nnf)]
sigma_prior_U2 <- solve(-fitOptim_prior_u$hessian)[(k+1):(k+nnf),(k+1):(k+nnf)]

mu_prior_Uall <- rep(mu_prior_U1,num_fish_all)
sigma_prior_Uall <- diag(rep(sigma_prior_U1^2,num_fish_all))

fsh_d <- unique(data$FishNet)
mu_prior_Uall[fsh_d] <- mu_prior_U2
sigma_prior_Uall[fsh_d,fsh_d] <- sigma_prior_U2

## temporal random effect
num_year <- length(unique(data$Year))
mu_prior_t <- 0
sigma_prior_t <- unname(1/sqrt(exp(fitOptim_prior$par["log_prec_t"])))

## to integrate out in the full data likelihood
rand_u <- matrix(runif(B1*num.knots),num.knots,B1)
rand_U <- matrix(runif(B1*num_fish_all),num_fish_all,B1) 
rand_t <- matrix(runif(B1*num_year),num_year,B1) 

## set up for data generation in utility approximation (should be posterior from historical data)
# model parameters 
theta_prior <- rmvnorm(B, mu_prior, sigma_prior)  #### Should be posterior from historical data
u_all <- rmvnorm(B, mu_prior_u, sigma_prior_u)
u_all <- as(u_all, "dgCMatrix")
U_all <- rmvnorm(B, mu_prior_Uall, sigma_prior_Uall)
U_all <- as(U_all, "dgCMatrix")
t_all <- rmvnorm(B, mu_prior_t, as.matrix(sigma_prior_t^2))
t_all <- as(t_all, "dgCMatrix")

## to get deterministic y values
n_all <- 50*ntrans #+ nrow(data)
rand_y <- matrix(runif(n_all*B),nrow = B , ncol = n_all, byrow = TRUE)

init_tr_designs <- lapply(1:10, function(x) sample(nt, ntrans, replace = FALSE))
tr_design <- init_tr_designs[[indexI]]

para_path <- paste0("codes/designs/GAMM/para_all_",year,".RData")
save(list=c("rand_U","rand_u","rand_y","rand_t","theta_prior","u_all","U_all","t_all","tr_design"), file=para_path)
# load(para_path)

fol <- paste0("designs/GAMM/",year,"/",indexI,"/")
fol1 <- paste0("outputs/",fol)
if (!dir.exists(fol1)) {
  dir.create(fol1, recursive = TRUE)
}

curr_u <- mean(exp_crit(d=tr_design,B),na.rm=TRUE)
iname <- paste0("outputs/",fol,"Iteration_",0,".RData")
save(list = c("tr_design","curr_u"), file = iname)

stopCluster(cl)

