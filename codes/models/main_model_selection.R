## ---------------------------------------------------------------------------------------------
## Authors: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  Run models using Laplace approximation for model selection
## ---------------------------------------------------------------------------------------------

#################################  HPC-INDEX ###################################################
indexM <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))   ## index of the model
indexY <- as.integer(Sys.getenv("indexY"))            ## index for the data

################################ CODE #########################################################

library(HRW)
library(VGAM)
library(coda)
library(foreach)
library(doParallel)
library(LaplacesDemon)
library(plyr)
library(readxl)
library(mgcv)

source("codes/functions/basic_functions.R")
source("codes/functions/automated_Laplace_tp.R")

BE_data <- read_excel("data/BE_data/Barracouta_East_Rdata_FS.xlsx")[,-1]

covari <- c("asp", "depth", "plan", "prof", "slp", "hyp05", "hyp50", "rng05","rng50", "std05", "std50") ## all covariates
year <- c(2010,2011,2013,"yre")
yre <- c(rep(FALSE,3), TRUE)

BE_data <- read_excel("data/BE_data/Barracouta_East_Rdata_FS.xlsx")
if(year[indexY]=="yre"){
  data <- BE_data
}else{
  data <- subset(BE_data,BE_data$Year==year[indexY])
}

data[,covari] <- apply(data[,covari], 2, norm_01)
apply(data[,covari], 2, range)

FishNet <- unname(unlist(as.vector(data[,"FishNet50"])))
data$FishNet <- FishNet 

load(paste0("outputs/model_sets/models_",year[indexY],".RData"))
model <- models[indexM,2]

## log gamma priors on log precision
out_Lap <- Laplace_approx_gamm_OSullivan(f=models[indexM,2],
                                         data=data,B=500,alpha=FALSE,
                                         fre=TRUE,yre=yre[indexY],
                                         opt_method = "L-BFGS-B",
                                         prior="loggamma", beta_start = 0,
                                         HesP = TRUE, HesR = FALSE) 

model_evi <- Bayes_model_selec_lap(out_Lap[[1]])

if (!dir.exists(paste0("outputs/model_selection/",year[indexY],"/") ) ) {
  dir.create(paste0("outputs/model_selection/",year[indexY],"/"), recursive = TRUE, showWarnings = FALSE)
}

v_ut <- paste0("outputs/model_selection/",year[indexY],"/model_lg_",indexM)
v_fileName <- paste(v_ut,"RData",sep=".")
save(list=c("out_Lap","model_evi","model"),file=v_fileName)

