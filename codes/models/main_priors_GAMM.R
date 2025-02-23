## ---------------------------------------------------------------------------------------------
## Authors: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  Run GAMM models for data from year 2010, 2011, 2103, 2016 and all data 
##           with year random effect to use as priors for design
## ---------------------------------------------------------------------------------------------

################################ HPC-INDEX #####################################################
indexY <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))     # index for data
indexF <- as.numeric(Sys.getenv("indexF"))              # index for fishnet id
indexY <- 4
indexF <- 1

################################ CODE ##########################################################

library(NCmisc)
list.functions.in.file("codes/models/main_priors_GAMM.R", alphabetic = TRUE)

library(readxl)

source("codes/functions/basic_functions.R")
source("codes/functions/automated_Laplace.R")
source("codes/functions/cal_z.R")

set.seed(4654)

model <- as.formula(HardCoral ~  s(depth ,k=25))
year <- c(2010,2011,2013,"yre")
yre <- c(rep(FALSE,3), TRUE)

BE_data <- read_excel("data/BE_data/Barracouta_East_Rdata_FS.xlsx")
if(year[indexY]=="yre"){
  data <- BE_data
}else{
  data <- subset(BE_data,BE_data$Year==year[indexY])
}
range(data$depth)
unique(data$Year)
data$depth <- norm_01(data$depth)

fshname <- paste0("FishNet",seq(50,250,by=50))
FishNet <- unname(unlist(as.vector(data[,fshname[indexF]])))
data$FishNet <- FishNet 

## equally spaced
t1 <- Sys.time()
out_Lap <- Laplace_approx_gamm_OSullivan(f=model, data=data,B=1000,alpha=TRUE,
                                         fre=TRUE,yre=yre[indexY],prior="loggamma",
                                         ms=NULL,sds=NULL,mr=NULL,sdr=NULL,mt=NULL,sdt=NULL,
                                         HesP = TRUE, HesR = TRUE,opt_method = "L-BFGS-B",beta_start = 0,
                                         knot_type = "es", basis = "OS",
                                         intKnots=NA,boundKnots=NA)

out_Lap[[1]]$par

t2 <- Sys.time()
time <- difftime(t2, t1, units = "secs")

v_ut <- paste0("outputs/priors/GAMM/GAMM_prior_es_",fshname[indexF],"_",year[indexY])
v_fileName <- paste(v_ut,"RData",sep=".")
save(list=c("out_Lap","time"),file=v_fileName)

## quantiles
t1 <- Sys.time()
out_Lap <- Laplace_approx_gamm_OSullivan(f=model, data=data,B=1000,alpha=TRUE,
                                         fre=TRUE,yre=yre[indexY],prior="loggamma",
                                         ms=NULL,sds=NULL,mr=NULL,sdr=NULL,mt=NULL,sdt=NULL,
                                         HesP = TRUE, HesR = TRUE,opt_method = "L-BFGS-B",beta_start = 0,
                                         knot_type = "q", basis = "OS",
                                         intKnots=NA,boundKnots=NA)

out_Lap[[1]]$par

t2 <- Sys.time()
time <- difftime(t2, t1, units = "secs")

if (!dir.exists( paste0("outputs/priors/GAMM/") ) ) {
  dir.create(paste0("outputs/priors/GAMM/"), recursive = TRUE, showWarnings = FALSE)
}

v_ut <- paste0("outputs/priors/GAMM/GAMM_prior_q_",fshname[indexF],"_",year[indexY])
v_fileName <- paste(v_ut,"RData",sep=".")
save(list=c("out_Lap","time"),file=v_fileName)
