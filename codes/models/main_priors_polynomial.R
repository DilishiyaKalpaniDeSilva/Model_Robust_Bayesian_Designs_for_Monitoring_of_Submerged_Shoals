## ---------------------------------------------------------------------------------------------
## Authors: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  Run polynomial models for data from year 2010, 2011, 2103 and all data 
##           with year random effect to use as priors for design
## ---------------------------------------------------------------------------------------------

################################ HPC-INDEX #####################################################
indexY <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))  ## sampling year of the historical data
indexO <- as.integer(Sys.getenv("indexO"))           ## order of the polynomial
indexF <- as.integer(Sys.getenv("indexF"))           ## index for different fishnet sizes

################################ CODE ##########################################################

library(NCmisc)
list.functions.in.file("codes/models/main_priors_polynomial.R", alphabetic = TRUE)

library(readxl)

source("codes/functions/cal_z.R")
source("codes/functions/basic_functions.R")
source("codes/functions/automated_Laplace.R")

set.seed(4654)

models <- list()
models[[1]] <- as.formula(HardCoral ~  poly(depth ,1))
models[[2]] <- as.formula(HardCoral ~  poly(depth ,2))
models[[3]] <- as.formula(HardCoral ~  poly(depth ,3))

year <- c(2010,2011,2013,"yre")
yre <- c(rep(FALSE,3), TRUE)

BE_data <- read_excel("data/BE_data/Barracouta_East_Rdata_FS.xlsx")
if(year[indexY]=="yre"){
  data <- BE_data
}else{
  data <- subset(BE_data,BE_data$Year==year[indexY])
}

range(data$depth)
data$depth <- scaleP(data$depth,min(data$depth),max(data$depth),-1,1)
range(data$depth)
unique(data$Year)

fshname <- paste0("FishNet",seq(50,250,by=50))
FishNet <- unname(unlist(as.vector(data[,fshname[indexF]])))
data$FishNet <- FishNet 

## log gamma priors on log precision
t1 <- Sys.time()
out_Lap <- Laplace_approx_gamm_OSullivan(f=models[[indexO]], data=data,B=1000,alpha=TRUE,
                                         fre=TRUE,yre=yre[indexY],prior="loggamma",
                                         ms=NULL,sds=NULL,mr=NULL,sdr=NULL,mt=NULL,sdt=NULL,
                                         HesP = TRUE, HesR = TRUE,opt_method = "L-BFGS-B",beta_start = 0.001,
                                         knot_type = NA, basis = NA,
                                         intKnots=NA,boundKnots=NA)
t2 <- Sys.time()
time <- difftime(t2, t1, units = "secs")

if ( !dir.exists(paste0("outputs/priors/polynomial/") ) ) {
  dir.create(paste0("outputs/priors/polynomial/"), recursive = TRUE, showWarnings = FALSE)
}
v_ut <- paste0("outputs/priors/polynomial/poly_prior_lg_",fshname[indexF],"_",year[indexY],"_deg",indexO)
v_fileName <- paste(v_ut,"RData",sep=".")
save(list=c("out_Lap","time"),file=v_fileName)

