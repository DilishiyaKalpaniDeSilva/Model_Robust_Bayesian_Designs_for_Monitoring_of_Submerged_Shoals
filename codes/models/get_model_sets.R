# ---------------------------------------------------------------------------------------------
## Authors: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  Generate model sets with all the covariates for model comparison
## ---------------------------------------------------------------------------------------------

################################ CODE ##########################################################

library(readxl)
source("codes/functions/generate_model_sets.R")

BE_data <- data.frame(read_excel("data/BE_data/Barracouta_East_Rdata_FS.xlsx"))
colnames(BE_data)
range(BE_data$depth)

num.knots <- 25
limit <- c(0.5,-0.5)         ## maximum correlation between two covariates 
year <- c(2010,2011,2013,"yre")
covari <- c("asp", "depth", "plan", "prof", "slp", "hyp05", "hyp50", "rng05","rng50", "std05", "std50") ## all covariates
respon <- "HardCoral"

covari %in% colnames(BE_data)

## model sets when all covariates are considered as splines
covariates_lin <- NULL       ## linear terms in the model
covariates_spl <- covari     ## spline terms in the model

if (!dir.exists(paste0("outputs/model_sets/") ) ) {
  dir.create(paste0("outputs/model_sets/"), recursive = TRUE, showWarnings = FALSE)
}

for(indexY in 1:4){
  if(year[indexY]=="yre"){
    data <- BE_data
  }else{
    data <- subset(BE_data,BE_data$Year==year[indexY])
  }
  
  cor_data <- data[ , covari]
  models <- interac_models_1(limit,respon,covari, covariates_lin,covariates_spl,num.knots,cor_data,interac_type="te")
  save(list = c("models"),file = paste0("outputs/model_sets/models_",year[indexY],".RData"))
}
