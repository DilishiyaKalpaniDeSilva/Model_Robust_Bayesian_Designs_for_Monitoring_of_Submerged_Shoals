# ---------------------------------------------------------------------------------------------
## Authors: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Last updated: 26 May 2023 - Finalised
## ---------------------------------------------------------------------------------------------
## Purpose:  Generate model sets with all the covariates for model comparison
## ---------------------------------------------------------------------------------------------

################################ HPC-INDEX #####################################################

# library(NCmisc)
# ## check the packages used in the R script
# list.functions.in.file("codes/models/model_selection/generate_model_sets.R", alphabetic = TRUE) 

library(gtools)
library(stringr)
library(data.table)

source("codes/functions/basic_functions.R")

interac_models_1 <- function(limit,respon, covari, covariates_lin, covariates_spl,num.knots,cor_data,interac_type){
  
  covariates_s <- c(covariates_lin,paste0("s(",covariates_spl,", k=",num.knots,")")) ## all the terms
  cor_data <- apply(cor_data, 2, norm_01)
  core <- cor(cor_data) # correlation between covariates
  com_pair <- combinations(n=length(covari),r=2,v=covari,repeats.allowed=F) ## all possible combinations of two covariates
  cor_pair <- lapply(seq_len(nrow(com_pair)), function(x){core[com_pair[x,1],com_pair[x,2]]}) # correlation between each pair
  high_cor <- com_pair[which(cor_pair>=limit[1] | cor_pair<=limit[2]),] ## covariate pairs which exceed correlation limits
  
  ## all the possible models
  com_s <- lapply(seq_len(length(covariates_s)),
                  FUN = function(x)
                  {
                    c <- combinations(n=length(covariates_s),r=x,v=covariates_s,repeats.allowed=F)
                    f1 <-lapply(seq_len(nrow(c)),
                                FUN = function(y)
                                  paste(paste(respon, " ~ ", paste(c[y,], collapse= " + ")))
                                
                    )}
  )
  models <- unlist(com_s)
  
  ## searching index of the models that contains high correlated pairs for each high_corr pair
  models1 <- lapply(seq_len(nrow(high_cor)), function(x){
    patern <- paste(high_cor[x,],collapse = "|")
    more <- models[str_count(models, patern) > 1]
    w <- which(str_count(models, patern) > 1)
    #str_view_all(more, patern)
    return(w)
  })
  
  models2 <- unlist(models1)
  models_hc <- sort(unique(models2))
  models_wi <- models[-models_hc] ##  removing models with high correlated pairs
  
  # number of terms in a model
  terms <- str_count(models_wi, "\\+")
  table <- data.table(table(terms))
  freq <- table$N
  
  ## adding tensor product if both covariates are present , number of terms greater than 1
  terms1 <- lapply(c((freq[1]+1):length(models_wi)), function(x){
    f <- as.formula(models_wi[x])
    d <- all.vars(f)[-1]
    com1 <- list()
    com1 <- (combinations(n=length(d),r=2,v=d,repeats.allowed=F)) 
    
    int <- lapply(seq_len(nrow(com1)), function(y){
      e <- paste0(paste0(interac_type,"(", paste0(com1[y,], collapse= ", ")),", k=",num.knots,")")
    })
    
    if(length(int)>1){
      lapply(seq_len(nrow(com1)), function(z){
        com2 <- combinations(n=(length(int)),r=z,v=unlist(int),repeats.allowed=F)
        lapply(seq_len(nrow(com2)),function(w)
          f <- paste(paste(models_wi[x]," + "),paste(com2[w,], collapse=" + "))
        )
      })
    }else{
      f <- paste(models_wi[x],unlist(int),sep =" + ")
    }
  })
  models_int <- unlist(terms1)
  print(c(length(models_wi),length(models_int)))
  nm <- length(models_wi)+ length(models_int)
  models_all <- data.frame(index =c(1:nm),model=c(models_wi,models_int))
  return(models_all)
}

## model sets when slp is considered as a linear term
# covariates_lin <- "slp"          ## linear terms in the model
# covariates_spl <- covari[-5]     ## spline terms in the model
# 
# for(indexY in 1:6){
#   if(indexY==5){
#     data <- subset(BE_data,BE_data$Year!=2010)
#   }else if(indexY==6){
#     data <- BE_data
#   }else{
#     data <- subset(BE_data,BE_data$Year==year[indexY])
#   }
#   cor_data <- data[ , covari]
#   models <- interac_models_1(limit,respon,covari, covariates_lin,covariates_spl,num.knots,cor_data,interac_type="te")
#   save(list = c("models"),file = paste0("outputs/model_sets/models_slp_lin",year[indexY],".RData"))
# }

