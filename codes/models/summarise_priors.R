# ---------------------------------------------------------------------------------------------
## Authors: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  Summarise priors: Results for Supplement A, Table S4
## ---------------------------------------------------------------------------------------------

year <- c(2010,2011,2013,"yre")

for (y in year){
  v_ut <- paste0("outputs/priors/GAMM/GAMM_prior_es_FishNet50_",y)
  v_fileName <- paste(v_ut,"RData",sep=".")
  load(v_fileName)
  print(paste("Prior:",y))
  print("Mean:")
  print(out_Lap[[1]]$par)
  print("s.d")
  print(diag(solve(-out_Lap[[1]]$hessian)))
}


