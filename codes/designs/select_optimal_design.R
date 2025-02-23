## ---------------------------------------------------------------------------------------------
## Authors: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  Save the best optimal design out of all the optimal designs under different 
##           initial designs
##           
## ---------------------------------------------------------------------------------------------

##########################################  CODE  ################################################

year_all <- c(2010,2011,2013, "yre") 

## save GAMM best optimal design
for(year in year_all){
  loc_u <- vector()
  for (indexI in 1:5) {
    fol<-paste0("designs/GAMM/",year,"/",indexI,"/")
    fname <- paste0("outputs/",fol,"Optimal.RData")
    if(file.exists(fname)){
      load(fname)
      loc_u[indexI] <- curr_u
    }
  }
  print(loc_u)
  ind <- which.max(loc_u) 
  print(loc_u[ind[1]])
  file.copy(from = paste0("outputs/",paste0("designs/GAMM/",year,"/",ind[1],"/"),"Optimal.RData"), 
            to = paste0("outputs/",paste0("designs/GAMM/",year,"/"),"Optimal.RData"), 
            overwrite = TRUE)
}

## save polynomial best optimal design
for(year in year_all){
  for(indexO in 1:3){
    loc_u <- vector()
    for (indexI in 1:5) {
      fol <- paste0("designs/polynomial/",year,"/",indexO,"/",indexI,"/")
      fname <- paste0("outputs/",fol,"Optimal.RData")
      if(file.exists(fname)){
        load(fname)
        loc_u[indexI] <- curr_u
      }
    }
    print(loc_u)
    ind <- which.max(loc_u) 
    print(loc_u[ind[1]])
    file.copy(from = paste0("outputs/",paste0("designs/polynomial/",year,"/",indexO,"/",ind[1],"/"),"Optimal.RData"), 
              to = paste0("outputs/",paste0("designs/polynomial/",year,"/",indexO,"/"),"Optimal.RData"), 
              overwrite = TRUE)
  }
}
  
