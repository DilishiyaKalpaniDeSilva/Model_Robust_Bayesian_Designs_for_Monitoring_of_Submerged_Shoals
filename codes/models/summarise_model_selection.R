## ---------------------------------------------------------------------------------------------
## Author: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  Summarise model selection
## ---------------------------------------------------------------------------------------------

library(xlsx)

year <- c(2010,2011,2013,"yre")
app <- c(FALSE,rep(TRUE,3))

log_model_evi <- vector()

for(indexY in 1:length(year)){
  load(paste0("outputs/model_sets/models_",year[indexY],".RData"))
  for(indexM in 1:nrow(models)){
    v_ut <- paste0("outputs/model_selection/",year[indexY],"/model_lg_",indexM)
    v_fileName <- paste(v_ut,"RData",sep=".")
    if(file.exists(v_fileName)){
      load(v_fileName)
      log_model_evi[indexM] <- model_evi
    }else{
      print(indexM)
      log_model_evi[indexM] <- NA
    }
  }
  out <- data.frame(models,log_model_evi)
  outsort <- out[order(log_model_evi, decreasing = TRUE),]
  
  write.xlsx(outsort, file="outputs/model_selection/model_selection_summary.xlsx", sheetName = as.character(year[indexY]),
             append = app[indexY], row.names=FALSE)
  
}


