##########################################  HPC-INDEX  ##########################################

indexL <- 1
indexY <- as.integer(Sys.getenv("indexY"))
indexO <- as.integer(Sys.getenv("indexO"))
iter <- as.integer(Sys.getenv("iter"))
indexTr <- as.integer(Sys.getenv("indexTr"))
indexI <- as.integer(Sys.getenv("indexI"))

##########################################  CODE  ################################################

year_all <- c(2010,2011,2013) 
year <- year_all[indexY]

tlength_all <- c(500,1000,1500)
tlength <- tlength_all[indexL]

fol<-paste0("designs/polynomial/",year,"/",indexO,"/",indexI,"/")

iname <- paste0("outputs/",fol,"Iteration_",iter-1,".RData")
load(iname)
util_itO <- curr_u

# optimal answer of last transect is equals to the initial input of next iteration
load(paste0("outputs/",fol,"Iteration",iter,"_transect_opt_",indexTr,".RData"))
iname <- paste0("outputs/",fol,"Iteration_",iter,".RData")
save(list = c("tr_design","curr_u"), file = iname)
util_itN <- curr_u

if((util_itN-util_itO)==0){
  save(list = c("tr_design","curr_u"), file = (paste0("outputs/",fol,"Optimal.RData")))
}
