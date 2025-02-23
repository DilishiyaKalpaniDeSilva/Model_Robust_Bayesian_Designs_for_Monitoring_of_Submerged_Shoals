##########################################  HPC-INDEX  ##########################################

indexL <- 1
indexY <- as.integer(Sys.getenv("indexY"))
iter <- as.integer(Sys.getenv("iter"))
indexTr <- as.integer(Sys.getenv("indexTr"))
indexP <- as.integer(Sys.getenv("indexP"))
B <- as.integer(Sys.getenv("B"))            # no of MC draws for the utility approximation
B1 <- as.integer(Sys.getenv("B1"))          # no of MC draws in Laplace approximation 
ncl <- as.integer(Sys.getenv("ncl"))
indexI <- as.integer(Sys.getenv("indexI"))

##########################################  CODE  ################################################

year_all <- c(2010,2011,2013) 
year <- year_all[indexY]
print(year)

tlength_all <- c(500,1000,1500)
tlength <- tlength_all[indexL]
print(tlength)

ntrans_all <- c(18,9,6)
ntrans <- ntrans_all[indexL]

fol<-paste0("designs/GAMM/",year,"/",indexI,"/")

Q1 <- 20
index <- 1

rm <- vector()
for(indexA in 1:Q1){
  fnameA <- paste0(paste0("outputs/",fol,"Iteration_",iter,"/Transect_",indexTr,"/point",indexP,"_angle",indexA),".RData")
  if(!file.exists(fnameA)){
    rm[index] <- indexA
    index <- index + 1
    print(indexA)
  }
}

fnameP <- paste0(paste0("outputs/",fol,"Iteration_",iter,"/Transect_",indexTr,"/point_",indexP,"_rest"),".RData")
save(list = c("rm"), file = fnameP)
