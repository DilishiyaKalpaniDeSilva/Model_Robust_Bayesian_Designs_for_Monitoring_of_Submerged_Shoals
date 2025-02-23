## ---------------------------------------------------------------------------------------------
## Authors: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose: Check correlation between polynomial terms for different scaling of data
## ---------------------------------------------------------------------------------------------

################################ CODE ###########################################################

library(corrplot)
library(readxl)

source("codes/functions/basic_functions.R")

BE_data <- read_excel("data/BE_data/Barracouta_East_Rdata_FS.xlsx")
year <- c(2010,2011,2013,"yre")

pdf(file = "plots/generate_data/correlation_plots.pdf")
par(mfrow=c(2,2))

for(indexY in 1:length(year)){
  if(year[indexY]=="yre"){
    data <- BE_data
  }else{
    data <- subset(BE_data,BE_data$Year==year[indexY])
  }
  print(paste("year:",year[indexY]))
  data$depth <- (data$depth)
  print(range(data$depth))
  Xd_hist <- outer(X = data$depth,Y = (1:5),FUN = "^")
  
  M <- cor(Xd_hist)
  corrplot(M, method = 'number',mar=c(0,0,1,0),
           type = 'lower',title = paste("Original data","year:",year[indexY]))
}

for(indexY in 1:length(year)){
  if(year[indexY]=="yre"){
    data <- BE_data
  }else{
    data <- subset(BE_data,BE_data$Year==year[indexY])
  }
  print(paste("year:",year[indexY]))
  data$depth <- norm_01(data$depth)
  print(range(data$depth))
  Xd_hist <- outer(X = data$depth,Y = (1:5),FUN = "^")
  
  M <- cor(Xd_hist)
  corrplot(M, method = 'number',mar=c(0,0,1,0),
                 type = 'lower',title = paste("range 0 1 ","year:",year[indexY]))
}

for(indexY in 1:length(year)){
  if(year[indexY]=="yre"){
    data <- BE_data
  }else{
    data <- subset(BE_data,BE_data$Year==year[indexY])
  }
  print(paste("year:",year[indexY]))
  data$depth <- scaleP(data$depth,min(data$depth),max(data$depth),-1,1)
  print(range(data$depth))
  Xd_hist <- outer(X = data$depth,Y = (1:5),FUN = "^")
  
  M <- cor(Xd_hist)
  corrplot(M, method = 'number',mar=c(0,0,1,0),
           type = 'lower',title = paste("range -1 1 ","year:",year[indexY]))
}

for(indexY in 1:length(year)){
  if(year[indexY]=="yre"){
    data <- BE_data
  }else{
    data <- subset(BE_data,BE_data$Year==year[indexY])
  }
  print(paste("year:",year[indexY]))
  data$depth <- data$depth - mean(data$depth)
  print(range(data$depth))
  Xd_hist <- outer(X = data$depth,Y = (1:5),FUN = "^")
   
  M <- cor(Xd_hist)
  corrplot(M, method = 'number',mar=c(0,0,1,0),
           type = 'lower',title = paste("depth-mean(depth) ","year:",year[indexY]))
}

for(indexY in 1:length(year)){
  if(year[indexY]=="yre"){
    data <- BE_data
  }else{
    data <- subset(BE_data,BE_data$Year==year[indexY])
  }
  print(paste("year:",year[indexY]))
  data$depth <- stand(data$depth)
  print(range(data$depth))
  Xd_hist <- outer(X = data$depth,Y = (1:5),FUN = "^")
   
  M <- cor(Xd_hist)
  corrplot(M, method = 'number',mar=c(0,0,1,0),
           type = 'lower',title = paste("standadize ","year:",year[indexY]))
}

dev.off()
