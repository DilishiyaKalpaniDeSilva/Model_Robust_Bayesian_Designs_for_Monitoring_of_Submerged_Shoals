## ---------------------------------------------------------------------------------------------
## Authors: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  1. Generate fishnet ids for different fishnet sizes
##           2. Get subsets of Barracouta_East_Rdata to use as historical data for designs
## ---------------------------------------------------------------------------------------------

################################ CODE ###########################################################

library(NCmisc)
## check the packages used in the R script
list.functions.in.file("codes/generate_data/generate_covariate_data.R", alphabetic = TRUE) 

library(readxl)
library(dplyr)
library(xlsx)

load("data/BE_raster_data/fishnets_diff_sizes.RData")
source("codes/functions/find_fishnet.R")

# Download Barracouta East R data from https://doi.org/10.25845/p1sr-s997 and save it as Barracouta_East_Rdata.xlsx
Barracouta_East_Rdata <- data.frame(read_excel("data/BE_data/Barracouta_East_Rdata.xlsx"))
unique(Barracouta_East_Rdata$Year)
dim(Barracouta_East_Rdata)
table(Barracouta_East_Rdata$Year)
range(Barracouta_East_Rdata$depth) 
aggregate(Barracouta_East_Rdata$depth ~ Barracouta_East_Rdata$Year, 
                          data = Barracouta_East_Rdata, function(x) c(min = min(x), max = max(x)))

## Generate fishnet ids for different fishnet sizes
name <- paste0("FishNet",seq(50,250,50))
for(i in 1:5){
  fishnet <- find_fishnet(fishnet_data[[i]],cbind(Barracouta_East_Rdata$Easting,Barracouta_East_Rdata$Northing))
  Barracouta_East_Rdata[,(20+i)] <- fishnet
}

colnames(Barracouta_East_Rdata)[21:25] <- name
Barracouta_East_Rdata_FS <- Barracouta_East_Rdata

write.xlsx(x=Barracouta_East_Rdata_FS, file="data/BE_data/Barracouta_East_Rdata_FS.xlsx", row.names = FALSE)

## Get subsets of Barracouta_East_Rdata to use as historical data for designs
BE_data_2010 <- filter(Barracouta_East_Rdata_FS, Year==2010)
dim(BE_data_2010)
BE_data_2010 <- BE_data_2010[,c("Year", "HardCoral", "depth", name)]

BE_data_2011 <- filter(Barracouta_East_Rdata_FS, Year==2011)
dim(BE_data_2011)
BE_data_2011 <- BE_data_2011[,c("Year","HardCoral", "depth", name)]

BE_data_2013 <- filter(Barracouta_East_Rdata_FS, Year==2013)
dim(BE_data_2013)
BE_data_2013 <- BE_data_2013[,c("Year","HardCoral", "depth", name)]

BE_data_yre <- rbind(BE_data_2010,BE_data_2011,BE_data_2013)
dim(BE_data_yre)
table(BE_data_yre$Year)

range(BE_data_2010$depth)
range(BE_data_2011$depth)
range(BE_data_2013$depth)
range(BE_data_yre$depth)

write.xlsx(x=BE_data_2010, file="data/BE_data/BE_data_2010.xlsx", row.names = FALSE)
write.xlsx(x=BE_data_2011, file="data/BE_data/BE_data_2011.xlsx", row.names = FALSE)
write.xlsx(x=BE_data_2013, file="data/BE_data/BE_data_2013.xlsx", row.names = FALSE)
write.xlsx(x=BE_data_yre, file="data/BE_data/BE_data_yre.xlsx", row.names = FALSE)
