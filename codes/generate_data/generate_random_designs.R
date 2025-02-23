## ---------------------------------------------------------------------------------------------
## Authors: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose: 1. Generate random designs for utility and time comparisons
##          2. Plot generated random designs
## ---------------------------------------------------------------------------------------------

################################ CODE ##########################################################

library(NCmisc)
## check the packages used in the R script
list.functions.in.file("codes/generate_data/generate_random_designs.R", alphabetic = TRUE) 

library(DSpat)
library(raster)
library(sf)
library(RColorBrewer)
library(spatstat)

gen_random_designs <- function(nd, md, n, numt, files, filep){
  load(filep)
  ## remove empty transects which are outside the shoal
  j <- 1; ind <- vector()
  for(i in 1:length(coordp)){
    if(length(coordp[[i]])==1){
      ind[j] <- i
      j <- j+1
    }
  }
  
  nt <- (1:length(coordp))[-ind]
  
  fll <-list()
  fl <- vector()
  for (ind_fish_size in 1:5) {
    for(i in nt){
      fl[i] <- (length(unique(fishnet_all[[ind_fish_size]][[i]])))
    } 
    fll[[ind_fish_size]] <- fl 
    print(table(fl))
  }
  
  set.seed(55)
  design <- list()
  j <-1
  shallow_t <- vector()
  for(i in nt){
    if(!any(depthp[[i]]<(md))){
      if(length(depthp[[i]])==n){
        shallow_t[j] <- i
        j <- j+1
      }
    }
  }
  
  for(i in 1:nd){
    design[[i]] <- sample(shallow_t, size =  numt,replace = FALSE)
  }
  
  j <-1
  deep_t <- vector()
  for(i in nt){
    if(!any(depthp[[i]]>(md))){
      if(length(depthp[[i]])==n){
        deep_t[j] <- i
        j <- j+1
      }
    }
  }
  
  for(i in (nd+1):(nd*2)){
    design[[i]] <- sample(deep_t,size =  numt,replace = FALSE)
  }
  
  middle_t <- vector()
  middle_tt <- c(nt)[-c(shallow_t,deep_t)]
  j <- 1
  for(i in middle_tt){
    if(length(depthp[[i]])==n){
      middle_t[j] <- i
      j <- j+1
    }
  }
  
  for(i in (2*nd+1):(3*nd)){
    design[[i]] <- sample(middle_t,size =  numt,replace = FALSE)
  }
  
  clst <- nt[seq(1,length(nt),2)]
  clust_t <- vector()
  j <-1
  for(i in clst){
    if(length(depthp[[i]])==n){
      clust_t[j] <- i
      j <- j+1
    }
  }
  
  j <- 1
  for(i in (3*nd+1):(4*nd)){
    design[[i]] <- clust_t[j:(j+(numt-1))]
    j <- j+numt
  }
  
  save(list = c("design"), file = files)
  return(design)
}

nd <-25; md <- (-30); n <-50; numt<- 18

files <- "data/transect_data/random_designs_fs.RData"
filep <- paste0("data/transect_data/data_for_design_fshsize_predang100_n50_dep50_len500.RData")

gen_random_designs(nd, md, n, numt, files, filep)
  
## plot random designs
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- rep(col_vector,50)

files <- "data/transect_data/random_designs_fs.RData"
filep <- paste0("data/transect_data/data_for_design_fshsize_predang100_n50_dep50_len500.RData")

load(filep)
load("data/BE_raster_data/BE_raster.RData") 
load(file = files)
source("codes/functions/plot_design.R")

fll <-list()
fl <- vector()
for (ind_fish_size in 1:5) {
  for(i in 1:length(design)){
    fl[i] <- (length(unique(unlist(fishnet_all[[ind_fish_size]][design[[i]]]))))
  } 
  fll[[ind_fish_size]] <- fl 
  print(table(fl))
}

fl <- list()
for(i in 1:length(design)){
  fl[[i]] <- range(unlist(depthp[design[[i]]]))
} 
fll[[ind_fish_size]] <- fl 
print(table(unlist(depthp[design[[i]]])))

pdf(file = "plots/generate_data/random_designs.pdf")   
for(j in 1:length(design)){
  tr <- sort(design[[j]])
  
  plot(covP3,  xlim=c(610000,614200),main=j, xlab="Easting", ylab="Northing",legend.only=FALSE,
       legend.args=list(text='Depth (m)',side=4,font=2, line=2.5, cex=0.8))
  abline(v=xd)
  abline(h=yd)
  
  ntt <- length(lnsp[[1]]$transects)
  point <- ceiling(tr/ntt)
  print(point)
  rem <- tr - (point-1)*ntt
  
  for(i in 1:length(tr)){
    plot(owin(poly=lnsp[[point[i]]]$transects[[rem[i]]]),col=col_vector[i],add=TRUE)
  }
}

dev.off()

