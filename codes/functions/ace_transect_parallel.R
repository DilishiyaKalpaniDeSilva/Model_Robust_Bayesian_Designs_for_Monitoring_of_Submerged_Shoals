library(DiceKriging)

f.pred<-function(x,fit)  #### Change how you predict from GP
{
  nd<-matrix(x,nrow=1)  
  #Predicted values and (marginal of joint) conditional variances based on a km model. 95 % confidence intervals are given, based on strong assumptions: Gaussian process assumption, specific prior distribution on the trend parameters, known covariance parameters. This might be abusive in particular in the case where estimated covariance parameters are plugged in
  ans<-predict(fit,newdata=data.frame(z=nd),type="UK",cov.compute=TRUE,se.compute=TRUE)
  out <- ans$mean
  
  return(out)
}

ace_transect <- function(utility,point,B,Q1,Q2,lower,upper, by){
  
  z <- seq(lower,upper,by=(upper/Q1))[-(Q1+1)] # angles to find the utility, 20 fixed angles, similar angles as generated data
  uz_all <- list()
  uz <- vector() #  to save mean utility for 20 transects
  zp <- seq(lower,upper,by=(upper/Q2))[-(Q2+1)]  ## new angles to predict for the GP emulator
  up <- vector() # to save predicted utility values from zp
  
  # calculating utility when each of 20 transects are replaced with the t transect of the design
  print("tr design:")
  for(indexA in 1:Q1){
    fnameA <- paste0(paste0("outputs/",fol,"/Iteration_",iter,"/Transect_",indexTr,"/point",indexP,"_angle",indexA),".RData")
    if(file.exists(fnameA)){
      load(fnameA)
      uz_all[[indexA]] <- uz_allA
      uz[indexA] <- uzA
      print(tr_design)
      print(indexA)
    }else{
      tr_no <- (5*(Q1*(point-1) + indexA))-4 # from 1 to 300
      tr_design[indexTr] <- tr_no
      print(tr_design)
      
      cl <- makePSOCKcluster(ncl)
      clusterCall(cl, RcppFunctions)
      registerDoParallel(cl)
      
      uz_allA <- uz_all[[indexA]] <- utility(d=tr_design,B)
      uzA <- uz[indexA] <- mean(uz_all[[indexA]],na.rm=TRUE)
      save(list = c("tr_design","uz_allA","uzA"), file = fnameA)
      
      stopCluster(cl)
    }
  }
  
  print(uz)
  ind_c <- which(uz==max(uz,na.rm = TRUE))#index of maximum utility
  tr_no <- (5*(Q1*(point-1) + ind_c[1]))-4 
  curr_x <- tr_no #x value which maximizes the utility, for each x value
  curr_u <- mean(uz_all[[ind_c[1]]],na.rm=TRUE)
  
  tryCatch({
    if(!any(is.na(uz))){
      fit <- km(formula=~1,design=data.frame(z=z),response=uz,covtype="gauss",control=list(trace=FALSE),nugget.estim=TRUE)
      print("no NA values")
    }else{
      fit <- km(formula=~1,design=data.frame(z=z[-which(is.na(uz))]),response=uz[!is.na(uz)],covtype="gauss",control=list(trace=FALSE),nugget.estim=TRUE)
      print("contains NA values")
    }
  }, error = function(e) { 
    out <- list(curr_x,curr_u)
    return(out)})
  
  #### Emulation, #distribution of utility
  for(k in 1:Q2){
    tr_no1 <- (Q2*(point-1)) + k
    up[k] <- f.pred(zp[k],fit) # predicting for new angles in zp
  }
  
  ind <- which(up==max(up,na.rm = TRUE))#index of maximum utility
  tr_no1 <- (Q2*(point-1)) + ind[1]
  tr_design[indexTr] <- tr_no1
  prop_x <- tr_no1 #x value which maximizes the utility, for each x value

  cl <- makePSOCKcluster(ncl)
  clusterCall(cl, RcppFunctions)
  registerDoParallel(cl)
  prop_u <- mean(utility(tr_design,B),na.rm=TRUE)
  stopCluster(cl)
  
  if(prop_u >= curr_u){ #compare with cu
    curr_x <- prop_x
    curr_u <- prop_u
  }
  
  out <- list(curr_t=curr_x,curr_u=curr_u)
  return(out)
}

