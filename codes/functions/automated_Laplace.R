## ---------------------------------------------------------------------------------------------
## Authors: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  Create an automated Laplace approximation to approximate posterior for generalised additive 
##           mixed model class 
## ---------------------------------------------------------------------------------------------

# library(NCmisc)
# list.functions.in.file("codes/functions/automated_Laplace.R", alphabetic = TRUE)

library(dplyr)
library(formula.tools)
library(plyr)
library(tidyverse)
library(VGAM)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(RcppGSL)
library(RcppParallel)

sourceCpp("codes/functions/matrix_calc.cpp")
sourceCpp("codes/functions/cpp_functions.cpp")

## to approximate posterior of theta 
log_posterior <- function(param,Y,X,rand_u ,z.spline,rand_U,rand_Y,B,random,Year,prior,
                          ms,sds,mr,sdr,mt,sdt){
  
  ############# spline ###############
  if(!is.null(z.spline)){
    if(prior=="loggamma"){
      log_prec_s <- param["log_prec_s"]
      ps <- dlgamma(log_prec_s,location=-33,
                    scale = 5,shape = 10,log = TRUE)
      u <- qnorm(rand_u,0, 1/(sqrt(exp(log_prec_s))))
    }else{
      log_sigma_s <- param["log_sigma_s"]
      ps <- dnorm(log_sigma_s, ms, sds,log = TRUE) 
      u <- qnorm(rand_u,0, exp(log_sigma_s))
    }
    spl <- eigenMapMatMult2(z.spline,u, n_cores = 8)
  }else{
    spl <- 0; ps <- 0
    
  }
  
  ############# fishnet random effect ###############
  if(!is.null(random)){
    if(prior=="loggamma"){
      log_prec_r <- param["log_prec_r"]
      pr <- dlgamma(log_prec_r,location=-21,
                    scale = 5, shape = 10,log = TRUE)  
      U <- qnorm(rand_U,0, 1/(sqrt(exp(log_prec_r))))
    }else{
      log_sigma_r <- param["log_sigma_r"]
      pr <- dnorm(log_sigma_r, mr, sdr, log = TRUE)
      U <- qnorm(rand_U,0, exp(log_sigma_r))
    } 
    fish_re <- U[random,]
  }else{
    fish_re <- 0; pr <- 0
  }
  
  ############# Year random effect ###############
  if(!is.null(Year)){
    if(prior=="loggamma"){
      log_prec_t <- param["log_prec_t"]
      pt <- dlgamma(log_prec_t,location=-21, 
                    scale = 5, shape = 10,log = TRUE)  
      Yr <- qnorm(rand_Y,0, 1/(sqrt(exp(log_prec_t))))
    }else{
      log_sigma_t <- param["log_sigma_t"]
      pt <- dnorm(log_sigma_t,  mt, sdt ,log = TRUE)
      Yr <- qnorm(rand_Y,0, exp(log_sigma_t))
    } 
    Yr1 <- rbind(0,Yr)
    Year_re <- Yr1[Year,]
  }else{
    Year_re <- 0; pt <- 0
  }
  
  x_beta <- as.vector(X %*% param[1:ncol(X)] )
  sum_expit <- x_beta + spl + fish_re + Year_re
  prob <-  expitCpp(as.matrix(sum_expit)) 
  like <- colSums( (dbinomCpp( Y, prob) ) )
  log_like <- mean(like)
  
  pb <- sum(dnorm(x=param[1:ncol(X)],
                  mean=rep(0,ncol(X)),sd=rep(5,ncol(X)),log=TRUE))
  log_prior <- pb + ps + pr + pt
  
  log_post<- log_like+log_prior
  if(is.infinite(log_post)||is.na(log_like)){
    #print(-1000000)
    return(-1000000)
  }else{
    #print(paste("ni:",log_post))
    return(log_post)
  }
}


log_likelihood <- function(param1,Y,X,z.spline,random,mu_prior,num.knots,m, Year,prior){
  
  ############# spline  ###############
  if(!is.null(z.spline)){
    u <- param1[1:sum(num.knots)]
    
    if(prior=="loggamma"){
      sigma_s <- 1/(sqrt(exp(mu_prior["log_prec_s"])))
    }else{
      sigma_s <- exp(mu_prior["log_sigma_s"])
    }
    
    prior_u <- sum(dnorm(x=u, mean=0,sd=sigma_s,log=TRUE))
    s <- eigenMapMatMult2(z.spline , u, n_cores = 4)
    
  }else{
    s <- 0; prior_u <- 0
    num.knots <- 0
  }
  
  ############# fishnet random effect ###############
  if(!is.null(random)){
    U <- param1[(sum(num.knots)+1):(sum(num.knots)+m)]
    if(prior=="loggamma"){
      sigma_r <- 1/(sqrt(exp(mu_prior["log_prec_r"])))
    }else{
      sigma_r <- exp(mu_prior["log_sigma_r"])
    }
    f_re <- U[random]
    prior_r <- sum(dnorm(x=U, mean=0,sd=sigma_r,log=TRUE))
  }else{
    f_re <- 0; prior_r <- 0
  }
  
  ############# Year random effect ###############
  if(!is.null(Year)){
    Yr <- tail(param1,(length(levels(Year))-1))
    Yr1 <- c(0,Yr)
    if(prior=="loggamma"){
      sigma_t <- 1/(sqrt(exp(mu_prior["log_prec_t"])))
    }else{
      sigma_t <- exp(mu_prior["log_sigma_t"])
    }
    y_re <- Yr1[Year]
    prior_t <- sum(dnorm(x=Yr, mean=0,sd=sigma_t,log=TRUE))
  }else{
    y_re <- 0; prior_t <- 0
  }
  
  #probability of success for each data point
  prob  <- expit(X%*%mu_prior[1:ncol(X)] + s + f_re + y_re)
  
  #likelihood of the data
  log_like <- sum(dbinom(x=Y,size = 20,prob =  prob ,log=TRUE)) + prior_u + prior_r + prior_t
  
  if(!is.infinite(log_like)|| !is.na(log_like)){
    #print(log_like)
    return(log_like)    
  }else{
    #print(-10000000)
    return(-10000000)
  }
  
}

## automated function for Laplace approximation
Laplace_approx_gamm_OSullivan <- function(f,data,B,fre,yre,prior,alpha,ms,sds,mr,sdr,mt,sdt,HesP,HesR,
                                          opt_method,beta_start,knot_type,intKnots,boundKnots, basis ){
  suppressMessages(attach(data))
  f <- as.formula(f)
  Y <- data %>% pull(lhs.vars(f))
  var_names <- all.vars(f)[-1]
  terms <- rhs.vars(f)
  poly_terms <- var_names[which(grepl( "poly(",terms,  fixed = TRUE)==TRUE)]
  deg <- readr::parse_number(terms[which(grepl( "poly(",terms,  fixed = TRUE)==TRUE)])
  spline_terms <- var_names[which(grepl( "s(",terms,  fixed = TRUE)==TRUE)]
  ind <- c(which(grepl( "poly(",terms,  fixed = TRUE)==TRUE),
           which(grepl( "s(",terms,  fixed = TRUE)==TRUE))
  linear_terms <- var_names[!var_names %in% c(spline_terms,poly_terms)]
  
  if(length(poly_terms)!=0){
    X1 <- outer(X = data %>% dplyr::pull(poly_terms),Y = 1:deg,FUN = "^")
    X <- as.matrix(cbind(1,data %>% dplyr::select(c(linear_terms,spline_terms)),X1))
  }else{
    X <- as.matrix(cbind(1,data %>% dplyr::select(c(linear_terms,spline_terms))))
  }
  set.seed(5)
  param <- c("beta"=rep(beta_start,ncol(X)))
  param1 <- NULL
  
  if(length(spline_terms)!=0){  ## spline terms
    XsplPreds <- (as.matrix(data %>% dplyr::select(spline_terms)))
    # XsplPreds <- apply(XsplPreds, 2, norm_01)
    num.knots <- as.numeric(gsub(')','',
                                 sub(".*= ", "", 
                                     terms[which(grepl( "s(",terms,  fixed = TRUE)==TRUE)])))
    if(any(is.na(num.knots))){
      num.knots <- nknots
    }
    #k <- sum(num.knots)
    zz <- cal_z(numBasis = num.knots,XsplPreds,knot_type ,intKnots,boundKnots, basis )
    z.spline <- zz$Z
    knot_loc <- unlist(zz$knot_loc)
    rand_u <- matrix(runif(B*sum(num.knots)),sum(num.knots),B)
    
    if(prior=="loggamma"){
      param <- c(param,"log_prec_s"=-1)
    }else{
      param <- c(param,"log_sigma_s"=ms+0.1)
    }
    param1 <- c(param1,rep(0,sum(num.knots)))
  }else{
    rand_u <- NULL; z.spline <- NULL
    knot_loc <- NULL
  }
  
  if(fre){ ## fishnet random effect
    random <- factor(data$FishNet)
    m <- length(levels(random))
    random <- plyr::mapvalues(random, from =levels(random), to=c(1:m))
    rand_U <- matrix(runif(B*m),m,B) # generating probabilities
    
    if(prior=="loggamma"){
      param <- c(param,"log_prec_r"=-1)
    }else{
      param <- c(param,"log_sigma_r"=mr-0.1)
    }
    param1 <- c(param1,rep(0,m))
  }else{
    rand_U <- NULL; random <- NULL
  }
  
  if(yre){ ## Year random effect
    Year <- factor(data$Year)
    Year <- plyr::mapvalues(Year, from =levels(Year), to=c(1:length(levels(Year))))
    rand_Y <- matrix(runif(B*(length(levels(Year))-1)),(length(levels(Year))-1),B)
    
    if(prior=="loggamma"){
      param <- c(param, "log_prec_t"=-1)
    }else{
      param <- c(param, "log_sigma_t"=mt-0.1)
    }
    param1 <- c(param1,rep(0,(length(levels(Year))-1)))
  }else{
    rand_Y <- NULL; Year <- NULL
  }
  
  set.seed(5)
  t1 <- Sys.time()
  fitOptim <- optim(par=param, fn=log_posterior, method = opt_method,
                    control = list(fnscale = -1,maxit=100000,factr=1e8,
                                   trace=1,REPORT=1),
                    hessian = HesP, X=X, Y=Y, rand_u=rand_u, rand_U=rand_U, 
                    rand_Y=rand_Y ,z.spline=z.spline, B=B,random=random, Year=Year,
                    prior=prior,ms=ms,sds=sds,mr=mr,sdr=sdr,mt=mt,sdt=sdt)
  t2 <- Sys.time()
  time1 <- difftime(t2, t1, units = "secs")
  print(time1)
  print(fitOptim$convergence)
  
  if(!is.null(param1) && alpha){
    mu_prior <- fitOptim$par
    t1 <- Sys.time()
    fitOptim_u <- optim(par=param1, fn=log_likelihood, method = "L-BFGS-B",
                        control = list(fnscale = -1,maxit=10000000,factr=1e11,
                                       trace=1,REPORT=1), 
                        hessian = HesR,X=X,Y=Y,z.spline=z.spline,num.knots=num.knots,
                        m=m,random=random,mu_prior=mu_prior,Year=Year,prior=prior)
    t2 <- Sys.time()
    time2 <- difftime(t2, t1, units = "secs")
    print(fitOptim_u$convergence)
    print(time2)
    out <- list(fitOptim,fitOptim_u,knot_loc,time1,time2)
  }else{
    out <- list(fitOptim,knot_loc,time1)
  }
  
  return(out)
  
}

