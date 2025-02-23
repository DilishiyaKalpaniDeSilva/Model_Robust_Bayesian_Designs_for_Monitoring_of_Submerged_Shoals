## ---------------------------------------------------------------------------------------------
## Authors: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  Create an automated function to evaluate posterior for generalised additive 
##           mixed model class (with tensor product)
## ---------------------------------------------------------------------------------------------

library(NCmisc)
list.functions.in.file("codes/functions/automated_Laplace.R", alphabetic = TRUE)

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

source("codes/functions/cal_z.R")
sourceCpp("codes/functions/matrix_calc.cpp")
sourceCpp("codes/functions/cpp_functions.cpp")

## to approximate posterior of theta 
log_posterior <- function(param,X, Y, 
                          rand_u, num.knots,z.spline,
                          rand_te,x1x2,num.knotst,
                          rand_U, rand_Y , B,random, Year,
                          prior){
  
  ############# spline ###############
  if(!is.null(num.knots)){
    if(prior=="loggamma"){
      log_prec_s <- param[paste0("log_prec_s",1:length(num.knots))]
      ps <- dlgamma(log_prec_s,location=-33,
                    scale = 5,shape = 10,log = TRUE)
      inds <- c(0,cumsum(num.knots))
      u <- NULL
      for(i in 1:length(num.knots)){
        u <- rbind(u,qnorm(rand_u[(inds[i]+1):inds[i+1],],0, 1/(sqrt(exp(log_prec_s[i])))))
      }
    }
    
    spl <- eigenMapMatMult2(z.spline,u, n_cores = 4)
  }else{
    spl <- 0; ps <- 0
    
  }
  
  ############# tensor ###############
  if(!is.null(num.knotst)){
    if(prior=="loggamma"){
      log_prec_te <- param[paste0("log_prec_te",1:length(num.knotst))]
      pte <- dlgamma(log_prec_te,location=-21,
                    scale = 5,shape = 10,log = TRUE)
      indte <- c(0,cumsum(num.knotst^2))
      te <- NULL
      for(i in 1:length(num.knotst)){
        te <- rbind(te,qnorm(rand_te[(indte[i]+1):indte[i+1],],0, 1/(sqrt(exp(log_prec_te[i])))))
      }
    }
    
    tens <- eigenMapMatMult2(x1x2,te, n_cores = 4)
  }else{
    tens <- 0; pte <- 0
    
  }
  
  ############# fishnet random effect ###############
  if(!is.null(random)){
    if(prior=="loggamma"){
      log_prec_r <- param["log_prec_r"]
      pr <- dlgamma(log_prec_r,location=-21,
                    scale = 5, shape = 10,log = TRUE)  
      U <- qnorm(rand_U,0, 1/(sqrt(exp(log_prec_r))))
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
  sum_expit <- x_beta + spl + tens + fish_re + Year_re
  prob <-  expitCpp(as.matrix(sum_expit)) 
  like <- colSums( (dbinomCpp( Y, prob) ) )
  log_like <- mean(like)
  
  pb <- sum(dnorm(x=param[1:ncol(X)],
                  mean=rep(0,ncol(X)),sd=rep(5,ncol(X)),log=TRUE))
  log_prior <- pb + sum(ps) + sum(pte) + pr + pt
  
  log_post <- log_like + log_prior
  if(is.infinite(log_post)||is.na(log_like)){
    return(-10000000)
  }else{
    return(log_post)
  }
}


log_likelihood <- function(param1,Y,X,z.spline,x1x2,num.knotst,random,mu_prior,num.knots,m, Year,prior){
  
  ############# spline  ###############
  if(!is.null(num.knots)){
    u <- param1[paste0("u",1:sum(num.knots))]
    
    if(prior=="loggamma"){
      sigma_s <- 1/(sqrt(exp(mu_prior[paste0("log_prec_s",1:length(num.knots))])))
    }
    
    inds <- c(0,cumsum(num.knots))
    prior_u <- sapply(1:length(num.knots),function(i){
      sum(dnorm(x=u[(inds[i]+1):inds[i+1]], mean=0,sd=sigma_s[i],log=TRUE))
    })
    s <- eigenMapMatMult2(z.spline , u, n_cores = 4)
    
  }else{
    s <- 0; prior_u <- 0
  }
  
  ############# tensor  ###############
  if(!is.null(num.knotst)){
    te <- param1[paste0("te",1:ncol(x1x2))]
    
    if(prior=="loggamma"){
      sigma_te <- 1/(sqrt(exp(mu_prior[paste0("log_prec_te",1:length(num.knotst))])))
    }
    
    indte <- c(0,cumsum(num.knotst^2))
    prior_te <- sapply(1:length(num.knotst),function(i){
      sum(dnorm(x=te[(indte[i]+1):indte[i+1]], mean=0,sd=sigma_te[i],log=TRUE))
    })
    tens <- eigenMapMatMult2(x1x2 , te, n_cores = 4)
    
  }else{
    tens <- 0; prior_te <- 0
  }
  
  
  
  ############# fishnet random effect ###############
  if(!is.null(random)){
    U <- param1[paste0("U",1:m)]
    if(prior=="loggamma"){
      sigma_r <- 1/(sqrt(exp(mu_prior["log_prec_r"])))
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
    }
    y_re <- Yr1[Year]
    prior_t <- sum(dnorm(x=Yr, mean=0,sd=sigma_t,log=TRUE))
  }else{
    y_re <- 0; prior_t <- 0
  }
  
  #probability of success for each data point
  prob  <- expit(X%*%mu_prior[1:ncol(X)] + s + tens + f_re + y_re)
  
  #likelihood of the data
  log_like <- sum(dbinom(x=Y,size = 20,prob =  prob ,log=TRUE)) + sum(prior_u) + sum(prior_te) + prior_r + prior_t
  
  if(!is.infinite(log_like)|| !is.na(log_like)){
    #print(log_like)
    return(log_like)    
  }else{
    #print(-10000000)
    return(-10000000)
  }
  
}

## automated function for Laplace approximation
Laplace_approx_gamm_OSullivan <- function(f,data,B,fre,yre,prior,alpha,
                                          HesP,HesR,opt_method,beta_start){
  suppressMessages(attach(data))
  f <- as.formula(f)
  var_names <- all.vars(f)[-1]
  terms <- rhs.vars(f)
  
  poly_all <- terms[which(grepl( "poly(",terms,  fixed = TRUE)==TRUE)]
  if(length(poly_all)!=0){
    pp <- sapply(base::strsplit(poly_all, ","), "[", 1)
    poly_terms <- sapply(base::strsplit(pp, "\\("), "[", 2)
    deg <- readr::parse_number(terms[which(grepl( "poly(",terms,  fixed = TRUE)==TRUE)])  
  }else{
    poly_terms <- NULL
  }
  
  spline_all <- terms[which(grepl( "s(",terms,  fixed = TRUE)==TRUE)]
  spline_terms <- sapply(base::strsplit(sapply(base::strsplit(spline_all, ","), "[", 1), "\\("), "[", 2)
  
  tensor_all <- terms[which(grepl( "te(",terms,  fixed = TRUE)==TRUE)]
  
  linear_terms <- terms[!(terms %in% c(poly_all,spline_all,tensor_all))]
  
  Y <- data %>% pull(lhs.vars(f))
  if(length(poly_terms)!=0){
    X1 <- outer(X = data %>% dplyr::pull(poly_terms),Y = 1:deg,FUN = "^")
    X <- as.matrix(cbind(1,data %>% dplyr::select(c(linear_terms,spline_terms)),X1))
  }else{
    X <- as.matrix(cbind(1,data %>% dplyr::select(c(linear_terms,spline_terms))))
  }
  set.seed(5)
  param <- c(rep(beta_start,ncol(X))) ## population parameters
  namespb <- paste0("beta",1:ncol(X))
  param1 <- NULL                             ## random effect parameters
  
  if(length(spline_terms)!=0){  ## spline terms
    XsplPreds <- as.matrix(data %>% dplyr::select(spline_terms))
    num.knots <- as.numeric(gsub(')','',
                                 sub(".*= ", "", 
                                     terms[which(grepl( "s(",terms,  fixed = TRUE)==TRUE)])))
    zz <- cal_z(numBasis = num.knots,XsplPreds,knot_type="es" ,intKnots=NA,boundKnots=NA, basis="OS" )
    z.spline <- zz$Z
    rand_u <- matrix(runif(B*sum(num.knots)),sum(num.knots),B)
    
    if(prior=="loggamma"){
      param <- c(param,rep(-1,length(spline_terms)))
      namesps <- paste0("log_prec_s",1:length(spline_terms))
    }
    param1 <- c(param1,"u"=rep(0,sum(num.knots)))
  }else{
    rand_u <- NULL; z.spline <- NULL; num.knots <- NULL; 
    namesps <- NULL
  }
  
  if(length(tensor_all)!=0){
    x1x2 <-  NULL; num.knotst <- vector()
    for(i in 1:length(tensor_all)){
      x1 <- strsplit(strsplit(tensor_all[i],",")[[1]][1], "\\(")[[1]][2]
      x2 <- strsplit(tensor_all[i],",")[[1]][2]
      XsplPreds1 <- as.matrix(data %>% dplyr::select(gsub(" ", "", x1)))
      XsplPreds2 <- as.matrix(data %>% dplyr::select(gsub(" ", "", x2)))
      num.knotst[i] <- as.numeric(gsub(')','',
                                    sub(".*= ", "", 
                                        terms[which(grepl( "te(",terms,  fixed = TRUE)==TRUE)])))[i]
      z.spline1 <- cal_z(numBasis = num.knots,XsplPreds1, knot_type="es" ,intKnots=NA,boundKnots=NA, basis="OS" )$Z
      z.spline2 <- cal_z(numBasis = num.knots,XsplPreds2, knot_type="es" ,intKnots=NA,boundKnots=NA, basis="OS" )$Z
      x1x2 <- cbind(x1x2, z.spline1 %.% z.spline2)
    }
    rand_te <- matrix(runif(B*ncol(x1x2)),ncol(x1x2),B) 
    
    if(prior=="loggamma"){
      param <- c(param,rep(-1,length(tensor_all)))
      namespte <- paste0("log_prec_te",1:length(tensor_all))
    }
    param1 <- c(param1,"te"=rep(0,ncol(x1x2)))
  }else{
    rand_te <- NULL; x1x2 <- NULL; num.knotst <- NULL; 
    namespte <- NULL
  }
  names(param) <- c(namespb,namesps,namespte)
  
  if(fre){ ## fishnet random effect
    random <- factor(data$FishNet)
    m <- length(levels(random))
    random <- plyr::mapvalues(random, from =levels(random), to=c(1:m))
    rand_U <- matrix(runif(B*m),m,B) # generating probabilities
    
    if(prior=="loggamma"){
      param <- c(param,"log_prec_r"=-1)
    }
    param1 <- c(param1,"U"=rep(0,m))
  }else{
    rand_U <- NULL; random <- NULL
  }
  
  if(yre){ ## Year random effect
    Year <- factor(data$Year)
    Year <- plyr::mapvalues(Year, from =levels(Year), to=c(1:length(levels(Year))))
    rand_Y <- matrix(runif(B*(length(levels(Year))-1)),(length(levels(Year))-1),B)
    
    if(prior=="loggamma"){
      param <- c(param, "log_prec_t"=-1)
    }
    param1 <- c(param1,rep(0,"Y"=(length(levels(Year))-1)))
  }else{
    rand_Y <- NULL; Year <- NULL
  }
  
  set.seed(5)
  t1 <- Sys.time()
  fitOptim <- optim(par=param, fn=log_posterior, method = opt_method,
                    control = list(fnscale = -1,maxit=100000,factr=1e8,
                                   trace=1,REPORT=1),
                    hessian = HesP, X=X, Y=Y, 
                    rand_u=rand_u, num.knots=num.knots,z.spline=z.spline,
                    rand_te=rand_te,x1x2=x1x2,num.knotst=num.knotst,
                    rand_U=rand_U, 
                    rand_Y=rand_Y , B=B,random=random, Year=Year,
                    prior=prior
                    )
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
                        x1x2=x1x2,num.knotst=num.knotst,
                        m=m,random=random,mu_prior=mu_prior,Year=Year,prior=prior)
    t2 <- Sys.time()
    time2 <- difftime(t2, t1, units = "secs")
    print(fitOptim_u$convergence)
    print(time2)
    out <- list(fitOptim,fitOptim_u)
  }else{
    out <- list(fitOptim)
  }
  
  return(out)
  
}
