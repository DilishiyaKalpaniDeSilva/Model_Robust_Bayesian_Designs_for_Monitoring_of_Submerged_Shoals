library(matrixcalc)
library(VGAM)

## log posterior for fixed effect and the variance parameters
log_posterior <- function(theta,X,Y,year,rand_u,rand_U,rand_t,z.spline,fishnet){
  
  pb <- sum(dnorm(x=theta[1:ncol(X)],mean=rep(0,ncol(X)),sd=rep(5,ncol(X)),log=TRUE))
  ps <- dlgamma(theta[3],location=-33,scale = 5, shape = 10,log = TRUE) 
  pr <- dlgamma(theta[4],location=-21,scale = 5, shape = 10,log = TRUE)  
  pt <- dlgamma(theta[5],location=-21,scale = 5, shape = 10,log = TRUE)  
  log_prior <- pb + ps + pr + pt
  
  u <- qnorm(rand_u,0, 1/(sqrt(exp(theta[3]))) )
  U <- qnorm(rand_U,0, 1/(sqrt(exp(theta[4]))) )
  tre <- qnorm(rand_t,0, 1/(sqrt(exp(theta[5]))) )
  tre1 <- rbind(0,tre)
  
  rand <- z.spline %*% u + U[fishnet,] + tre1[year,]
  x_beta <- as.vector(X %*% theta[1:ncol(X)] ) #### Historical plus simulated data
  prob <-  expit(x_beta + rand ) 
  like <- colSums( (dbinom( x=Y, size=20,  prob = prob, log= TRUE) ) )
  log_like <- mean(like)
  
  log_post<- log_like+log_prior
  if(is.infinite(log_post) || is.na(log_like)){
    return(10000000)
  }else{
    return(-log_post)
  }
  
}

kld_util <- function(mu_prior_all,sigma_prior_all,mu_post_all,sigma_post_all){
  
  num_par <- length(mu_prior_all)
  kld <- 0.5 * (  matrix.trace( (matrix.inverse(sigma_prior_all)) %*% sigma_post_all ) + 
                    t(mu_prior_all - mu_post_all) %*% (matrix.inverse(sigma_prior_all)) %*% (mu_prior_all-mu_post_all) -
                    num_par + 
                    (as.vector(determinant(sigma_prior_all)$modulus) - as.vector(determinant(sigma_post_all)$modulus) ) ) 
  return(kld)
  
}

## function to approximate KLD utility 
eval_utility <- function(theta,X,Y,year,rand_u,rand_U,rand_t,z.spline,fishnet){
  
  # minimise negative log posterior
  fitOptim_post <- optim(par=theta, fn =log_posterior, method = "L-BFGS-B",
                         control = list(maxit=10000000), hessian = TRUE,
                         X=X,rand_u=rand_u,rand_U=rand_U,z.spline=z.spline,
                         Y=Y,fishnet=fishnet,year=year, rand_t=rand_t
  )
  
  mu_post <- c(fitOptim_post$par)
  sigma_post <- solve(fitOptim_post$hessian)
  
  kld <- kld_util(mu_prior,sigma_prior,mu_post,sigma_post)
  
  return(kld)
  
}

exp_crit <- function(d,B){
  
  Xd <- na.omit(unlist(X_all[d]))  
  Xd <- cbind(1,Xd)
  fishnetd <- na.omit(unlist(fishnetp[d]))
  zz <- cal_z(numBasis = num.knots,XsplPreds=as.matrix((Xd[,2])),knot_type="es" ,intKnots,
              boundKnots=range(data$depth), basis ="OS")
  z.spline <- zz$Z
  rand1 <- (z.spline %*% t(u_all)) + t(U_all[,fishnetd]) 
  rand <- sweep(rand1, 2, as.vector(t_all), "+") 
    
  x_beta <- (Xd %*% t(theta_prior[,1:ncol(Xd)]) )
  pred_p <- expit(x_beta + rand)
  
  Xd_all <- rbind(Xd,Xd_hist)
  zz_all <- cal_z(numBasis = num.knots,XsplPreds=as.matrix(Xd_all[,2]),knot_type="es" ,intKnots=NA,
                  boundKnots=range(Xd_all[,2]), basis ="OS")
  z.spline_all <- zz_all$Z
  fishnet_all <- c(fishnetd,data$FishNet)
  y_all <- as.numeric(as.factor(c(data$Year,rep(2024,nrow(Xd)))))
  
  # apply parallel processing 
  crit  <- foreach(j = 1:B,
                   .packages=c("plyr","matrixcalc","Matrix","mvtnorm","doParallel","LaplacesDemon",
                               "HRW","stats","VGAM"),
                   .combine = c,
                   .noexport = c("expitCpp","dbinomCpp"),
                   #.errorhandling = 'remove',
                   .export = ls(globalenv())) %dopar%  {
                     
                     HardCoral_all <- c(qbinom(p=rand_y[j,1:nrow(pred_p)],size=20, prob=pred_p[,j]), data$HardCoral )
                     
                     #t1 <- Sys.time()
                     util <- eval_utility(X=Xd_all,Y=HardCoral_all,theta=c(5,5,-1,-1,-1),
                                          rand_u=rand_u,rand_U=rand_U,rand_t=rand_t,z.spline=z.spline_all, 
                                          year=y_all,
                                          fishnet=fishnet_all)
                     
                     util 
                   }
  crit
}
