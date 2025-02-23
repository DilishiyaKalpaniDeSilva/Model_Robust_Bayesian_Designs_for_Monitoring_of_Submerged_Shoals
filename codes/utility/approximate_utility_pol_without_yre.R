library(matrixcalc)
library(VGAM)

## log posterior for fixed effect and the variance parameters
log_posterior <- function(theta,X,Y,rand_U,fishnet){
  
  pb <- sum(dnorm(x=theta[1:ncol(X)],mean=rep(0,ncol(X)),sd=rep(5,ncol(X)),log=TRUE))
  pr <- dlgamma(theta[indexO+2],location=-21,scale = 5, shape = 10,log = TRUE)  
  log_prior <- pb + pr
  
  U <- qnorm(rand_U,0, 1/(sqrt(exp(theta[indexO+2]))) )
  
  rand <- U[fishnet,] 
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
eval_utility <- function(X,Y,theta,rand_U,B1,fishnet){
  
  # minimise negative log posterior
  fitOptim_post <- optim(par=theta, fn =log_posterior, method = "L-BFGS-B",
                         control = list(maxit=10000000), hessian = TRUE,
                         X=X,rand_U=rand_U,Y=Y,fishnet=fishnet
  )
  
  mu_post <- c(fitOptim_post$par)
  sigma_post <- solve(fitOptim_post$hessian)
  
  kld <- kld_util(mu_prior,sigma_prior,mu_post,sigma_post)
  return(kld)
  
}

exp_crit <- function(d,B){
  
  Xd <- na.omit(unlist(X_all[d]))  
  Xd <- cbind( 1,outer(X = Xd,Y = 1:indexO,FUN = "^"))
  fishnetd <- na.omit(unlist(fishnetp[d]))
  rand <- t(U_all[,fishnetd])
  x_beta <- (Xd %*% t(theta_prior[,1:ncol(Xd)]) )
  pred_p <- expit(x_beta + rand)
  
  Xd_all <- rbind(Xd,Xd_hist)
  fishnet_all <- c(fishnetd,data$FishNet)
  
  # apply parallel processing 
  crit  <- foreach(j = 1:B,
                   .packages=c("plyr","matrixcalc","Matrix","mvtnorm","doParallel","LaplacesDemon",
                               "HRW","stats","VGAM"),
                   .combine = c,
                   .noexport = c("expitCpp","dbinomCpp"),
                   #.errorhandling = 'remove',
                   .export = ls(globalenv())) %dopar%  {
                     
                     HardCoral_all <- c(qbinom(p=rand_y[j,1:nrow(pred_p)],size=20, prob=pred_p[,j]), data$HardCoral )
                     
                     util <- eval_utility(X=Xd_all,Y=HardCoral_all,theta=c(rep(5, indexO+1),-1),
                                          rand_U=rand_U,fishnet=fishnet_all)
                     
                     util 
                   }
  crit
}
