expit<-function(a){
  exp(a)/(1+exp(a))
}

logit <- function(a){
  log(a/(1-a))
}

norm_01 <- function(a){
  min <- min(a)
  max <- max(a)
  ans <-lapply(seq_len(length(a)), function(x)
    (a[x] - min ) / ( max - min)
  )
  return(unlist(ans))  
}

stand <- function(a){
  return((a-mean(a))/sd(a))
}

stand_m_sd <- function(a,mean,sd){
  return((a-mean)/sd)
}

norm_01_mm <- function(a,min,max){
  ans <-lapply(seq_len(length(a)), function(x)
    (a[x] - min ) / ( max - min)
  )
  return(unlist(ans))  
}

scaleP <- function(x, oldmin, oldmax, newmin, newmax){
  old_range <- oldmax - oldmin
  new_range <- newmax - newmin
  newval <- (((x -oldmin) * new_range) / old_range) + newmin
  return(newval)
}


trace <- function(M){
  tr <- sum(diag(M))
  tr
}

## Calculate euclidean distance of vector
euc_dist <- function(x) 
  sqrt(sum(x^2))

vector.is.empty <- function(x) 
  return(length(x) ==0 )

# convert Binomial data into binary data
as_binary <- function(variable) {
  vec1 <- NULL
  for(i in 1:length(variable)){
    vec <- c(rep(1,variable[i]),rep(0,(20-variable[i])))
    vec <- matrix(vec,nrow = 20,ncol = 1,byrow = TRUE)
    vec1 <- rbind(vec1,vec)
  }
  return(vec1)
}

# convert the full dataset into binary
convert_data_to_binary <- function(data,response){
  variables <- names(data)
  assign(response, as_binary(data[,response]))
  HardCoral <- get(response)
  other_variables <- variables[-which(variables==response)] 
  out <- data.frame(data[,other_variables]) 
  out1 <- out[rep(1:nrow(data), each=20),]
  data_bin <- cbind(out1,HardCoral)
  return(data_bin)
}

#plotting x vs logit p for original data by dividing data into categories
plot_x_vs_logit_p <- function(Mdata_cov,variable){
  
  t <- with(Mdata_cov, table(variable,Mdata_cov$HardCoral))
  rs<-rowSums(t, na.rm=TRUE)
  cs<-colSums(t,na.rm = TRUE)
  
  t <- cbind(t,rs)
  t <- rbind(t,"cs1"=c(cs,sum(rs)))
  
  #print(t) # contigency table for HardCoral and the independet variable
  
  p <- matrix(NA,nrow = length(rs),ncol =1)
  logit <- matrix(NA,nrow = length(rs),ncol = length(cs))
  mean_p <- vector()
  
  for (i in 1:length(rs)) {
    p[i] <- (rs[i])/(sum(rs)*20)
    logit[i] <- log(p[i]/(1-p[i]))
    mean_p[i] <- logit[i]
  }
  
  logit_mean_p <- mean_p
  data <- data.frame(depth=names(rs), logit_p = logit_mean_p, p = p)
  
  return(data)
  
}

f.pred<-function(x,fit)  #### Change how you predict from GP
{
  nd<-matrix(x,nrow=1)  
  #Predicted values and (marginal of joint) conditional variances based on a km model. 95 % confidence intervals are given, based on strong assumptions: Gaussian process assumption, specific prior distribution on the trend parameters, known covariance parameters. This might be abusive in particular in the case where estimated covariance parameters are plugged in
  ans<-predict(fit,newdata=data.frame(z=nd),type="UK",cov.compute=TRUE,se.compute=TRUE)
  out <- ans$mean
  
  return(out)
}

## function for Bayesian model probability using Laplace approximation
Bayes_model_selec_lap <- function(fitOptim){
  N <- length(fitOptim$par)
  H <- -fitOptim$hessian             ## using positive log posterior in optim, therefore use negative
  log_pDm <- fitOptim$value + (N/2)*log(2*pi) + (1/2)*log(det(solve(H)))
  return(log_pDm)
}