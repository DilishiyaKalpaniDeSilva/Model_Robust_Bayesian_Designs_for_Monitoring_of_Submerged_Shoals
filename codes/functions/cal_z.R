## ---------------------------------------------------------------------------------------------
## Authors: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose: Calculate model matrix for different spline basis
## ---------------------------------------------------------------------------------------------

library(HRW)
library(splines)

# intKnots: interior knots without boundaries
# numBasis: number of knots
# XsplPreds: covariate matrix
# knot_type: equally sapced or quantiles
# boundKnots: boundary knots
# basis: spline basis type


ZOSull1 <- function (x, range.x, intKnots, drv = 0) {
  if (drv > 2) 
    stop("splines not smooth enough for more than 2 derivatives")
  if (missing(range.x)) 
    range.x <- c(1.05 * min(x) - 0.05 * max(x), 1.05 * max(x) - 
                   0.05 * min(x))
  if (missing(intKnots)) {
    numIntKnots <- min(length(unique(x)), 35)
    intKnots <- quantile(unique(x), seq(0, 1, length = (numIntKnots + 
                                                          2))[-c(1, (numIntKnots + 2))])
  }
  numIntKnots <- length(intKnots)
  allKnots <- c(rep(range.x[1], 4), intKnots, rep(range.x[2], 
                                                  4))
  K <- length(intKnots)
  L <- 3 * (K + 8)
  xtilde <- (rep(allKnots, each = 3)[-c(1, (L - 1), L)] + rep(allKnots, 
                                                              each = 3)[-c(1, 2, L)])/2
  wts <- rep(diff(allKnots), each = 3) * rep(c(1, 4, 1)/6, 
                                             K + 7)
  Bdd <- spline.des(allKnots, xtilde, derivs = rep(2, length(xtilde)), 
                    outer.ok = TRUE)$design
  Omega <- t(Bdd * wts) %*% Bdd
  # eigOmega <- eigen(Omega)
  eigOmega <- with(svd(Omega), list(values = d, vectors = u))
  indsZ <- 1:(numIntKnots + 2)
  UZ <- eigOmega$vectors[, indsZ]
  LZ <- t(t(UZ)/sqrt(eigOmega$values[indsZ]))
  indsX <- (numIntKnots + 3):(numIntKnots + 4)
  UX <- eigOmega$vectors[, indsX]
  Lmat <- cbind(UX, LZ)
  stabCheck <- t(crossprod(Lmat, t(crossprod(Lmat, Omega))))
  if (sum(stabCheck^2) > 1.0001 * (numIntKnots + 2)) 
    print("WARNING: NUMERICAL INSTABILITY ARISING\\\n              FROM SPECTRAL DECOMPOSITION")
  B <- spline.des(allKnots, x, derivs = rep(drv, length(x)), 
                  outer.ok = TRUE)$design
  Z <- crossprod(t(B), LZ)
  Z <- Z[, order(eigOmega$values[indsZ])]
  attr(Z, "range.x") <- range.x
  attr(Z, "intKnots") <- intKnots
  return(Z)
}

## O'Sullivan splines 
OSullivan_basis <- function(numBasis ,XsplPreds, knot_type, intKnots, boundKnots){
  
  numSplCompons <- ncol(XsplPreds)
  ncZspl_p <- NULL
  Zspl_p <- NULL
  Knots.list <- list()
  
  if(any(is.na(intKnots))){ ## interior knots not given
    for (jSpl in 1:numSplCompons) {
      
      xCurr <- XsplPreds[, jSpl]
      range.xVal <- c(min(xCurr), max(xCurr))
      
      if(knot_type=="es"){
        equally_spaced_knots <- seq(range.xVal[1], range.xVal[2], length.out = numBasis[jSpl])
        intKnots <- equally_spaced_knots[-c(1, length(equally_spaced_knots))]
      }else{
        numIntKnotsVal <- numBasis[jSpl] - 2
        intKnots <- as.numeric(quantile(unique(xCurr), 
                                        seq(range.xVal[1], range.xVal[2], 
                                            length = numIntKnotsVal + 2)[-c(1, numIntKnotsVal + 2)]))
      }
      
      Zcurr <- ZOSull1(x=xCurr, intKnots = intKnots, range.x = range.xVal)
      Knots.list[[jSpl]] <- c(range.xVal[1],intKnots,range.xVal[2])
      ncZspl_p <- c(ncZspl_p, ncol(Zcurr))
      Zspl_p <- cbind(Zspl_p, Zcurr)
    }
  }else{ ## interior knots given
    
    for (jSpl in 1:numSplCompons) {
      xCurr <- XsplPreds[, jSpl]
      Zcurr <- ZOSull1(xCurr, intKnots = intKnots, range.x = boundKnots)
      # Knots.list[[jSpl]] <- c(range.xVal[1],intKnots,range.xVal[2])
      ncZspl_p <- c(ncZspl_p, ncol(Zcurr))
      Zspl_p <- cbind(Zspl_p, Zcurr)
    }
  }
  
  return(list("Z"=Zspl_p,"k"=ncZspl_p, "knot_loc"=Knots.list))
}

## Generate the B-spline (cubic, degree=3) basis matrix for a polynomial spline.
Bspline_basis <- function(numBasis ,XsplPreds, knot_type, intKnots, boundKnots){
  
  numSplCompons <- ncol(XsplPreds)
  ncZspl_p <- NULL
  Zspl_p <- NULL
  Knots.list <- list()
  
  if(any(is.na(intKnots))){ ## interior knots not given
    for (jSpl in 1:numSplCompons) {
      
      xCurr <- XsplPreds[, jSpl]
      range.xVal <- c(min(xCurr), max(xCurr))
      
      if(knot_type=="es"){
        equally_spaced_knots <- seq(range.xVal[1], range.xVal[2], length.out = numBasis[jSpl])
        intKnots <- equally_spaced_knots[-c(1, length(equally_spaced_knots))]
      }else{
        numIntKnotsVal <- numBasis[jSpl] - 2
        intKnots <- as.numeric(quantile(unique(xCurr), 
                                        seq(range.xVal[1], range.xVal[2], 
                                            length = numIntKnotsVal + 2)[-c(1, numIntKnotsVal + 2)]))
      }
      
      Zcurr <- splines::bs(xCurr, knots = intKnots, degree = 3,intercept = FALSE)
      Knots.list[[jSpl]] <- c(range.xVal[1],intKnots,range.xVal[2])
      ncZspl_p <- c(ncZspl_p, ncol(Zcurr))
      Zspl_p <- cbind(Zspl_p, Zcurr)
    }
  }else{ ## interior knots given
    
    for (jSpl in 1:numSplCompons) {
      xCurr <- XsplPreds[, jSpl]
      Zcurr <- splines::bs(xCurr, knots = intKnots, degree = 3,intercept = FALSE)
      # Knots.list[[jSpl]] <- c(range.xVal[1],intKnots,range.xVal[2])
      ncZspl_p <- c(ncZspl_p, ncol(Zcurr))
      Zspl_p <- cbind(Zspl_p, Zcurr)
    }
  }
  return(list("Z"=Zspl_p,"k"=ncZspl_p, "knot_loc"=Knots.list))
}

## Generate the B-spline basis matrix for a natural cubic spline.
natural_cubic_Bspline_basis <- function(numBasis ,XsplPreds, knot_type, intKnots, boundKnots){
  
  numSplCompons <- ncol(XsplPreds)
  ncZspl_p <- NULL
  Zspl_p <- NULL
  Knots.list <- list()
  
  if(any(is.na(intKnots))){ ## interior knots not given
    for (jSpl in 1:numSplCompons) {
      
      xCurr <- XsplPreds[, jSpl]
      range.xVal <- c(min(xCurr), max(xCurr))
      
      if(knot_type=="es"){
        equally_spaced_knots <- seq(range.xVal[1], range.xVal[2], length.out = numBasis[jSpl])
        intKnots <- equally_spaced_knots[-c(1, length(equally_spaced_knots))]
      }else{
        numIntKnotsVal <- numBasis[jSpl] - 2
        intKnots <- as.numeric(quantile(unique(xCurr), 
                                        seq(range.xVal[1], range.xVal[2], 
                                            length = numIntKnotsVal + 2)[-c(1, numIntKnotsVal + 2)]))
      }
      
      Zcurr <- splines::ns(xCurr, knots = intKnots, Boundary.knots = range.xVal)
      Knots.list[[jSpl]] <- c(range.xVal[1],intKnots,range.xVal[2])
      ncZspl_p <- c(ncZspl_p, ncol(Zcurr))
      Zspl_p <- cbind(Zspl_p, Zcurr)
    }
  }else{ ## interior knots given
    
    for (jSpl in 1:numSplCompons) {
      xCurr <- XsplPreds[, jSpl]
      Zcurr <- splines::ns(xCurr, knots = intKnots, Boundary.knots = boundKnots)
      # Knots.list[[jSpl]] <- c(range.xVal[1],intKnots,range.xVal[2])
      ncZspl_p <- c(ncZspl_p, ncol(Zcurr))
      Zspl_p <- cbind(Zspl_p, Zcurr)
    }
  }
  
  return(list("Z"=Zspl_p,"k"=ncZspl_p, "knot_loc"=Knots.list))
  
}

cal_z <- function(numBasis ,XsplPreds, knot_type, intKnots, boundKnots, basis){
  
  if(basis=="OS"){
    OSullivan_basis(numBasis ,XsplPreds, knot_type, intKnots, boundKnots)
  }else if(basis=="BS"){
    Bspline_basis(numBasis ,XsplPreds, knot_type, intKnots, boundKnots)
  }else if(basis=="NCBS"){
    natural_cubic_Bspline_basis(numBasis ,XsplPreds, knot_type, intKnots, boundKnots) 
  }else{
    print("Not Specified")
  }
  
}