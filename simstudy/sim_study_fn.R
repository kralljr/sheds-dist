getx <- function(ny, na, argvals1, typex = "shift", mean1 = 15, sd1 = 1.5, rate1 = 1/2, ns1 = 10) {
  # Specify number of arguments/quantiles
  # Specify total number of x (quantiles X days)
  N <- na * ny
  #x1 <- matrix(nrow = na, ncol = ny)

  # For shift in distribution
  if(typex == "shift") {
    # Define shift, varies for each day
    #rates1 <- rep(rexp(ny, rate = rate1), each = na)
    rates1 <- rep(runif(ny, min = -5, max = 20), each = na)
    # Base dist truncnorm, plus rates vary for each day
    x0 <- rtruncnorm(N, a = 0, mean = mean1, sd = sd1)
    xall <- matrix(x0 + rates1, nrow = na, byrow = F)
  
  # For long right tail
  } else if (typex == "longr") {
    # Find base distribution
    x1 <- rtruncnorm(N, a = 0, mean = mean1, sd = sd1)
    med <- median(x1)
    med <- 15
    # Add in right tail
    # q3 <- quantile(x1, probs = 0.75)
    # What to scale right tail by
    #adds <- rexp(ny, rate = rate1)
    adds <- runif(ny, min = 2, max = 6)
    # Same for each day 
    x0 <- rep(adds, each = na)
    # Do not scale lower values
    x0 <- (x1 > med) * ((x1 - med) * x0)
    # Find xall
    xall <- matrix(x0 + x1, nrow = na, byrow = F)



  # For shifted left tail
  } else if (typex == "longl") {


    # Find base distribution
    x1 <- rtruncnorm(N, a = 0, mean = mean1, sd = sd1)
    med <- median(x1)
    med <- 15
    # Add in left tail
    # q3 <- quantile(x1, probs = 0.25)
    # What to scale left tail by
    #adds <- rexp(ny, rate = 1 / rate1)
    adds <- runif(ny, min = 0, max = 2)
    # Same for each day 
    x0 <- rep(adds, each = na)
    # Do not scale lower values
    x0 <- (x1 < med) * ((-x1 + med) * -x0)
    # Find xall
    xall <- matrix(x0 + x1, nrow = na, byrow = F)

  # For increased variance
  } else if (typex == "wide") {

    # Find standard deviations
    sd2 <- rep(runif(ny, min = 0.3, max = 6), each = na) 
    # Get x as truncnorm
    x0 <- rtruncnorm(N, a = 0, mean = 15, sd = sd2)
    # Make matrix
    xall <- matrix(x0, nrow = na, byrow = F)

  } else {

    stop("typex not recognized")
  
  }
  
  
  # Find mean
 xM <- apply(xall, 2, mean)
 # xall <- sweep(xall, 2, xM)
  
  # Find quantiles
  x1 <- apply(xall, 2, quantile, probs = argvals1)
  x1 <- abs(x1)

  # Get functional x
  xfn1 <- getxfn(x1, argvals1, ns1)
  xfn <- xfn1$xfn
  basis1 <- xfn1$basis1

  # Get outcome 
  return(list(x1 = x1, xall = xall, xfn = xfn, basis1 = basis1, xM = xM))

}



# Function to get functional x
getxfn <- function(xvar1, argvals1, ns1 = 15) {
  # Get bspline basis over 0,1 
  basis1 <- create.bspline.basis(c(0, 1), norder = ns1)
  # Get smooth of x
  xfn <- smooth.basis(argvals1, xvar1, basis1)$fd
  list(xfn = xfn, basis1 = basis1)
}




# beta function of x
getbeta <- function(type, val = 0) { function(x) {
  # Beta constant over quantiles
  if(type == "constant") {
    b1 <- rep(val, length = length(x))
  
  # Beta increases for lower & higher quantiles
  } else if (type == "x2") {
    b1 <- val + 1/4 * (x - 0.5)^2 

    #b1 <- val + 1 / 3 *  (x - 0.5)^2 
  # Beta larger for low quantiles
  } else if (type == "low") {
    b1 <- val + 1 / 10 * exp(x * -7)

    #b1 <- val + 1 / 5 * exp(x * -7)
  # Beta larger for high quantiles
  } else if (type == "high") {
    b1 <- val + 1 / 10000 * exp(x * 7) 	

    #b1 <- val + 1 / 5000 * exp(x * 7) 	
  # Beta not specified  
  } else {
    stop("Beta type not recognized")
  }

  #rescale for appropriately sized beta
  b1 <- b1 / 10
}}


gety <- function(argvals1, betaM, betaf, x1, disttype, sd1 = 0.01) {
  xvar1 <- x1$x1
  xM <- x1$xM
  # Get values of beta at x

  # If FDA is truth
  if(class(betaf) == "function") {
    beta1 <- betaf(argvals1)


    # find linear function of x and beta	
    linf <- rowSums(sweep(t(xvar1), 2, beta1, "*"))
    linf <- linf * 1 / length(beta1)
    #linf <- apply(linf, 2, function(x) auc(argvals1, x))

  # If other is truth
  } else{
    nums <- as.numeric(names(betaf))/100
    
    xhold <- apply(x1$xall, 2, quantile, probs = nums)
    if(length(nums) > 1) {
      xhold <- t(xhold)
      
      linf <- rowSums(sweep(xhold, 2, betaf, "*"))
      linf <- linf * 1 / length(beta1)
    }else{
      linf <- betaf * xhold
    } 

  }


  # Add in median
  linf <- linf + xM * betaM


  # For normally dist outcome
  if(disttype == "norm") {
    # get additive error
    eps <- rnorm(ncol(xvar1), sd = sd1)
    # compute y
    y1 <- linf + eps

  # For count outcome
  } else if(disttype == "pois") {
    # Find mean
    mu <- exp(linf)

    # Get poisson
    y1 <- rpois(length(mu), mu)
  }  
  y1
}



simout <- function(x1, argvals1, betaM, typeb, disttype = "norm", sd1 = 0.01, argvalslr = argvals1, val1 = 1, std = F, quants = F,...) {
  # Get function of beta
  if(class(typeb) != "numeric") {
    betaf <- getbeta(typeb, val = val1)
  } else {
    betaf <- typeb
  }
  # Generate y
  y1 <- gety(argvals1, betaM, betaf, x1, disttype, sd1)
  # Get functional x
  xfn <- x1$xfn
  ns1 <- x1$basis1$nbasis

  # Traditional univariate and multivariate regression
  # Set up data frame
  if(!quants) {
	  xmat <- apply(x1$xall, 2, quantile, probs = argvalslr )
	  if(length(argvalslr) > 1) {
		  xmat <- t(xmat)
	  }
  }else{
	  xmat <- t(x1$x1)
  }

  # Standardize?
  if(std) {
   mn1 <- apply(xmat, 2, mean)
   sd1 <- apply(xmat, 2, sd)
   xmat <- sweep(xmat, 2, mn1, "-")
   xmat <- sweep(xmat, 2, sd1, "/")  
  }


  dat1 <- data.frame(y1, xmat)

  # do multivariate regression
  colnames(dat1) <- c("y", paste0("x", seq(1, ncol(xmat))))
  eqn1 <- paste0("y ~", paste(colnames(dat1)[-1], collapse = "+"))

  beta2 <- matrix(nrow = (ncol(dat1) - 1), ncol = 4)
    
  # Depending on type of regression
  if(disttype == "norm") {  
    fmod1 <- flm(x1, y1)

    beta3 <- summary(lm(eval(eqn1), data = dat1))$coef[-1, ] 


    # Do univariate regression
    for(i in 2 : ncol(dat1)) {
      eqn1 <- paste("y ~", colnames(dat1)[i])
      beta2[i- 1, ] <- summary(lm(eval(eqn1), data = dat1))$coef[-1, ]
    }


  } else if (disttype == "pois") {
    #fmod1 <- fglm1(x1, y1, argvals1, ns1)
    fmod1 <- NULL
    beta3 <- summary(glm(eval(eqn1), data = dat1, family = "poisson"))$coef[-1, ] 
  
  

    # Do univariate regression
    for(i in 2 : ncol(dat1)) {
      eqn1 <- paste("y ~", colnames(dat1)[i])
      beta2[i- 1, ] <- summary(glm(eval(eqn1), data = dat1, family = "poisson"))$coef[-1, ]
    }
  
  }


  rownames(beta2) <- argvalslr
  rownames(beta3) <- argvalslr 



  #freg1 <- fRegress(y1 ~ x1$xfn)

  # Save output
  list(y1 = y1, betaf = betaf, fmod1 = fmod1, beta2 = beta2, beta3 = beta3, basis1 = x1$basis1)
}



flm <- function(x1, y1) {
  xfn <- x1$xfn
  xM <- x1$xM
  ny <- length(y1)

  # Get beta
  # The intercept must be constant for a scalar response
  betabasis1 <- create.constant.basis(c(0, 1))
  betafd1    <- fd(0, betabasis1)
  betafdPar1 <- fdPar(betafd1)

  betafd2     <- with(xfn, fd(basisobj=basis, fdnames=fdnames))
  # convert to an fdPar object
  betafdPar2  <- fdPar(betafd2)
  betalist <- list(const=betafdPar1, xM = betafdPar1, xfn=betafdPar2)

  # Get x
  xfdlist <- list(const=rep(1, ny), xM = xM, xfn=xfn)
  
  # Do functional regression
  fd1 <- fRegress(y1, xfdlist, betalist)

  # Find corresponding CIs
  yhatfdobj   <- fd1$yhatfdobj
  errmat <- y1 - yhatfdobj
  sigmae <- as.numeric(var(errmat))
  diag1 <- diag(1, length(y1))
  std1 <- fRegress.stderr(fd1, diag1, diag1 * sigmae)

  list(freg = fd1, betafstd = std1)
}



fglm1 <- function(x1, y1, argvals1, ns1) {
  xM <- x1$xM
  form1 <- formula(y1 ~ x + xM)
 

  basisx <- create.bspline.basis(c(0, 1), norder = ns1)
  basisb <- create.bspline.basis(c(0, 1), norder = ns1)

  #basis1 <- x1$basis1
  basx <- list(x = basisx)
  basb <- list(x = basisb)

  xfn <- fdata(t(x1$x1), argvals = argvals1)

  y1 <- data.frame(y1, xM)
  dat1 <- list("x" = xfn, "df" = y1) 
  
  
  fre1 <- fregre.glm(form1, family = "poisson", data = dat1,
    basis.x = basx, basis.b = basb, CV = F)
  fre1
}




fglm <- function(x1, y1, argvals1, ns1) {
  pfr1 <- pfr(y1, funcs = t(x1$x1), kz = ns1, nbasis = ns1, kb = ns1, family = "quasipoisson" )
  pfr1
}
