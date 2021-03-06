# File to fit Howard's biostatistics model
# Adapted from AutoReg.R in WindowsRcode

#' Function to perform MCMC
#'  
#' @param y vector of outcomes (counts) length n
#' @param x matrix of quantiles (n X p)
#' @param guessvec list of items in mcmc
#' @param tunes list of tuning parameters
#' @param hyperp list of hyperparameters
#' @param niter number of MCMC iterations
#' @param burnin length of burnin period
#' @param thin Amount to thin out
mcmcout <- function(y, x, quants, guessvec = NULL, tunes = NULL, hyperp = NULL,
  niter = 100, burnin = 0, thin = 1, type = "pois", theta = T, xmean = NULL) {

  # Get guessvec
  if(is.null(guessvec)) {
    guessvec$beta0 <- 1
    guessvec$betam <- 0
    b1 <- t(matrix(rep(1, ncol(x))))
    guessvec$beta1 <- b1

    guessvec$phi <- 1
    # works!
    guessvec$theta <- 0.00035698
    guessvec$sigma2 <- 0.000001
    if(!theta) { 
      guessvec$theta <- 0
    }
  }
  
  # Add in data
  guessvec$y <- y
  guessvec$x <- x
  guessvec$xmean <- xmean

  # Add in dist
  guessvec$type <- type


  # Get distance matrix
  np <- ncol(x)
  seq1 <- rep(quants, each = np)
  rows <- matrix(seq1, nrow = np, byrow = T)
  cols <- matrix(seq1, nrow = np, byrow = F)
  guessvec$Dists <- abs(rows - cols)

  # Get tuning parameters
  if(!is.null(tunes)) {
    beta0.tune <- tunes$beta0.tune
    betam.tune <- tunes$betam.tune
    beta1.tune <- tunes$beta1.tune
    phi.tune <- tunes$phi.tune
  } else{
    # for simulated data
    beta0.tune <- 0.1
    beta1.tune <- 0.8
    phi.tune <- 1.5
    
    beta0.tune <- 0.001
    betam.tune <- 0.001
    beta1.tune <- .0001
    phi.tune <- 0.5
      
    # From howard 
    #phi.tune <- 0.5


  }

  # Keep track of acceptance
  guessvec$accept <- c(0, 0, 0, rep(0, ncol(x)))

  # Get hyperparameters
  if(!is.null(hyperp)) {
    sd.beta0 <- hyperp$sd.beta0
    sd.betam <- hyperp$sd.betam
    a.sig <- hyperp$a.sig
    b.sig <- hyperp$b.sig
    a.phi <- hyperp$a.phi
    b.phi <- hyperp$b.phi
    mu.theta <- hyperp$mu.theta
    sig2.theta <- hyperp$sig2.theta
  } else {
    sd.beta0 <- 100 
    sd.betam <- 100
    #sd.beta0 <- 0.001
    
    # Howard: too diffuse!
    a.sig <- 0.001
    b.sig <- 0.001
    # Worked for me
    #a.sig <- 100
    #b.sig <- 10

    # Howard: distance units 1/10??
    #a.phi <- 0.03
    #b.phi <- 0.005
    # Worked for me
    a.phi <- 9
    b.phi <- 10
  
    # Try lower cor
    #a.phi <- 3
    #b.phi <- 10

    mu.theta <- 0
    sig2.theta <- 100^2
  }
  


  # Setup output
  # Find number of output samples
  nsamp <- ceiling((niter - burnin) / thin)
  beta0.out <- vector(, length = nsamp)
  betam.out <- vector(, length = nsamp)
  beta1.out <- matrix(nrow = nsamp, ncol = length(guessvec$beta1)) 
  sigma2.out <- vector(, length = nsamp)
  phi.out <- vector(, length = nsamp)
  theta.out <- vector(, length = nsamp)

  # For niter iterations, do MCMC
  # where to save
  k <- 1
  l <- 1
  for (i in 1 : niter) {
    # Update beta0
    guessvec <- beta0f(guessvec, sd.beta0, beta0.tune, quants)
    # Update betam
    if(!is.null(xmean)) {
      guessvec <- betamf(guessvec, sd.betam, betam.tune, quants)
    }
    # Update beta1
    guessvec <- beta1f(guessvec, beta1.tune, quants)
    # Update sigma2
    guessvec <- sigma2f(guessvec, a.sig, b.sig)
    # Update phi	  
    guessvec <- phif(guessvec, a.phi, b.phi, phi.tune)
    # Update theta
    if(theta) {
      guessvec <- thetaf(guessvec, mu.theta, sig2.theta) 
    }
    #save lth iteration of guesses
    if(i > burnin) {
      
      # If multiple of thin, save output
      if(((k - 1) %% thin) == 0) { 
        beta0.out[l] <- guessvec$beta0
        betam.out[l] <- guessvec$betam
        beta1.out[l, ] <- guessvec$beta1
        sigma2.out[l] <- guessvec$sigma2   
        phi.out[l] <- guessvec$phi
        theta.out[l] <- guessvec$theta
	l <- l + 1	
      }

      # Print count every 100
      if(i %% 1000 == 0) {
        print(i)
      }

      k <- k + 1
    }
  }
  # Get acceptance proportion
#  accept <- guessvec$accept / niter
  accept <- guessvec$accept
  names(accept) <- c("beta0", "phi", paste0("beta1_",1 : ncol(x) ))


  # Get history of samples
  out <- list(beta0 = beta0.out, betam = betam.out, beta1 = beta1.out, sigma2 = sigma2.out, 
    phi = phi.out, theta = theta.out, accept = accept)
  return(out)

}




#' Function to update beta0
#'
#' @param guessvec list of items in mcmc
#' @param sd.beta0 standard deviation for normal prior
#' @param beta0.tune tuning parameter for random walk 
beta0f <- function(guessvec, sd.beta0, beta0.tune, quants) {
  # Get guesses
  beta0 <- guessvec$beta0

  # Propose new beta
  beta0.prop <- rnorm(1, beta0, beta0.tune)

  # Get new guess
  guessvec.new <- guessvec
  guessvec.new$beta0 <- beta0.prop

  # Get likelihoods
  llhood.old <- llhood.beta0(guessvec, sd.beta0, quants)
  llhood.new <- llhood.beta0(guessvec.new, sd.beta0, quants)

  # Select old vs. new guess
  guessvec <- mhstep(guessvec, guessvec.new, llhood.old, llhood.new, j = 1)
  
  return(guessvec)
}




#' Function to update betam
#'
#' @param guessvec list of items in mcmc
#' @param sd.betam standard deviation for normal prior
#' @param betam.tune tuning parameter for random walk 
betamf <- function(guessvec, sd.betam, betam.tune, quants) {
  # Get guesses
  betam <- guessvec$betam

  # Propose new beta
  betam.prop <- rnorm(1, betam, betam.tune)

  # Get new guess
  guessvec.new <- guessvec
  guessvec.new$betam <- betam.prop

  # Get likelihoods: same lhood as beta0
  llhood.old <- llhood.beta0(guessvec, sd.betam, quants)
  llhood.new <- llhood.beta0(guessvec.new, sd.betam, quants)

  # Select old vs. new guess
  guessvec <- mhstep(guessvec, guessvec.new, llhood.old, llhood.new, j = 2 + length(guessvec$beta1) + 1)
  
  return(guessvec)
}


#' Function to update beta1
#'
#' @param guessvec list of items in mcmc
#' @param beta1.tune tuning covariance for random walk 
beta1f <- function(guessvec, beta1.tune, quants) {
  # Get guesses
  beta1 <- guessvec$beta1
  phi <- guessvec$phi
  sigma2 <- guessvec$sigma2
  Dists <- guessvec$Dists  

  # Set tune based on correlation
#  Sigma <- beta1.tune * exp(-1 / phi * Dists)
  #beta1.tune <- Sigma

  # Propose new beta
  # beta1.prop <- rmvnorm(1, beta1, Sigma)


  # Get new guess
  guessvec.new <- guessvec 
  # guessvec.new$beta1 <- beta1.prop
  for(i in 1 : length(beta1)) {
    # update beta1
    beta1 <- guessvec$beta1
	  
    # find conditional mean and variance
#    iSig <- chol2inv(chol(Sigma[-i, -i]))
#    s1 <- Sigma[i, -i] %*% iSig
#    mean1 <- s1 %*% matrix(beta1[-i])
    #sd1 <- sqrt(Sigma[i, i] - s1 %*% Sigma[-i, i])
    beta1.prop <- rnorm(1, mean = (beta1[i]), sd = beta1.tune)
   
    guessvec.new$beta1[i] <- (beta1.prop)
   
    # Get likelihoods
    llhood.old <- llhood.beta1.out(guessvec, quants)
    llhood.new <- llhood.beta1.out(guessvec.new, quants)

    # Select old vs. new guess
    guessvec <- mhstep(guessvec, guessvec.new, llhood.old, llhood.new, j = 2 + i)
   
    guessvec.new <- guessvec 

  
  }




 return(guessvec)
}





# Function to sample theta
#' @param guessvec list if items in mcmc
#' @param mu.theta prior mean for theta
#' @param sig2.theta prior variance for theta
thetaf <- function(guessvec, mu.theta, sig2.theta) {
   
  # First get guessvec of importance
  Dists <- guessvec$Dists
  phi <- guessvec$phi
  beta1 <- guessvec$beta1
  sigma2 <- guessvec$sigma2

  # Get covariance matrix and inverse
  Sigma1 <- sigma2 * exp(-1/phi * Dists)
  sSigma1 <- chol2inv(chol(Sigma1))


  # Find posterior variance
  postvar <-  1/ (1 / sig2.theta + sum( sSigma1))
  
  # Find posterior mean
  meanp <- (mu.theta / sig2.theta) + sum(beta1 %*% sSigma1)
  mean1 <- postvar * meanp

  # Sample new theta
  ntheta <- rnorm(1, mean1, sqrt(postvar))

  guessvec$theta <- ntheta

  return(guessvec)
}


# Function to sample sigma2
#' @param guessvec list of items in mcmc
#' @param a.sig prior shape for gamma for (sigma^2)^{-1}
#' @param b.sig prior rate for gamma for (sigma^2)^{-1}
sigma2f <- function(guessvec, a.sig, b.sig) {
  
  # First get guessvec of importance
  Dists <- guessvec$Dists
  phi <- guessvec$phi
  beta1 <- guessvec$beta1
  sigma2 <- guessvec$sigma2
  theta <- guessvec$theta

  # Get correlation,covariance matrix
  C1 <- exp(-1/phi * Dists)


  # Find scaling factor
  C2 <- chol(chol2inv(chol(C1)))

  # Find scaled beta1
  gamma1 <-  t(C2) %*% t((beta1))
  mu1 <- rowSums(t(C2)) * theta

  n <- length(gamma1)

  # Sample posterior (normal lhood, gamma prior)
  gam1 <- rgamma(1, a.sig + n / 2, b.sig + sum((gamma1 - mu1)^2)/2)
  # Get inverse gamma
  sig2 <- 1 / gam1
  # Update guessvec
  guessvec$sigma2 <- sig2
  
  # Return guessvec
  return(guessvec)
}






#' Function to update phi
#'
#' @param guessvec list of items in MCMC
#' @param a.phi prior shape for gamma for phi
#' @param b.phi prior rate for gamma for phi
#' @param phi.tune tuning parameter for random walk 
phif <- function(guessvec, a.phi, b.phi, phi.tune) {

  # Get guessvec
  phi <- guessvec$phi

  # Propose lognormal random walk
  phi.prop <- rlnorm(1, log(phi), phi.tune)

  # Get new guess
  guessvec.new <- guessvec
  guessvec.new$phi <- phi.prop

  # Get likelihoods
  llhood.old <- llhood.phi(guessvec, a.phi, b.phi)
  llhood.new <- llhood.phi(guessvec.new, a.phi, b.phi)

  # Add in proposal (bc not symmetric)
  llhood.old <- llhood.old + log(phi)
  llhood.new <- llhood.new + log(phi.prop)

  # Select old vs. new guess
  guessvec <- mhstep(guessvec, guessvec.new, llhood.old, llhood.new, j = 2)
  return(guessvec)

}






#' Log likelihood function for phi
#' 
#' @param guessvec list of items in MCMC
#' @param a.phi prior shape for gamma for phi
#' @param b.phi prior rate for gamma for phi
llhood.phi <- function(guessvec, a.phi, b.phi) {
  # Get guesses
  phi <- guessvec$phi

  # Get gamma llhood
  gam1 <- dgamma(phi, a.phi, b.phi, log = T)  
  # Get normal llhood
  beta1 <- llhood.beta1(guessvec)

  # Get llhood
  llhood <- gam1 + beta1
  return(llhood)
}


#' Log likelihood function for beta0
#' 
#' @param guessvec list of items in MCMC
#' @param sd.beta0 standard deviation for normal prior
llhood.beta0 <- function(guessvec, sd.beta0, quants) {
  
  # get guesses
  beta0 <- guessvec$beta0

  # beta0 ~ N(0, sd.beta0^2)
  norm1 <- dnorm(beta0, 0, sd.beta0, log = T)
  ly <- llhood.y(guessvec, quants)

  llhood <- norm1 + ly
  return(llhood)
}


#' Log likelhood function for beta1 prior
#'
#' @param guessvec list of items in MCMC
llhood.beta1 <- function(guessvec) {
  
  # Get guesses
  beta1 <- guessvec$beta1
  Dists <- guessvec$Dists
  sigma2 <- guessvec$sigma2
  phi <- guessvec$phi
  theta <- guessvec$theta

  # Get Sigma covariance
  Sigma <- sigma2 * exp(-(1/phi) * Dists)
  # Get normal log likelihood
  llhood <- dmvnorm((beta1), rep(theta, length(beta1)), Sigma, log = T)
  return(llhood)
}


#' Log likelhood function for beta1 total
#'
#' @param guessvec list of items in MCMC
llhood.beta1.out <- function(guessvec, quants) {

  # Get normal likelihood of beta1
  lbeta1 <- llhood.beta1(guessvec)
  # Get likelihood for y | beta1
  ly <- llhood.y(guessvec, quants)

  # Get log likelihood
  llhood <- lbeta1 + ly
  return(llhood)
}


#' Function to get likelihood of y
#'
#' @param guessvec list of items in MCMC
llhood.y <- function(guessvec, quants) {
  # Get guesses
  y <- guessvec$y
  beta0 <- guessvec$beta0
  beta1 <- guessvec$beta1
  x <- guessvec$x

  # normal or poisson?
  type <- guessvec$type

  # Get mean of  distribution
  beta1b <- sweep(x, 2, beta1, "*")
  beta1b <-  rowSums(beta1b) * 1/length(beta1)
  #beta1b <- apply(beta1b, 1, function(x) auc(quants, x))
  
  if(type == "norm") {
    mu <- beta0 + beta1b
    llhood <- sum(dnorm(y, mu, sd = 0.01, log = T))
  }else if(type == "pois") {
    mu <- exp(beta0 + beta1b)

    # Get log likelihood
    llhood <- sum(dpois(y, mu, log = T)) 
  
  }  
  return(llhood)
}





#' Function to choose new vs. old

#' @param guessvec list of items in MCMC
#' @param guessvec.new list of items in MCMC with updated guess
#' @param llhood.old log likelihood for old guess
#' @param llhood.new log likelihood for new guess
#' @param j 1 = beta0, 2 = beta1, 3 = phi
mhstep <- function(guessvec, guessvec.new, llhood.old, llhood.new, j) {
  # Find ratio
  rat <- llhood.new - llhood.old

  # Get compare
  alpha <- log(runif(1))

  # Accept with probability alpha
  accept <- ifelse(alpha < rat, 1, 0)

  if(is.null(accept) || is.na(accept) || is.infinite(accept)) {
    accept <- 0
  }

  # If accept, update guessvec
  if(accept) {
    guessvec <- guessvec.new
  }

  # Update acceptance sum
  guessvec$accept[j] <- guessvec$accept[j] + accept

  return(guessvec)
}


