
# File to see what beta functions will look like for different associations

# Function to get outcome
gety <- function(beta1, xvar1, sd1 = 0.01, ...) {
  # get additive error
  eps <- rnorm(ncol(xvar1), sd = sd1)
  # compute y
  y1 <- colSums(sweep(xvar1, 1, beta1, "*")) + eps
  y1
}


getx <- function(ny, argvals, N = 100, sd2 = 1, ...) {
  na <- length(argvals)
  xvar1 <- matrix(nrow = na, ncol = ny)
  seq1 <- sample(c(0, 1), ny, replace = T)
  mns <- 0 + seq1 * 1

  for(i in 1 : ny) {
    xvars <- (rnorm(N, mean = mns[i], sd = sd2))
    xvar1[, i] <- quantile(xvars, probs = argvals) 
  }
  xvar1
}

getreg <- function(xvar1, y1, argvals, na, ns1 = 5, ...) {  
  # Create basis
  basis1 <- create.bspline.basis(c(0, 1), ns1)
  # Get smooth of x
  xfn <- smooth.basis(argvals1, xvar1, basis1)$fd  
  freg1 <- fRegress(y1 ~ xfn)
  freg1
}


simx <- function(ny, beta1, na = 10, ...) {
  argvals <- seq(0.1, 0.9, length = na)
  xvar1 <- getx(ny, ...)
  y1 <- gety(beta1, xvar1, ...)
  f1 <- getreg(xvar1, y1, na, ...)
  list(xvar1 = xvar1, y1 = y1, f1 = f1)
}






library(dplyr)
library(tidyr)
library(fda)




ny <- 100
na <- 20
beta1 <- rep(c(0, 1), each = na / 2)
#beta1 <- 0.5 - dnorm(seq(0, 1, length = na), mean = 0.5)

#f1 <- simx(ny, beta1, ns1 = 10, sd1 = 0.01)




# Specify argument values (%iles)
argvals1 <- seq(0, 1, length = na)


xvar1 <- getx(ny, argvals = argvals1)
y1 <- gety(beta1, xvar1, sd1 = 1)


# Create basis
ns1 <- 5
#basis1 <- create.polygonal.basis(c(min(argvals1), max(argvals1)), argvals = argvals1)

basis1 <- create.bspline.basis(c(0, 1), norder = ns1)

# Get smooth of x
xfn <- smooth.basis(argvals1, xvar1, basis1)$fd

# The intercept must be constant for a scalar response
betabasis1 <- create.constant.basis(c(0, 1))
betafd1    <- fd(0, betabasis1)
betafdPar1 <- fdPar(betafd1)

betafd2     <- with(xfn, fd(basisobj=basis, fdnames=fdnames))
# convert to an fdPar object
betafdPar2  <- fdPar(betafd2)
betalist <- list(const=betafdPar1, tempfd=betafdPar2)


xfdlist <- list(const=rep(1, ny), xfn=xfn)

fd1 <- fRegress(y1, xfdlist, betalist)



#fd1 <- f1$f1
plot1 <- fd1$betaestlist[[2]]
plot(plot1$fd)
#std1 <- fRegress.stderr(freg1, diag(1, nrow(y1)), )


# try standard regression
txvar1 <- t(xvar1)
y1 <- f1$y1
colnames(txvar1) <- paste0("x", seq(1, ncol(txvar1)))
eqn1 <- paste0("y1 ~", paste(colnames(txvar1), collapse = "+")) 
dat <- data.frame(y1, txvar1)
s1 <- summary(lm(eval(eqn1), data = dat))$coef %>% round(., 2)

print(s1)
plot(s1[, 1] )

