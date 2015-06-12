# File to test mcmc_corr.R
library(mvtnorm)
library(truncnorm)
library(MESS)

source("mcmc_corr.R")

# CHECK! Works for nothing
ny <- 99
na <- 10
set.seed(10)
x1 <- matrix(rtruncnorm(ny * na, mean = 15, a = 0, sd = 1.5), nrow = ny)
quants1 <- c(0.5, 0.9)
quants1 <- seq(0.01, 0.99, length = 20)
x <- apply(x1, 1, quantile, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
x <- apply(x1, 1, quantile, probs = quants1)
x <- t(x)
lmu <- 0


# CHECK!  Works for beta0
lmu <- 5 


# Try for beta1
# Issues:
# 1. Low acceptance probability (incorporate correlation?)
beta1 <- rep(0.1, ncol(x))
lmu <- rowSums(sweep(x, 2, beta1, "*"))
lmu <- sweep(x, 2, beta1, "*")
lmu <- apply(lmu, 1, function(x) auc(quants1, x))
lmu <- lmu + 3

#phi <- 0
#beta0 <- 0
#beta1 <- rep(10, 10)
#sigma2 <- 1
#guessvec1 <- list(beta0 = beta0, beta1 = beta1, phi = phi, sigma2 = sigma2)




set.seed(10)
y <- rpois(ny, exp(lmu))
m1 <- mcmcout(y, x, quants = quants1, niter = 250000, burnin = 0, thin = 1)
save(m1, file = "test_mcmc.RData")


m1$accept
apply(m1$beta1, 2, mean)


plot(m1$beta0)
plot(m1$sigma2)
plot(m1$phi)

plot(m1$beta1[, 1])
plot(m1$beta1[, 2])
