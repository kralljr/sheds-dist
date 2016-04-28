library(rstan)
library(truncnorm)
library(mvtnorm)

# Set up data
ny <- 99
na <- 100
nquants <- 10
set.seed(10)
# previous mean =15
x1 <- matrix(rtruncnorm(ny * na, mean = 10, a = 0, sd = 1.5), nrow = ny)
means <- runif(ny, 0, 5)
x1 <- sweep(x1, 1, means, "+")
quants1 <- seq(0.01, 0.99, length = nquants)
x <- apply(x1, 1, quantile, probs = quants1)
x <- t(x)


Q <- ncol(x)
D <- expand.grid(quants1, quants1)
D <- abs(D[, 1] - D[, 2]) %>% matrix(., nrow = Q)


#seeds1 <- sample(seq(1, 10000), 1)
#print(seeds1)
seeds1 <- 8420
set.seed(seeds1)
beta1 <- rep(0.01, ncol(x))
sigma2 <- 0.0001
phi <- 10000
sigma1 <- sigma2 * exp(-1 / phi * D) 
beta1 <- rmvnorm(1, beta1, sigma = sigma1)
beta1
set.seed(15)
beta0 <- rnorm(1, 0, 10)
lmu <- beta0 + rowSums(sweep(x, 2, beta1, "*"))

set.seed(10)
#y <- rpois(ny, exp(lmu))
y <- rnorm(ny, mean = lmu, sd = 0.001)
summary(y)

truth <- list(beta0 = beta0, beta1 = beta1, Sigma = sigma1, phi = phi, tau2 = 1 / sigma2)


dat <- (x) 
Y <- y
sheds <- list(T = nrow(dat), Q = Q, Qp = Q, X = dat, Y = Y)
              #, D = D, mu0 = rep(0, Q), beta0 = beta0)

# Fit stan model  
fit <- stan(file = "sheds.stan", data = sheds, seed = 7027, iter = 5000, chains = 2, thin = 10, control = list(max_treedepth = 30))
#control = list(adapt_delta = 0.7, max_treedepth = 30))
save(fit, file = "sheds-1.RData")


# Plots

g1 <- function(pars, wchain = 1, lims) {
  dat1 <- extract(fit, pars, permuted = F) 
  iters <- dim(dat1)[1]
  nchains <- dim(dat1)[2]
  varl <- dim(dat1)[3]
  dimnames(dat1) <- list(1 : iters, 1 : nchains, 1 : varl)
  ldat1 <- melt(dat1)
  colnames(ldat1) <- c("iteration", "chain", "variable", "value")
  ldat1 <- mutate(ldat1, chain = factor(chain))

  #plot

  t1 <- truth[[pars1[i]]] %>% as.numeric()

  s1 <- dplyr::filter(ldat1, chain == wchain) %>% spread(., variable, value) %>% select(., -iteration, -chain)
  s1 <- apply(s1, 2, mean) 

  t1 <- cbind(unique(ldat1$variable), t1, s1) %>% data.frame()
  colnames(t1) <- c("variable", "truth", "postmean")
  t1 <- gather(t1, type, value, -variable)

  ggplot(ldat1) + 
    theme_bw() + ylim(lims[1], lims[2]) + 
    geom_point(alpha =0.3, aes(x = iteration, y = value, colour = chain)) +
    geom_hline(data = t1, aes(yintercept = value, color = type)) + 
    facet_wrap(~ variable)
}

pars1 <- extract(fit) %>% names()
lims1 <- rbind(c(2.2, 2.8), c(0, .03), c(0, 1))
for(i in 1 : (length(pars1) - 1)) {

  g1(pars1[i], wchain = 1, lims = lims1[i, ])

}


# issues with tau2 and phi


# Check beta1 output
getdat <- function(pars, wchain = 1) {
  for(i in 1 : length(pars)) {
    dat1 <- extract(fit, pars[i], permuted = F) 
    iters <- dim(dat1)[1]
    nchains <- dim(dat1)[2]
    varl <- dim(dat1)[3]
    dimnames(dat1) <- list(1 : iters, 1 : nchains, paste0(pars[i], 1 : varl))
    ldat1 <- melt(dat1)
    colnames(ldat1) <- c("iteration", "chain", "variable", "value")
    ldat1 <- mutate(ldat1, chain = factor(chain))
    
    
    hold <- dplyr::filter(ldat1, chain %in% wchain)  %>% 
      spread(., variable, value) 
    if(i == 1) {
      output <- hold
    } else{
      output <- full_join(output, hold)
    }
  }
  output
}
s1 <- getdat(pars1[1 : 3])

cor(s1[, -c(1, 2)])
apply(s1, 2, mean) %>% round(., 4)




