# Function to plot information about x

plotx <- function(x, argvals1) {
  # Specify quantiles of x, and all x
  x1 <- x$x1
  xall <- x$xall
  colnames(xall) <- paste0("x", seq(1, ncol(xall)))
  colnames(x1) <- colnames(xall)  
  # Get in form for ggplot
  mxall <- melt(xall)

  # Get in form for ggplot
  d1 <- data.frame(argvals1, x1)
  mx <- melt(d1, id.vars = "argvals1")
  mx$variable <- tolower(mx$variable)

  # Select one day to highlight
  samp1 <- paste0("x", sample(seq(1, ncol(xall)), 1))
  mx <- mutate(mx, select1 = ifelse(variable == samp1, "Yes", "No")) 
  mxall <- mutate(mxall, select1 = ifelse(Var2 == samp1, "Yes", "No"))

  # Order so plot highlighted last
  lev1 <- mxall$Var2[order(mxall$select1)]
  mxall$Var2 <- factor(mxall$Var2, levels = unique(lev1))  

  lev1 <- mx$variable[order(mx$select1)]
  mx$variable <- factor(mx$variable, levels = unique(lev1))  
  
  # Create density plots of all x's
  g1 <- ggplot(data = mxall, aes(x = value, group = Var2, colour = select1)) + 
    geom_density(aes(size = select1)) + theme_bw() +
    scale_colour_manual(guide = F, values = c("grey80", "red")) +
    scale_size_manual(guide = F, values = c(1, 1.4)) +
    theme(legend.position = "none") +
    xlab("Concentration") + ggtitle("Density plots of personal exposure")

  # Create quantile plots of x quantiles
  g2 <- ggplot(data = mx, aes(x = argvals1, y = value, group = variable, colour = select1)) + 
    geom_line(aes(size = select1)) + theme_bw() +

    scale_colour_manual(guide = F, values = c("grey80", "red")) +
    scale_size_manual(guide = F, values = c(1, 1.4)) +
	
    theme(legend.position = "none") + xlab("Quantile") +
    ylab("Concentration") + ggtitle("Plot of quantiles of concentration")

  # Arrange output for plots 
  grid.arrange(g1, g2, ncol = 2)
}


# Function to plot functional x
plotxfn <- function(x, argvals1) {
  # Side by side plots
  par(mfrow = c(1, 2))

  # Specify x quantiles and functional x
  x1 <- x$x1
  xfn <- x$xfn

  # Plot quantiles of x
  plot(argvals1, x1[, 1], type = "l", xlim = c(0, 1), ylim = c(min(x1), max(x1)), 
    xlab = "Quantiles", ylab = "Concentration",
    main = "Sample quantiles")
  for(i in 2 : ncol(x1)) {
    points(argvals1, x1[, i], type = "l")
  }

  # Plot functional x
  p1 <- plot(xfn, col = 1, lty = 1, ylim = c(min(x1), max(x1)), xlab = "Quantiles", ylab = "Concentration",
    main = "Smoothed quantiles using b-splines")
}


# Function to plot betas
plotbeta1 <- function(sim1, argvals1, cols = NULL, disttype = "norm", main1 = NULL, argvalslr = argvals1) {
  if(is.null(cols)) {
    cols <- brewer.pal(4, "Dark2")[c(2, 1, 3, 4)]
    cols[2] <- "black"
    cols[3 : 4] <- brewer.pal(5, "Blues")[c(5, 3)]

  }



  # Specify traditional regression beta
  beta2 <- sim1$beta2[, 1]
  beta3 <- sim1$beta3[, 1]
  # Fix multivariate beta
  argvals2 <- as.numeric(names(beta3))
  
  # Get beta function
  betaf <- sim1$betaf

  if(disttype == "norm") {
    #lmf <- sim1$fmod1
    # Specify functional beta
    #fd1 <- lmf$freg

    # Get std for fun beta
    #betafstd <- lmf$betafstd
    #betafstd <- betafstd$betastderrlist[[2]]

    # Specify functional beta and plot
    #betaest <- fd1$betaestlist[[2]]
    #betaest <- betaest$fd
    #p1 <- plotbeta(betaest, betafstd, argvals1) 
    
    plot(1, 1, type = "n", ylim = c(-.1, .3), xlim = c(0, 1),
      xlab = "Quantile", ylab = "Beta", 
      main = "Beta function", col = cols[1])
  } else if(disttype == "pois") {
    basis1 <- sim1$basis1
    #lmf <- sim1$fmod1
    #coef1 <- summary(lmf)$coef[-c(1,2), ]
    #fd1 <- fd(coef1[, 1], basis1)
    #fdlb <- fd(coef1[, 1] - 1.96 * coef1[, 2], basis1)
    #fdub <- fd(coef1[, 1] + 1.96 * coef1[, 2], basis1)
    
    if(is.null(main1)) {
      main1 <- "Beta function"
    }
   cols[1] <- "white" 
    #p1 <- plot(fd1, xlab = "Quantile", ylab = "Beta",
    #  ylim = c(-0.5, 1), main = main1, col = cols[1])
   plot(.5, .5, xlim = c(0, 1), type = "n", xlab = "Quantile", ylab = "Beta", ylim = c(-0.5, 1), main = main1) 
   #plot(fdlb, add = T, col = cols[1], lty = 2)
    #plot(fdub, add = T, col = cols[1], lty = 2)

    #p1 <- plot(lmf$beta.l$x, xlab = "Quantile", ylab = "Beta",
     # main = "Beta function", col = cols[1])

    #matplot(argvals1, cbind(lmf$BetaHat[[1]], lmf$Bounds[[1]]), type = "l", 
     # col = cols[1], xlab = "Quantile", ylab = "Beta",
     # main = "Beta function", ylim = c(-1, 5 ), xlim = c(0, 1))
  }



  # Add true betas
  if(class(betaf) == "function") {
   
    beta1 <- betaf(argvals1)
    points(argvals1, beta1, type = "l", col = cols[2])
  }else{
    nums <- as.numeric(names(betaf)) / 100
    points(nums, betaf, col = cols[2], pch = 1, cex = 2.5)
  }

  
  # Add points for traditional regression betas
  points(argvals2 + 0.01, beta3, col = cols[3], pch = 8)
  lbs <- beta3 - sim1$beta3[, 2] * 1.96
  ubs <- beta3 + sim1$beta3[, 2] * 1.96
  segments(y0 = lbs, x0 = argvals2 + 0.01, y1 = ubs, col = cols[3], lty = 1)

  
  points(argvalslr - 0.01, beta2, col = cols[4], pch = 16)  
  lbs <- beta2 - sim1$beta2[, 2] * 1.96
  ubs <- beta2 + sim1$beta2[, 2] * 1.96
  segments(y0 = lbs, x0 = argvalslr - 0.01, y1 = ubs, col = cols[4], lty = 1)

  abline(h = 0, lty = 2, col = "grey70")  
  # Add legend
  legend("top", legend = c("Truth", "FunDat", "mvLM", "uvLM"), 
    col = cols[c(2, 1, 3, 4)], lty = 1, pch = 16)
}


# Plot all x (data + smoothed)
plotallx <- function(x, argvals2) {
  # Plot density plots and quantiles	  
  plotx(x, argvals2)
  # Plot quantiles and smoothed quantiles
  plotxfn(x, argvals2)
}

