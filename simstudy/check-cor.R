
load("true-sim-poisson.RData")

i <- 1
rs1 <- runsim1$simout1[[i]]



# Does subtracting out the mean/median help?
x1 <- t(rs1$x1$x1)
summary(as.vector(cor(x1)))

#try subtract out median
wh1 <- which(colnames(x1) == "50%")
xm <- x1[,wh1 ]
xm <- apply(x1, 1, mean)

xs <- sweep(x1, 1, xm)
summary(as.vector(cor(xs[, -wh1])))




# what about rescaling?
xall <- t(rs1$x1$xall)
xall <- data.frame(xall, stringsAsFactors = F)
xall[, 2] <- as.numeric(xall[, 2])
mn <- tapply(xall[, 2], xall[, 1], mean)

sd1 <- tapply(xall[, 2], xall[, 1], sd)
