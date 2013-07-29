# ================================================== #
# Canadian Bioinformatics Workshops Series           #
# Toronto, May 23 and 24 2013                        #
# Exploratory Analysis of Biological Data using R    #
#                                                    #
# Faculty: Boris Steipe <boris.steipe@utoronto.ca>   #
#                                                    #
# Contributions by:                                  #
#   Raphael Gottardo, FHCRC                          #
#   Sohrab Shah, UBC                                 #
#                                                    #
#                                                    #
# Module 6: Hypothesis testing                       #
#                                                    #
# ================================================== #



# =============================================

# One sample t-test
logHiv <- log(read.table(file="hiv.raw.data.24h.txt", sep="\t", header=TRUE))
# Compute the log ratios (its a ratio of log of one expression value divided by another expression value. Division of logs is subtraction. )
logRHiv <- (logHiv[, 1:4] - logHiv[, 5:8])
# t-test a single gene
gene1 <- t.test(logRHiv[1, ], mu=0)
x <- seq(-4, 4, 0.1)
f <- dt(x, df=3)
plot(x, f, xlab="x", ylab="density", type="l", lwd=5)
segments(gene1$stat, 0, gene1$stat, dt(gene1$stat, df=3), col=3, lwd=5)
segments(-gene1$stat, 0, -gene1$stat, dt(gene1$stat, df=3), col=2, lwd=5)
segments(x[abs(x)>abs(gene1$stat)], 0, 
         x[abs(x)>abs(gene1$stat)], f[abs(x)>abs(gene1$stat)], col=4, lwd=1)

# There are an insane amount of genes that you'd like to ttest for significance. 7680 to be precise. 
# Therefore, it is important to use the script to t-test for a large number of genes. 

# =============================================

# repeat t-test for gene 4 (this is a 1-tailed t-test)
gene4 <- t.test(logRHiv[4, ], mu=0)
x <- seq(-4, 4, 0.1)
f <- dt(x, df=3)
plot(x, f, xlab="x", ylab="density", type="l", lwd=5)
segments(gene4$stat, 0, gene4$stat, dt(gene4$stat, df=3), col=3, lwd=5)
segments(-gene4$stat, 0, -gene4$stat, dt(gene4$stat, df=3), col=2, lwd=5)
segments(x[abs(x)>abs(gene4$stat)], 0, x[abs(x)>abs(gene4$stat)], f[abs(x)>abs(gene4$stat)], col=4, lwd=1)

# =============================================

# Two sample t-test (this is a two sample / or two paired t-test) - find changes amongst genes. 
data <- log(read.table(file="hiv.raw.data.24h.txt", sep="\t", header=TRUE))
gene1 <- t.test(data[1, 1:4], data[1, 5:8], var.equal = TRUE)
x <- seq(-4, 4, 0.1)
f <- dt(x, df=6)
plot(x, f, xlab="x", ylab="density", type="l", lwd=5)
segments(gene1$stat, 0, gene1$stat, dt(gene1$stat, df=6), col=3, lwd=5)
segments(-gene1$stat, 0, -gene1$stat, dt(gene1$stat, df=6), col=2, lwd=5)
segments(x[abs(x)>abs(gene1$stat)], 0, x[abs(x)>abs(gene1$stat)], f[abs(x)>abs(gene1$stat)], col=4, lwd=1)

# =============================================

# Two sample t-test
data <- log(read.table(file="hiv.raw.data.24h.txt", sep="\t", header=TRUE))
gene4 <- t.test(data[4, 1:4], data[4, 5:8], var.equal=TRUE)
x <- seq(-4, 4, 0.1)
f <- dt(x, df=6)
plot(x, f, xlab="x", ylab="density", type="l", lwd=5)
segments(gene4$stat, 0, gene4$stat, dt(gene4$stat, df=6), col=3, lwd=5)
segments(-gene4$stat, 0, -gene4$stat, dt(gene4$stat, df=6), col=2, lwd=5)
segments(x[abs(x)>abs(gene4$stat)], 0, x[abs(x)>abs(gene4$stat)], f[abs(x)>abs(gene4$stat)], col=4, lwd=1)

# =============================================

# Two sample t-test (Welch's)
data <- log(read.table(file="hiv.raw.data.24h.txt", sep="\t", header=TRUE))
gene1 <- t.test(data[1, 1:4], data[1, 5:8])
x <- seq(-4, 4, 0.1)
f <- dt(x, df=6)
plot(x, f, xlab="x", ylab="density", type="l", lwd=5)
segments(gene1$stat, 0, gene1$stat, dt(gene1$stat, df=6), col=2, lwd=5)
segments(-gene1$stat, 0, -gene1$stat, dt(gene1$stat, df=6), col=2, lwd=5)
segments(x[abs(x)>abs(gene1$stat)], 0, x[abs(x)>abs(gene1$stat)], f[abs(x)>abs(gene1$stat)], col=2, lwd=1)
gene1

# =============================================

# Paired t-test
data <- log(read.table(file="hiv.raw.data.24h.txt", sep="\t", header=TRUE))
gene1 <- t.test(as.double(data[4, 1:4]), as.double(data[4, 5:8]), paired=TRUE)
x <- seq(-4, 4, 0.1)
f <- dt(x, df=6)
plot(x, f, xlab="x", ylab="density", type="l", lwd=5)
segments(gene1$stat, 0, gene1$stat, dt(gene1$stat, df=6), col=2, lwd=5)
segments(-gene1$stat, 0, -gene1$stat, dt(gene1$stat, df=6), col=2, lwd=5)
segments(x[abs(x)>abs(gene1$stat)], 0, x[abs(x)>abs(gene1$stat)], f[abs(x)>abs(gene1$stat)], col=2, lwd=1)
gene1

# =============================================
# Wilcoxon test

# Non-parametric test principle illustrated.

# generate two random data samples with slightly
# different means
set.seed(53)
n <- 25
M <- matrix(nrow = n+n, ncol=2)
for (i in 1:n) {
	M[i,1]   <- rnorm(1, 10, 1)
	M[i,2]   <- 1
	M[i+n,1] <- rnorm(1, 11, 1)
	M[i+n,2] <- 2
}
boxplot(M[1:25,1], M[26:50,1]) # this is how different the two groups are
plot(M[,1], col=M[,2]) # we can see the trend in the values, plotted by index
plot(sample(1:(2*n)), M[,1], col=M[,2]) # however, with random indices it is quite hard to tell that these are two different distributions

# The Wilcoxon test counts for each observation of one
# group how many observations from the other group are
# ranked below it.

# To illustrate:
o <- order(M[,1])
plot(M[o,1], col=M[o,2])

wilcox.test(M[1:n,1], M[(1:n)+n,1])
#conclusion, we have two different populations. 

# Exercise
wilcox.test(as.double(data[1, 1:4]), as.double(data[1, 5:8]))
wilcox.test(as.double(data[4, 1:4]), as.double(data[4, 5:8]))



# =============================================

# Permutation test gene 1
set.seed(100)
Np <- 100
t0 <- rep(0, Np)
# Here I use the difference in means, but any statistics could be used!
t <- mean(as.double(data[1, 1:4]))-mean(as.double(data[1, 5:8]))
t
for(i in 1:Np)
{
  perm <- sample(1:8)
  newdata <- data[1, perm]
  t0[i] <- mean(as.double(newdata[1:4]))-mean(as.double(newdata[5:8]))
}
pValue <- sum(abs(t0)>=abs(t))/Np
pValue

# Permutation test gene 4
set.seed(100)
Np <- 1000
t0 <- rep(0, Np)
# Here I use the difference in means, but any statistics could be used!
t <- mean(as.double(data[4, 1:4]))-mean(as.double(data[4, 5:8]))
t
for(i in 1:Np)
{
  perm <- sample(1:8)
  newdata <- data[4, perm]
  t0[i] <- mean(as.double(newdata[1:4]))-mean(as.double(newdata[5:8]))
}
pValue <- sum(abs(t0)>=abs(t))/Np
pValue
# given the data he has, he cant distinguish between two populations - this is for gene1 
# For gene4 this means that only very rarely does a mutation looks as extreme as the one he pulled to test, this means thaat the observation is typical and intrinsic to the data. 
# Bootstrapping = sampling data again and again until to find out where the observation lies wrt to the distribution. 

# =============================================

# Bootstrap
set.seed(100)
x <- rnorm(15)
muHat <- mean(x)
sigmaHat <- sd(x)
Nrep <- 100
muHatNew <- rep(0, Nrep)
for(i in 1:Nrep)
{
  xNew <- sample(x, replace=TRUE)
  muHatNew[i] <- mean(xNew)
}
se <- sd(muHatNew)
muHat
se

set.seed(100)
x <- rnorm(15)
muHat <- mean(x)
sigmaHat <- sd(x)
Nrep <- 100
muHatNew <- rep(0, Nrep)
for(i in 1:Nrep)
{
  xNew <- sample(x, replace=TRUE)
  muHatNew[i] <- median(xNew)
}
se <- sd(muHatNew)
se

# =============================================

### Power calculation in R i.e. that we found a true positive and not a false positive (i.e. a negative result) 
# How many poeple will I need to analyze before the statistical test becomes possible (i.e. valid)

power.t.test(n = 5, delta = 1, sd=2, alternative="two.sided", type="one.sample")

# =============================================

# Power caculation
# You can skip this
# PDF version - plots in consecutive pages in a PDF file
# alpha <- 0.05
# n <- 4
# mu0 <- 0
# mu1 <- 1
# s <- 2
# name <- paste("Power_", "dmu=", mu1-mu0, "_", "s=", s, ".pdf", sep="")
# pdf(name)
# nSeq <- 2^seq(1, 8, 1)+1
# power <- rep(0, length(n))
# for(n in nSeq)
# {
#   cut <- qt(1-alpha/2, n-1)
#   x <- seq(-4, 10, 0.1)
#   f0 <- dt(x, df=n-1)
#   f1 <- dt(x, df=n-1, ncp=(mu1-mu0)/(s/sqrt(n)))
#   plot(x, f0, xlab="x", ylab="density", type="l", lwd=5, main=paste("n=", n, " dmu=", mu1-mu0, " s=", s, sep=""), ylim=c(-0.03, .5))
#   lines(x, f1, xlab="x", ylab="density", type="l", lwd=5, lty=2)
#   segments(cut, 0, cut, dt(cut, df=n-1), col=2, lwd=5)
#   segments(-cut, 0, -cut, dt(-cut, df=n-1), col=2, lwd=5)
#   segments(x[abs(x)>cut], 0, x[abs(x)>cut], f0[abs(x)>cut], col=2, lwd=1)
#   segments(x[abs(x)<cut], 0, x[abs(x)<cut], f1[abs(x)<cut], col=3, lwd=1)
#   # Power is 1-beta
#   power=1-(pt(cut, df=n-1, ncp=(mu1-mu0)/(s/sqrt(n)))-pt(-cut, df=n-1, ncp=(mu1-mu0)/(s/sqrt(n))))
#   legend(4, .3, paste("power=", round(power, digit=2), sep=""), bty="n", cex=2)
#   points(mu0, 0, pch="|", cex=1.5)
#   legend(mu0-1.4, 0.0, "mu0", bty="n", cex=10.1)
#   points(mu1, 0, pch="|", cex=1.5)
#   legend(mu1-1.2, 0.0, "mu1", bty="n", cex=10.1)
#   
# }
# dev.off()
# 

# # Power calculation - lattice plot version
# Opar <- par(no.readonly = TRUE)
# par(mfrow = c(2,4)) 
# alpha <- 0.05
# n <- 4
# mu0 <- 0
# mu1 <- 1
# s <- 2
# nSeq <- 2^seq(1, 8, 1)+1
# power <- rep(0, length(n))
# LW <- 1
# for(n in nSeq)
# {
  # cut <- qt(1-alpha/2, n-1)
  # x <- seq(-4, 10, 0.1)
  # f0 <- dt(x, df=n-1)
  # f1 <- dt(x, df=n-1, ncp=(mu1-mu0)/(s/sqrt(n)))
  # plot(x, f0, xlab="x", ylab="density", type="l", lwd=LW, main=paste("n = ", n, sep=""), ylim=c(-0.03, .5))
  # lines(x, f1, xlab="x", ylab="density", type="l", lwd=LW, lty=2)
  # segments(cut, 0, cut, dt(cut, df=n-1), col=2, lwd=LW)
  # segments(-cut, 0, -cut, dt(-cut, df=n-1), col=2, lwd=LW)
  # segments(x[abs(x)>cut], 0, x[abs(x)>cut], f0[abs(x)>cut], col=2, lwd=1)
  # segments(x[abs(x)<cut], 0, x[abs(x)<cut], f1[abs(x)<cut], col=3, lwd=1)
  # # Power is 1-beta
  # power=1-(pt(cut, df=n-1, ncp=(mu1-mu0)/(s/sqrt(n)))-pt(-cut, df=n-1, ncp=(mu1-mu0)/(s/sqrt(n))))
  # legend(4, .3, paste("power=", round(power, digit=2), sep=""), bty="n", cex=0.67)
  # points(mu0, 0, pch="|", cex=1.5)
  # points(mu1, 0, pch="|", cex=1.5)  
# }
# par(Opar)

# Power caculation
# You can skip the above

# =============================================

# Multiple testing
# If i compare all of them, how often do I find signficiance in each when i use a pairwise t-test (the answer is that they should all be insigificant, because they were taken from a normal dist. )

set.seed(100)
y <- matrix(rnorm(100000), 1000, 5)
myt.test <- function(y){
  t.test(y, alternative="two.sided")$p.value
}
P <- apply(y, 1, myt.test)
sum(P<0.05)
# All values should be false, but the prob of occuring a true is a very large possibility - you will find 42 false positive.s 

# =============================================

# FDR 
set.seed(100)
N <- 10000
alpha <- 0.05
y1 <- matrix(rnorm(9000*4, 0, 1), 9000, 4)
y2 <- matrix(rnorm(1000*4, 5, 1), 1000, 4)
y <- rbind(y1, y2)

myt.test <- function(y){
  t.test(y, alternative="two.sided")$p.value
}
P <- apply(y, 1, myt.test)
sum(P<alpha)

Psorted <- sort(P)
plot(Psorted, xlim=c(0, 1000), ylim=c(0, .01))
abline(b=alpha/N, a=0, col=2)

p <- p.adjust(P, method="bonferroni")
sum(p<0.05)
p <- p.adjust(P, method="fdr")
sum(p<0.05)

# Calculate the true FDR
sum(p[1:9000]<0.05)/sum(p<0.05)

# =============================================

# Lab on HIV data
data <- log(read.table(file="hiv.raw.data.24h.txt", sep="\t", header=TRUE))
M <- data[, 1:4]-data[, 5:8] # M = log ratios
A <- (data[, 1:4]+data[, 5:8])/2 # A= averages 

# Here I compute the mean over the four replicates
M.mean <- apply(M, 1, "mean")
A.mean <- apply(A, 1, "mean")

### t-test

n <- nrow(M)

# Basic normalization
for(i in 1:4)
  M[, i] <- M[, i]-mean(M[, i])

p.val <- rep(0, n)
for(j in 1:n)
  p.val[j] <- t.test(M[j, ], alternative = c("two.sided"))$p.value

p.val.tilde <- p.adjust(p.val, method="bonferroni")

plot(A.mean, M.mean, pch=".")
points(A.mean[p.val<0.05], M.mean[p.val<0.05], col=3)
points(A.mean[p.val.tilde<0.05], M.mean[p.val.tilde<0.05], col=2)

p.val.tilde <- p.adjust(p.val, method="fdr")
# picks the best data-points based on changes in a p-value bar set at 10% 
# 

plot(A.mean, M.mean, pch=".")
points(A.mean[p.val<0.05], M.mean[p.val<0.05], col=3)
points(A.mean[p.val.tilde<0.1], M.mean[p.val.tilde<0.1], col=2)

M.sd <- apply(M, 1, "sd")

plot(M.mean, M.sd)
points(M.mean[p.val.tilde<0.1], M.sd[p.val.tilde<0.1], col=2)

# =============================================

library(samr)
y <- c(1, 1, 1, 1)
data.list <- list(x=M, y=y, geneid=as.character(1:nrow(M)), genenames=paste("g", as.character(1:nrow(M)), sep=""), logged2=TRUE)
res <- samr(data.list, resp.type=c("One class"), s0=NULL, s0.perc=NULL, nperms=100, center.arrays=FALSE, testStatistic=c("standard"))

delta.table <- samr.compute.delta.table(res, nvals=500)

kk <- 1
while(delta.table[kk, 5]>0.1)
  kk <- kk+1

delta <- delta.table[kk, 1]
siggenes.table <- samr.compute.siggenes.table(res, delta, data.list, delta.table) 
ind.diff <- sort(as.integer(c(siggenes.table$genes.up[, 3], siggenes.table$genes.lo[, 3])))
n1 <- dim(data)[1]
ind.log <- rep(FALSE, n1)
ind.log[ind.diff] <- TRUE

plot(A.mean, M.mean, pch=".")
points(A.mean[ind.log], M.mean[ind.log], col=3)

# =============================================

### Limma ###
source("http://bioconductor.org/biocLite.R")
biocLite("limma")
## Compute necessary parameters (mean and standard deviations)
library(limma)
fit <- lmFit(M)
## Regularize the t-test
fit <- eBayes(fit)
## Default adjustment is BH (FDR)
topTable(fit, p.value=0.05, number=100)

# [End]
