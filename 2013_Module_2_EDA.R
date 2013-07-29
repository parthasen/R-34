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
# Module 2: EDA (Exploratory Data Analysis)          #
#                                                    #
# ================================================== #


# ==================================================
# Digression: Generating synthetic data for analysis
# ==================================================

# Reminder: sequences of numbers in R
0:10 
seq(0, pi, 5*pi/180) 
rep(1:3, each=3, times=2)
for (i in 1:10) { print(i*i) }

# Synthetic data: usually a known function plus noise.
x <- seq(0.1, 1, 0.05)

## Function...
y <- 2*x - 0.1  #linear
plot(x,y)

y <- sin(2* pi*x) - 0.1  # sinusoidal
plot(x,y)

y <- exp(-4*x)  # exponential decay
plot(x,y)

y <- 1/(1+exp(-10 * (x-0.5)))  # Sigmoid (logistic) function
plot(x,y)

y <- dnorm(6*x - 3)   # normal (Gaussian) distribution
plot(x,y)

## "noise"...

# uniform
y <- runif(21)
plot(x,y)

y <- rnorm(21)
plot(x,y)

x <- seq(0.1, 1, 0.002)
y <- rnorm(length(x))
plot(x,y)

# function + noise
y <- rnorm(length(x), 0, 0.1) + 1/(1+exp(-10 * (x-0.5 + rnorm(length(x), 0, 0.01))))
plot(x,y)



# Random sequences
# A: uniformly distributed elements
rNucUnif <- function(n) {
	nuc <- c("A", "C", "G", "T")
	rSeq <- sample(nuc, n, replace=TRUE)
	return(rSeq[1])
}

# B: arbitrary target frequencies
rNucArb <- function(n, char, freq) { #char = nuc = c("A", "C", "G", "T")
	l <- length(char) #length of nuc , which is 4. 
	if (length(freq) < l) { #if length of freq = f = c(20,15,10,40) < 
		stop("Panic: insufficient number of target frequencies.")
	}
	freq <- freq[1:l]      # truncate vector if it is longer than l
	freq <- cumsum(freq)   # convert to interval boundaries
	freq <- freq / freq[l]    # normalize by dividing by the last element
	rSeq <- rep(0,n) # initialize the ouput vector
	for (i in 1:n) {
		x <- runif(1)
		# we use the fact that logical values are converted
		# into 0's and 1's to count the number of intervals
		# that are smaller than x.
		index <- max(cumsum( freq < x )) + 1 
		rSeq[i] <- char[index] # assign the appropriate nucleotide
	}
	return(rSeq)
}

# Try:
nuc <- c("A", "C", "G", "T") # define an "alphabet"
f <- c(20,15,10,40) # define target frequencies

seq <- rNucArb(200, nuc, f) # invoke the function

seq
paste(seq, collapse="")  # write the sequence as a single string
table(seq) # tabulate character frequencies
sort(table(seq)) # sort
barplot(sort(table(seq), de=TRUE)) # sort descending for the barplot

# ==================================================
# Permutations
# ==================================================

sums <- rep(0, 10)
digits <- 0:9
rsd <- vector()
for (i in 1:10000) {
	perm <- sample(digits, 10)
	sums <- sums + perm
	if (i %% 1000 == 0) {
		rsd <- c(rsd, sd(sums)/mean(sums))
	}
}

names(sums) <- 0:9
sums
barplot(sums)
plot(rsd)



# ==================================================
# Probability distributions
# ==================================================

# The binomial distribution
?dbinom
x <- 0:1
f <- dbinom(x, size=1, 0.1)
plot(x, f, xlab="x", ylab="density", type="h", lwd=5)

set.seed(100)
x <- rbinom(100, size=1, prob=.1)
hist(x)

# ==================================================
# The Normal distribution
# ==================================================
?dnorm
x <- seq(-4, 4, 0.1)
f <- dnorm(x, mean=0, sd=1)
plot(x, f, xlab="x", ylab="density", lwd=5, type="l")

# ==================================================
# Explore lines (Section 3 of PlottingReference.R)
# ==================================================


# =============================================
# Overlay a histogram with a line
set.seed(100)
x <- rnorm(100, mean=0, sd=1)
hist(x)
lines(seq(-3,3,0.1),50*dnorm(seq(-3,3,0.1)), col="red")
#lines is great to plot explanatory curves or a line into a graph. 
#lines()
#points()
#text()

# =============================================
# quantiles

q90 <- qnorm(0.90, mean = 0, sd = 1)
x <- seq(-4, 4, 0.1)
f <- dnorm(x, mean=0, sd=1)
plot(x, f, xlab="x", ylab="density", type="l", lwd=5)
abline(v=q90, col=2, lwd=5)

# =============================================
# empirical quantiles

set.seed(100)
x <- rnorm(100, mean=0, sd=1)
quantile(x)
quantile(x, probs=c(.1, .2, .9))

# =============================================

set.seed(100)
x <- rnorm(100, mean=0, sd=1)
mean(x)
median(x)
IQR(x)
var(x)
sd(x)
summary(x)

# =============================================

set.seed(100)
x <- rnorm(100, mean=0, sd=1)
boxplot(x)

x <- rnorm(100, mean=-1.5, sd=1)
x <- c(x, rnorm(100, mean=1.5, sd=1))
boxplot(x)

# Violin plot
install.packages("ggplot2")
library(ggplot2)
X <- as.data.frame(x)

p <- ggplot(X, aes(1,x))
p + geom_violin()
# See ggplot2 introductory tutorial at 
#     http://www.r-bloggers.com/basic-introduction-to-ggplot2/
# and violin plot documentation at
#     http://docs.ggplot2.org/current/geom_violin.html

# =============================================

set.seed(100)
x <- rnorm(100, mean=0, sd=1)
qqnorm(x)
qqline(x, col=2)

# =============================================
# Plotting lines and legends
# Example: compare the normal distribution with
# the t-distribution
x <- seq(-4, 4, 0.1)
f1 <- dnorm(x, mean=0, sd=1)
f2 <- dt(x, df=2)
plot(x, f1, xlab="x", ylab="density", lwd=5, type="l")
lines(x, f2, lwd=5, col=2)
legend(-4, .4, c("Normal", "t2"), col=1:2, lwd=5)

# =============================================
# use qqplot to compare normally distributed samples with
# t-distributed samples

set.seed(100)
t <- rt(100, df=2)
qqnorm(t)
qqline(t, col=2)

# =============================================
# Programming example: 
# Verify the Central Limit Theorem
# Generate random deviates: numbers that are the
# sum of many small uniformly distributed deviations.
# Contrary to naive expectations, these don't in fact
# all average out to their mean.

set.seed(101)
generateVariates <- function(n) {
	Nvar <- 10000    # number of small random numbers to add
	Vout <- c()      # initialize a vector
	for (i in 1:n) { 
		x <- runif(Nvar, -0.01, 0.01) # fill a vector with many small uniformly distributed random numbers
		Vout <- c(Vout, sum(x) ) # sum over the random numbers and append to vector
	}
	return(Vout)
}

x <- generateVariates(1000)
y <- rnorm(1000, mean=0, sd=1)
qqnorm(x)
qqline(x, col=2)

# =============================================
# QQ- plot: sample against sample
set.seed(100)
x <- rt(100, df=2)
y <- rnorm(100, mean=0, sd=1)
qqplot(x, y)

# =============================================

# GvHD flow cytometry data
gvhd <- read.table("GvHD.txt", header=TRUE)
# Only extract the CD3 positive cells
hist(gvhd[,5])
gvhdCD3p <- as.data.frame(gvhd[gvhd[, 5]>280, 3:6])
plot(gvhdCD3p[, 1:2])

# =============================================
# Explore scatter plots
# =============================================
# Scatter plots are extremely useful, but learning
# a bit about R's plotting capabilities will help
# tremendously to create INFORMATIVE plots.

# Topics:
# 6 - Plotting symbols and characters
# 2 - Colors
# 4 - Coordinates

# Some special packages
# The overlap in the GvHD data can obscure
# data relationships. Here are some alternatives
# for dense plots:

# load "hexbin" package from CRAN
install.packages("hexbin")
library("hexbin")
hb <- hexbin(gvhdCD3p[, 1:2], xbins = 20)
plot(hb, colramp = colorRampPalette(c("#FFFFDD", "#77EE99", "#3377AA", "#0033BB")))

# load "prada" package from BioConductor
source("http://www.bioconductor.org/biocLite.R")
biocLite("prada")

smoothScatter(gvhdCD3p[, 1:2], nrpoints=50, pch=20, cex=0.5, col="#6633BB55")

plot (gvhdCD3p[, c(1,5)], col=densCols(gvhdCD3p[, 1], gvhdCD3p[, 2]), pch=16)



# =============================================
# Trellis plots: all against all

plot(gvhdCD3p, pch=".")

# =============================================

boxplot(gvhdCD3p)

# =============================================

oPar <- par(mfrow=c(2, 2))
hist(gvhdCD3p[, 1], 50, main=names(gvhdCD3p)[1], xlab="fluorescent intensity", ylim=c(0, 120))
hist(gvhdCD3p[, 2], 50, main=names(gvhdCD3p)[2], xlab="fluorescent intensity", ylim=c(0, 120))
hist(gvhdCD3p[, 3], 50, main=names(gvhdCD3p)[3], xlab="fluorescent intensity", ylim=c(0, 120))
hist(gvhdCD3p[, 4], 50, main=names(gvhdCD3p)[4], xlab="fluorescent intensity", ylim=c(0, 120))
par(oPar)

# =============================================

data <- read.table(file="hiv.raw.data.24h.txt", sep="\t", header=TRUE)
summary(data)
boxplot(data)

# =============================================

oPar <- par(mfrow=c(2, 2))
hist(data[, 1], 50, main=names(data)[1], xlab="fluorescent intensity")
hist(data[, 2], 50, main=names(data)[2], xlab="fluorescent intensity")
hist(data[, 5], 50, main=names(data)[5], xlab="fluorescent intensity")
hist(data[, 6], 50, main=names(data)[6], xlab="fluorescent intensity")
par(oPar)

# =============================================

# 'apply' will apply the function to all rows of the data matrix
mean <- apply(data[, 1:4], 1, "mean")
sd <- apply(data[, 1:4], 1, "sd")
plot(mean, sd)
trend <- lowess(mean, sd)
lines(trend, col=2, lwd=5)

# =============================================

data <- log(read.table(file="hiv.raw.data.24h.txt", sep="\t", header=TRUE))
summary(data)
boxplot(data)

# =============================================

oPar <- par(mfrow=c(2, 2))
hist(data[, 1], 50, main=names(data)[1], xlab="fluorescent intensity")
hist(data[, 2], 50, main=names(data)[2], xlab="fluorescent intensity")
hist(data[, 5], 50, main=names(data)[5], xlab="fluorescent intensity")
hist(data[, 6], 50, main=names(data)[6], xlab="fluorescent intensity")
par(oPar)

# =============================================

# apply(M,1,f) will apply the function f to all rows of the data matrix M
mean <- apply(data[, 1:4], 1, "mean")
sd <- apply(data[, 1:4], 1, "sd")
plot(mean, sd)
trend <- lowess(mean, sd)
lines(trend, col=2, lwd=5)

# =============================================

# scatter plot
plot(data[, 1], data[, 5], xlab=names(data)[1], ylab=names(data)[5])

# =============================================

# MA plot
A <- (data[, 1]+data[, 5])/2
M <- (data[, 1]-data[, 5])
plot(A, M, xlab="A", ylab="M")

# MA plots per replicate
oPar <- par(mfrow=c(2, 2))
A1 <- (data[, 1]+data[, 5])/2
M1 <- (data[, 1]-data[, 5])
plot(A1, M1, xlab="A", ylab="M", main="rep 1")
trend <- lowess(A1, M1)
lines(trend, col=2, lwd=5)
A2 <- (data[, 2]+data[, 6])/2
M2 <- (data[, 2]-data[, 6])
plot(A2, M2, xlab="A", ylab="M", main="rep 2")
trend <- lowess(A2, M2)
lines(trend, col=2, lwd=5)
A3 <- (data[, 3]+data[, 7])/2
M3 <- (data[, 3]-data[, 7])
plot(A3, M3, xlab="A", ylab="M", main="rep 3")
trend <- lowess(A3, M3)
lines(trend, col=2, lwd=5)
A4 <- (data[, 4]+data[, 8])/2
M4 <- (data[, 4]-data[, 8])
plot(A4, M4, xlab="A", ylab="M", main="rep 4")
trend <- lowess(A4, M4)
lines(trend, col=2, lwd=5)
par(oPar)

# =============================================


# [End]
