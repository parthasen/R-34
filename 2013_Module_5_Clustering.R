# ================================================== #
# Canadian Bioinformatics Workshops Series           #
# Toronto, May 23 and 24 2013                        #
# Exploratory Analysis of Biological Data using R    #
#                                                    #
# Faculty:                                           #
#   Boris Steipe <boris.steipe@utoronto.ca>          #
#                                                    #
# Contributions by:                                  #
#   Raphael Gottardo, FHCRC                          #
#   Sohrab Shah, UBC                                 #
#                                                    #
#                                                    #
# Module 5: Clustering                               #
#                                                    #
# ================================================== #

# ==================================================
# Hierarchical clustering
# ==================================================

# select the expression data of the first 50 cell cycle data genes
# (We choose 50, only so that the plots don't get too dense)


cho.data<-as.matrix(read.table("./logcho_237_4class.txt",skip=1)[1:50,3:19])
D.cho<-dist(cho.data, method = "euclidean")

# single-linkage hierarchical clustering
hc.single<-hclust(D.cho, method = "single", members=NULL)

plot(hc.single)

# draw rectangles at different cut-levels
rect.hclust(hc.single,k=2)
rect.hclust(hc.single,k=3)
rect.hclust(hc.single,k=4)
rect.hclust(hc.single,k=5)
rect.hclust(hc.single,k=25)

# now retrieve the actual indices and use them to generate
# paralell coordinate plots

class.single<-cutree(hc.single, k = 4)

# Explain the output
class.single

oPar <- par(mfrow=c(2,2))
matplot(t(cho.data[class.single==1,]),type="l", xlab="time",ylab="log expression value")
matplot(t(cho.data[class.single==2,]),type="l", xlab="time",ylab="log expression value")
matplot(as.matrix(cho.data[class.single==3,]),type="l", xlab="time",ylab="log expression value") # use as.matrix(), R apparently does not distinguish row and column vectors otherwise
matplot(t(cho.data[class.single==4,]),type="l", xlab="time",ylab="log expression value")
par(oPar)



# Complete linkage clustering - same operations
# Complete linkage aims to find similar clusters
hc.complete<-hclust(D.cho, method = "complete", members=NULL)
plot(hc.complete)
rect.hclust(hc.complete,k=4)
class.complete<-cutree(hc.complete, k = 4)

rect.hclust(hc.complete,k=4)

class.complete<-cutree(hc.complete, k = 4)

oPar <- par(mfrow=c(2,2))
matplot(t(cho.data[class.complete==1,]),type="l", xlab="time",ylab="log expression value")
matplot(t(cho.data[class.complete==2,]),type="l", xlab="time",ylab="log expression value")
matplot(t(cho.data[class.complete==3,]),type="l", xlab="time",ylab="log expression value") 
matplot(t(cho.data[class.complete==4,]),type="l", xlab="time",ylab="log expression value")
par(oPar)


# Average linkage clustering - same operations
hc.average<-hclust(D.cho, method = "average", members=NULL)
plot(hc.average)
rect.hclust(hc.average,k=4)
class.average<-cutree(hc.average, k = 4)
oPar <- par(mfrow=c(2,2), mar=c(0,0,0,0))
matplot(t(cho.data[class.average==1,]),type="l",xlab="time",ylab="log expression value")
matplot(t(cho.data[class.average==2,]),type="l",xlab="time",ylab="log expression value")
matplot(as.matrix(cho.data[class.average==3,]),type="l",xlab="time",ylab="log expression value")
matplot(t(cho.data[class.average==4,]),type="l",xlab="time",ylab="log expression value")
par(oPar)


# test 20 clusters
hc.average<-hclust(D.cho, method = "average", members=NULL)
class.average<-cutree(hc.average, k = 20)
oPar <- par(mfrow=c(4,5), mar=c(0,0,0,0))
for (i in 1:20) {
matplot(t(cho.data[class.average==i,]),type="l",xlab="time",ylab="log expression value")
}
par(oPar)



# Ward's linkage clustering - same operations
# Ward's linkage minimzes the sum-of-squares error that is incurred
# when fusing two small clusters into a larger one
hc.ward<-hclust(D.cho, method = "ward", members=NULL)
plot(hc.ward)
rect.hclust(hc.ward,k=4)
class.ward <-cutree(hc.ward, k = 4)
oPar <- par(mfrow=c(2,2))
matplot(t(cho.data[class.ward ==1,]),type="l",xlab="time",ylab="log expression value")
matplot(t(cho.data[class.ward ==2,]),type="l",xlab="time",ylab="log expression value")
matplot(t(cho.data[class.ward ==3,]),type="l",xlab="time",ylab="log expression value")
matplot(t(cho.data[class.ward ==4,]),type="l",xlab="time",ylab="log expression value")
par(oPar)



# ==================================================
# Partitioning clustering
# ==================================================

# === K-means ======================================

# random, synthetic data
set.seed(100)
x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
           matrix(rnorm(100, sd = 0.3), ncol = 2)
           )
colnames(x) <- c("x", "y")
# two seperated clusters, one with a mean of 1 and one with a mean of 0. 

# = watch the individual iterations ... cluster centres
# and class assignments are optimized

for(i in 1:4) {
  set.seed(101)
  cl<-kmeans(x, matrix(runif(10,-.5,.5),5,2),iter.max=i)
  plot(x,col=cl$cluster)
  points(cl$centers, col = 1:5, pch = 8, cex=2)
  print("Next: >")
  readline()
}

# But: be aware ...
oPar <- par(mfrow=c(2,1))
  set.seed(101)
  cl<-kmeans(x, matrix(runif(10,-.5,.5),5,2),iter.max=10)
  plot(x,col=cl$cluster)
  points(cl$centers, col = 1:5, pch = 8, cex=2)

  set.seed(102)
  cl<-kmeans(x, matrix(runif(10,-.5,.5),5,2),iter.max=10)
  plot(x,col=cl$cluster)
  points(cl$centers, col = 1:5, pch = 8, cex=2)
par(oPar)

# ... K-means does not guarantee a globally optimal solution,
# merely a locally converged one.

# === K-medoids ======================================

# load library "cluster" for K-medoid partitioning
library(cluster)
  set.seed(101)
  cl<-pam(x, 5)
  plot(x, col=cl$cluster)
  plot(cl) # shows boundary and silhouette plots

# ==================================================
# Affinity propogation clustering
# ==================================================
# Sciences, 315(5814):972-976, 2007

install.packages("apcluster")

# ==================================================
# Density estimation
# ==================================================

set.seed(103)
x1<-array(c(runif(70, 0, 10)), c(35,2))
x2<-array(c(rnorm(30, 7, 0.7)), c(15,2))
plot(x1, col="black", xlab="x", ylab="y")
points(x2, col="red")


# density estimation provides a way to analyse such a situation
# consider the example of geysir eruption intervals
# see ?density

length(faithful$eruptions)
head(faithful$eruptions, 8)
hist(faithful$eruptions, col=rgb(0.9,0.9,0.9), main="")

par(new="T")
plot(density(faithful$eruptions, bw = "sj"),
    main="", xlab="", ylab="", axes="F", col="red", lwd=3)

# ==============

# Apply this to our synthetic data

x3<-rbind(x1, x2)

plot(x1, col="black", xlab="x", ylab="y")
points(x2, col="red")
densX <- density(x3[,2], bw = "sj")
par(new=T)
plot(densX, main="", xlab="",
    ylab="", axes="F", col="blue", lwd=3)
abline(v=7, col="firebrick")

# densX$x and densX$y contain the x and dens(x) coordinates, respectively
# extract the maximum
indexMaxY <- order(densX$y, decreasing=TRUE)[1]
densX$x[indexMaxY]
abline(v=densX$x[indexMaxY], lty=2, col="blue")

# Exercise: do the same thing for the Y-axis and order the indices
# of the (x,y) pairs according to their distance to the density maximum.



# [End]





























