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
# Module 1: Introduction to R                        #
#                                                    #
# ================================================== #

# ==================================================
# Using scripts.

# Put the cursor in any line, or select a block of code
# then press
#    <command><return> (Mac)
#    <ctrl><r>
# to execute it in the console.

# In the console, use <up-arrow>, <left-arrow> etc. to
# retrieve and edit previous commands.

# Use source("filename") to execute an entire script at once.



# ==================================================
# set up your workspace
getwd()
setwd("/path/of/your/directory") # the argument is a string
setwd("/path/of/your/course\ directory") # "escape" blank spaces
list.files()


# =============================================
# Scalars
# ==================================================

# Executing commands

2+2
exp(-2)
pi
# there is no definition of e - but it is easy enough to obtain ...
exp(1)
sin(2*pi) # shouldn't this be zero ???
round(sin(2*pi)) # there we go.
-1/0
X # ...doesn't exist
X <- NULL
X # now it exists, even though it has no contents
TRUE
true

# ==================================================
# Help
help (pi)
?pi
# the mathematical functions ...
?Arithmetic
?abs
?log
?round
?sin
?Special

format(ISOdate(2000, 1:12, 1), "%b")


# Keyword search
help.search ("trigonometry")
??input

# ... or browse:
# From the help page of any function of a package, 
# list an index of all functions by clicking 
# "Index" at the bottom of the page.

# ... or (preferred): Google

# ==================================================
# Assigning variables
x <- 1/sqrt(4)
y <- sin(pi/6)
x+y

# ==================================================
# Strings
# ==================================================
"Hello"
x <- paste("Hello", "World") # note: " " is default separator
m <- gregexpr("(\\b\\w{2})", x, perl=T)
# Tells you what the contents of the match is - (\\b\\w{2}) great for regular expressions
# (\\b\\w{2}) - matching pattern to Hello World is \\b = boundary , \\w = word boundary. {2} - there is a sequence to this - he wants exactly two characters from the string. 
 y<-regmatches(x,m)
 #Great for patient names in a list. 
paste(y, collapse='') # This is not what we want!
paste(y[[1]], collapse='') # List syntax ... (getting ahead of ourselves here).
paste(y[[1]][1:2], collapse='') 
paste(unlist(y), collapse='') 

# ==================================================
# Dates
# ==================================================
format(ISOdate(2000, 1:12, 1), "%b")  # see strptime() for formatting options ->
format(Sys.time(), "%W")              # Today's week of the year
format(ISOdate(2013, 07, 27), "%A")   # Weekday of my birthday this year
# see:
?format
?ISOdate
?strptime

# ==================================================
# Vectors
# =============================================

Weight <- c(60,72,75,90,95,72)
Weight[1]
Weight[2]
Weight[-2]
Weight[-2:-4]
Weight[-(2:4)]
Weight
Height <- c(1.75,1.80,1.65,1.90,1.74,1.91)
BMI <- Weight/Height^2 # vector based operation
BMI

# sequences 
1:5
1.1:5
-5:5
5:0
seq(-1,1, 0.1)
seq(0,2*pi,length.out=17)
rep(1:3, 4)
rep("waltz", 4)

# =============================================
# Vector operations
x <- 1:5
x+2 # operations with a constant are applied to each element
y <- 6:2
x+y # operations with two vectors are applied to each matching
    # pair. If the vectors have different lenghts, the shorter
    # one is recycled.
x <- rep(1,13)
y <- 0:2
x+y


# =============================================
# Vector types
x <- c(1, 1, 2, 3, 5, 8, 13)     # Numeric
x

x <- c(TRUE, TRUE, FALSE, TRUE)  # Logical
x

x <- c ("Hello","world")         # Character
x 

x <- c(1, TRUE, "Thursday")      # Mixed - R silently castes objects into different type. 
x
is.vector(x) # Remember: all elements in a vector 
             # must have the same type!



# =============================================

Weight[5] <- NA
mean(Weight)
mean(Weight,na.rm=TRUE)


# =============================================
# Matrices
# =============================================

x<-1:12
x
length(x)
dim(x)
dim(x)<-c(3,4) #(Row, Column)
x
a<-matrix(1:12,nrow=3,byrow=TRUE)
a
a<-matrix(1:12,nrow=3,byrow=FALSE)
a
rownames(a)<-c("A","B","C")
a
colnames(a)<-c("1","2","x","y")
a

# =============================================

x1 <- 1:4  # Define three vectors  
x2 <- 5:8 
y1 <- c(3,9) 
MyMatrix <- rbind(x1,x2) 
MyMatrix 
MyNewMatrix <- cbind(MyMatrix,y1) 
MyNewMatrix

# =============================================
# Factors

Pain <- c(0,3,2,2,1)
SevPain <- as.factor(c(0,3,2,2,1)) # or: as.factor(Pain)
SevPain
levels(SevPain) <- c("none","mild","medium","severe")
SevPain
is.factor(SevPain)
is.vector(SevPain)


# =============================================
# Lists
# =============================================

A<-c(31,32,40)
S<-as.factor(c("F","M","M","F"))
L<-c("London","School")
MyFriends<-list(age=A,sex=S,meta=L)
MyFriends
MyFriends$age


# =============================================
# Data frames
# =============================================

Probands <- data.frame(age=c(31,32,40,50), sex=S)
Probands
Probands[[1]]
Probands[["age"]]
Probands$age

# =============================================
# Names of objects - in vectors, lists, ...
# =============================================

x <- 1:3
names(x)
names(x) <- c("a", "b", "c")
x
names(Probands)
names(Probands) <- c("age", "gender") 
names(Probands)[1] <- c("Age")

# =============================================
# Indexing
# =============================================

?Extract

# Indexing a vector
Pain<-c(0,3,2,2,1)
Pain[1]
Pain[2]
Pain[1:2]
Pain[c(1,3)]
Pain[-5]

# Indexing a matrix
MyNewMatrix[1,1]
MyNewMatrix[1,] #first row
MyNewMatrix[,1] #first column
MyNewMatrix[,-2]

# Indexing a list
MyFriends[3]
MyFriends[[3]]
MyFriends[[3]][1]
# Indexing a data frame
Probands[1,]     # first row
Probands[2,]     # second row
Probands[,1]     # first column
Probands[[1]]    # first object in list: a vector # give object in a list or data frame. 
Probands["Age"]  # first column by name
Probands$Age     # first object by name

# =============================================

Pain; SevPain
Age <- c(45,51,45,32,90)
Pain[SevPain=="medium" | SevPain=="severe"]
Pain[Age>32] #extract data from various places into a single source. 
which(SevPain == "medium" | SevPain == "severe")

# =============================================
# Reading files
# =============================================

GvHD <-read.table("GvHD+.txt", header=TRUE) 
GvHD[1:10,] #flow-cytometry data. 

# there are many other "read" functions.
# refer to the R Data Import / Export manual
# http://cran.r-project.org/doc/manuals/R-data.html
# see:
?read.table # includes read.csv and read.delim
?read.fwf

# Excel spreadsheets should be converted to csv
# or tsv format. Alternatively the package
# xlsreadwrite is available via CRAN
# http://cran.r-project.org/web/packages/xlsReadWrite/

# =============================================
# Libraries
# =============================================

# The "survival" library is part of the standard R distribution
library(survival)
# Example
leukemia.surv <- survfit(Surv(time, status) ~ x, data = aml) 
plot(leukemia.surv, col=c(1,2), lty = 2:3) 
legend(100, .9, c("Maintenance", "No Maintenance"), col=1:2, lty = 2:3) 
title("Kaplan-Meier Curves\nfor AML Maintenance Study") 

# SAM (Significance Analysis for Microarrays) is a statistical
# analysis package from the Tibshirani lab at Stanford.
# It is not part of the normal R distribution.
library(samr)

install.packages("samr")
library(samr) # "install" is not enough - you must also LOAD the library
?samr

# =============================================
# Programming R
# =============================================

# Conditional statements
x <- -2
if (x>0) {
  print(x)
} else { 
  print(-x)
}

if(x>0) {
print(x) 
}else if (x==0) {
  print(0)
} else {
  print(-x)
}

# Logical vectors
?TRUE

# Explore conditional statements
if (TRUE)   {print("true")} else {print("false")}
if ("true") {print("true")} else {print("false")}
if ("t")    {print("true")} else {print("false")}
if (1)      {print("true")} else {print("false")} # can be used as internal true and false boolean.
if (0)      {print("true")} else {print("false")}
if (-1)     {print("true")} else {print("false")}
if (pi)     {print("true")} else {print("false")}
if (Inf)    {print("true")} else {print("false")}
if (NULL)   {print("true")} else {print("false")}
if (NA)     {print("true")} else {print("false")}
if (NaN)    {print("true")} else {print("false")}
if (Inf)    {print("true")} else {print("false")}

# Logical operators
TRUE
! TRUE
FALSE > TRUE
FALSE < TRUE
FALSE < -1
0 == FALSE
"x" == "u"
"x" >= "u"
"x" <= "u"
"x" != "u"
TRUE | FALSE
TRUE & FALSE

# equality and identity
?identical
a <- c(TRUE)
b <- c(FALSE)
a; b
a == b
identical(a, b)

b <- c(b, FALSE)
a; b
a == b # ! Ooops - no longer a single value!
identical(a, b)


# some other useful tests for conditional expressions
?all
?any
?duplicated
?exists
?is.character
?is.factor
?is.integer
?is.null
?is.numeric
?is.unsorted
?is.vector


# =============================================
# Loops

# for loop
n <- 1000000
x <- rnorm(n,10,1) #mean of 10, STD of 1
y <- sqrt(x)

# vectors can be extended in length as necessary ...
l <- c()
for (i in 1:100) {  # try this with successively larger values
  l[i] <- sqrt (x[i])
}

# ... at a great performance cost. Precreating them for the
# required length is much more efficient. 
l <- rep(0, n) # Try not to build complex data structures, build your loop with a beginning and end structure. 
for (i in 1:10000) {  # try this with successively larger values
  l[i] <- sqrt (x[i])
}


# while loop
Counter <- 1
while(Counter <= 100000) {
  y[Counter] <- sqrt(x[Counter])
  Counter <- Counter+1
}

mf <- function(a){
	for (i in 1:a) {
		print("Hello")
	}
}

ma <- mean (a){
	for (i in 1:a){
		print("Hello")
	}
}

# =============================================
# Explore functions

Oracle <- function() {
	WiseWords <- c(
		"Joy",
		"Plan",
		"Disappear",
		"Certainly",
		"Perhaps",
		"Past",
		"Respect",
		"Worry",
		"Struggle",
		"Play",
		"Regret",
		"Forget",
		"Remember",
		"Listen",
		"Elation",
		"Sacrifice",
		"Success",
		"Prosperity",
		"Apologize",
		"Strength",
		"Clarity",
		"Confusion",
		"Error",
		"Sorrow",
		"Hope",
		"Trust",
		"Change"
		)
	n <- sample(WiseWords, 2) #sample = randomly picks! 
	return(paste(n, collapse=" "))
}

Oracle()

# =============================================

MySqrt <- function(y) {
	x<-y/2
	z <- x
	while (abs(x*x - y) > 1e-10) {
		x <- (x + y/x)/2
		z <- c(z, x)
		}
	return(z)
}

initials <- function(str) {
	m <- gregexpr("(\\b\\w{1})",str, perl=T)
	q <- regmatches(str, m)
	paste(q[[1]], collapse='')
}

# [End]
