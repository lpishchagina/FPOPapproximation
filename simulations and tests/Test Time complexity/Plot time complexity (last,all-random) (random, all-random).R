#############################################################################################
#############################################################################################
##                                Plot Time complexity                                     ##
## (last, all), (random, all), (last, random), (random, random)                            ##
#############################################################################################
#############################################################################################

# Packages and libraries--------------------------------------------------------

library(base)
library(rstream)
library(tidyverse)
library("RColorBrewer")

set.seed(13)
#-------------------------------------------------------------------------------

# Different values of FPOP-method parameters (using approximation)

# Rectangle-approximation---------

# approximation = 'rectangle', intersection = 'last', exclusion = 'all'         
# approximation = 'rectangle', intersection = 'last', exclusion = 'random'
# approximation = 'rectangle', intersection = 'random', exclusion = 'all'
# approximation = 'rectangle', intersection = 'random', exclusion = 'random'
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Selecting parameters----------------------------------------------------------
Combination = c('(last, all)',' (random, all)', '(last, random)', '(random, random)' )

# Type Of Penalty-------------------------
TypeOfPenalty = 1                                                               # or 0

# Length Of Data-------------------------
LimDegree <- 3:6                                                               
LengthOfData <-10^(LimDegree)


# Dimension------------------------------
Dim <- c(2)      # 2  3

#Noise-----------------------------------
Noise <- 1

#Penalty---------------------------------
Penalty <- ifelse( TypeOfPenalty == 0, 2*log(length(LengthOfData) ), 2*Dim*log(length(LengthOfData) ))

#Result Vector--------------------------- 



#files for result------------------------

PlotResTimeComplexity <- NULL
PlotResTimeComplexity <- paste('FPOP-methods TimeComplexity Dim = ', Dim, ' TP = ', TypeOfPenalty,'Plot.png') 

PlotReslgTimeComplexity <- NULL
PlotReslgTimeComplexity <- paste('FPOP-methods lg TimeComplexity Dim = ', Dim, ' TP = ', TypeOfPenalty,'Plot.png') 

#-------------------------------------------------------------------------------
VectTimeComplexity <- NULL
MatrResTimeComplexity <- matrix(0, nrow = length(LengthOfData), ncol = length(Combination)+1)
MatrResTimeComplexity [, 1] <- LengthOfData
#-------------------------------------------------------------------------------

VectTimeComplexity <- readLines(con = 'Runtime Mean Dim = 2 A = rectangle I = last E = all  TP =  1 .txt', n = -1)
VectTimeComplexity <- strsplit(VectTimeComplexity,split = ' ')
VectTimeComplexity <- sapply(VectTimeComplexity, FUN = function(x) {as.double(unlist(x))})
MatrResTimeComplexity [, 2] <- VectTimeComplexity

VectTimeComplexity <- readLines(con = 'Runtime Mean Dim = 2 A = rectangle I = random E = all  TP =  1 .txt', n = -1)
VectTimeComplexity <- strsplit(VectTimeComplexity,split = ' ')
VectTimeComplexity <- sapply(VectTimeComplexity, FUN = function(x) {as.double(unlist(x))})
MatrResTimeComplexity [, 3] <- VectTimeComplexity

VectTimeComplexity <- readLines(con = 'Runtime Mean Dim = 2 A = rectangle I = last E = random  TP =  1 .txt', n = -1)
VectTimeComplexity <- strsplit(VectTimeComplexity,split = ' ')
VectTimeComplexity <- sapply(VectTimeComplexity, FUN = function(x) {as.double(unlist(x))})
MatrResTimeComplexity [, 4] <- VectTimeComplexity

VectTimeComplexity <- readLines(con = 'Runtime Mean Dim = 2 A = rectangle I = random E = random  TP =  1 .txt', n = -1)
VectTimeComplexity <- strsplit(VectTimeComplexity,split = ' ')
VectTimeComplexity <- sapply(VectTimeComplexity, FUN = function(x) {as.double(unlist(x))})
MatrResTimeComplexity [, 5] <- VectTimeComplexity

frameMeanTimeComplexity <- data.frame (MatrResTimeComplexity)

PlotMean <- ggplot(frameMeanTimeComplexity,  aes(MatrResTimeComplexity[,1])) + geom_line(aes(y = MatrResTimeComplexity[,2], col = Combination[1])) + geom_line(aes(y = MatrResTimeComplexity[,3], col = Combination[2]))  + geom_line(aes(y = MatrResTimeComplexity[,4], col = Combination[3]))  + geom_line(aes(y = MatrResTimeComplexity[,5], col = Combination[4])) + labs( x = "Time", y = "Runtime", title = paste('Time complexity')) + theme()
png(PlotResTimeComplexity ,  width = 1500, height = 1000)
print(PlotMean)
dev.off()
#-------------------------------------------------------------------------------
PlotlgMean <- ggplot(frameMeanTimeComplexity,  aes(log10(MatrResTimeComplexity[,1]))) + geom_line(aes(y = log10(MatrResTimeComplexity[,2]), col = Combination[1])) + geom_line(aes(y = log10(MatrResTimeComplexity[,3]), col = Combination[2]))  + geom_line(aes(y = log10(MatrResTimeComplexity[,4]), col = Combination[3]))  + geom_line(aes(y = log10(MatrResTimeComplexity[,5]), col = Combination[4])) + labs( x = "lg(Time)", y = "lg(Runtime)", title = paste('Time complexity')) + theme()
png(PlotReslgTimeComplexity ,  width = 1500, height = 1000)
print(PlotlgMean)
dev.off()
#############################################################################################
#############################################################################################
##                                        End Plot                                         ##
#############################################################################################
#############################################################################################