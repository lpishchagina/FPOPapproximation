#############################################################################################
#############################################################################################
##                                Plot Number of candidates                                ##
##             The number of candidates of change stored over time by FPOP                 ##
## (last, all), (random, all), (last, random), (random, random)     N = 100000             ##
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

#Type Of Penalty-------------------------
TypeOfPenalty = 1                                                               # or 0

# Length Of Data-------------------------                                                            # you can choose
LengthOfData <-1:100000

# Dimension------------------------------
Dim <- c(2)      # 2  3

#Noise-----------------------------------
Noise <- 1

#Penalty---------------------------------
Penalty <- ifelse( TypeOfPenalty == 0, 2*log(length(LengthOfData) ), 2*Dim*log(length(LengthOfData) ))

#Result Vector--------------------------- 



#files for result------------------------

PlotResNbCand <- NULL
PlotResNbCand <- paste('FPOP-methods Number of cadidates N = ', length(LengthOfData) , 'Dim = ', Dim, ' TP = ', TypeOfPenalty,'Plot.png') 

#-------------------------------------------------------------------------------
VectNbCands <- NULL
MatrResNbCands <- matrix(0, nrow = length(LengthOfData), ncol = length(Combination)+1)
MatrResNbCands [, 1] <- LengthOfData
#-------------------------------------------------------------------------------
VectNbCands <- readLines(con = 'Number of cadidates N =  1e+05 Dim =  2 A = rectangle I = last E = all  TP =  1 .txt', n = -1)
VectNbCands <- strsplit(VectNbCands,split = ' ')
VectNbCands <- sapply(VectNbCands, FUN = function(x) {as.double(unlist(x))})
MatrResNbCands [, 2] <- VectNbCands

VectNbCands <- readLines(con = 'Number of cadidates N =  1e+05 Dim =  2 A = rectangle I = random E = all  TP =  1 .txt', n = -1)
VectNbCands <- strsplit(VectNbCands,split = ' ')
VectNbCands <- sapply(VectNbCands, FUN = function(x) {as.double(unlist(x))})
MatrResNbCands [, 3] <- VectNbCands

VectNbCands <- readLines(con = 'Number of cadidates N =  1e+05 Dim =  2 A = rectangle I = last E = random  TP =  1 .txt', n = -1)
VectNbCands <- strsplit(VectNbCands,split = ' ')
VectNbCands <- sapply(VectNbCands, FUN = function(x) {as.double(unlist(x))})
MatrResNbCands [, 4] <- VectNbCands

VectNbCands <- readLines(con = 'Number of cadidates N =  1e+05 Dim =  2 A = rectangle I = random E = random  TP =  1 .txt', n = -1)
VectNbCands <- strsplit(VectNbCands,split = ' ')
VectNbCands <- sapply(VectNbCands, FUN = function(x) {as.double(unlist(x))})
MatrResNbCands [, 5] <- VectNbCands

frameMeanNbCands <- data.frame (MatrResNbCands)

PlotMean <- ggplot(frameMeanNbCands,  aes(MatrResNbCands[,1])) + geom_line(aes(y = MatrResNbCands[,2], col = Combination[1])) + geom_line(aes(y = MatrResNbCands[,3], col = Combination[2]))  + geom_line(aes(y = MatrResNbCands[,4], col = Combination[3]))  + geom_line(aes(y = MatrResNbCands[,5], col = Combination[4])) + labs( x = "Time", y = "Number of candidates being considered", title = paste('Number of change candidates (', NumberOfSimulations, 'iterations )')) + theme()
png(PlotResNbCand ,  width = 1500, height = 1000)
print(PlotMean)
dev.off()

#############################################################################################
#############################################################################################
##                                        End Plot                                         ##
#############################################################################################
#############################################################################################
