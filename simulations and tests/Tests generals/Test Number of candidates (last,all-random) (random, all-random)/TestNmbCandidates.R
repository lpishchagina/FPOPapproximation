#############################################################################################
#############################################################################################
##                                Test Number of candidates                                ##
##             The number of candidates of change stored over time by FPOP                 ##
#############################################################################################
#############################################################################################

# Packages and libraries--------------------------------------------------------

devtools::install_github("lpishchagina/FPOPapproximation")

library(FPOPapproximation)
library(base)
library(rstream)
library(tidyverse)
library("RColorBrewer")
library(iterators)
library(stats)
library(snow)
library(doSNOW)
library(foreach)

set.seed(13)
#-------------------------------------------------------------------------------

# Different values of FPOP-method parameters (using approximation)

# Rectangle-approximation---------

# approximation = 'rectangle', intersection = 'empty', exclusion = 'empty'      (PELT)
# approximation = 'rectangle', intersection = 'empty', exclusion = 'all'
# approximation = 'rectangle', intersection = 'all', exclusion = 'all'
# approximation = 'rectangle', intersection = 'all', exclusion = 'empty'
# approximation = 'rectangle', intersection = 'all', exclusion = 'random'
# approximation = 'rectangle', intersection = 'last', exclusion = 'all'         (using vector of links) +2,3
# approximation = 'rectangle', intersection = 'last', exclusion = 'allmodif'    (using list of disks)
# approximation = 'rectangle', intersection = 'last', exclusion = 'allmodif2'   (using full list of disks)
# approximation = 'rectangle', intersection = 'last', exclusion = 'ramdom'
# approximation = 'rectangle', intersection = 'random', exclusion = 'all'
# approximation = 'rectangle', intersection = 'random', exclusion = 'random'

# Simple-approximation---------

# approximation = 'sphere', intersection = 'last', exclusion = 'all'
# approximation = 'sphere', intersection = 'last', exclusion = 'random'
# approximation = 'sphere', intersection = 'random', exclusion = 'random'
# approximation = 'sphere', intersection = 'random', exclusion = 'all'
#-------------------------------------------------------------------------------

#  Runtime of one implementation of FPOP-method for data without changepoints---

#'@title ParameterCheck
#'
#' @description Check of the parameters for the function RuntimeFpopOneSimulation.
#' @param n is the length of data (positive)
#' @param dim is the dimension of data (positive)
#' @param typeOfPenalty is the type of the value of penalty: if 'typeOfPenalty = 0' then 'penality = 2log(n)'; if 'typeOfPenalty = 1' then 'penality = 2dim*log(n)' (by default, 1).
#' @param approximation is the type of approximation : 'rectangle' or 'sphere' (by default, 'rectangle') for the FPOP-method.
#' @param intersection is the type of intersection : 'empty', 'all', 'last' or 'random' (by default, 'last') for the FPOP-method.
#' @param exclusion is the type of intersection : 'empty', 'all'or 'random'(by default, 'all') for the FPOP-method.
#' @param func is the function of the implementation of FPOP-method (by  default, 'approxFPOP').
#' @param noise standard deviation of an additional normal noise (by default 'noise = 1').
#'
#' @return TRUE if the parameters are correct, else FALSE.
#'
#' @example
#' ParameterCheck(100, -2, 1, Approximation, Intersection, Exclusion,  func = "approxFPOP", 1)
#' ParameterCheck(100, 2, 1, Approximation, Intersection, Exclusion,  func = "approxFPOP", 1)
ParameterCheck <- function(n, p, typeOfPenalty, approximation, intersection, exclusion,  func = "approxFPOP", noise) {
  #---stop---#
  if ( !(is.double(n)) || !(is.double(p)) || !(is.double(typeOfPenalty))) {
    stop('The values of n, p, typeOfPenality  are integers!')
    return (FALSE)
  }
  if ((n < 1) || (p < 1)) {
    stop('The values of n, p  are positives!')
    return (FALSE)
  }
  if (noise < 0) {
    stop('The noise  is non-negative!')
    return (FALSE)
  }
  if ( (typeOfPenalty != 1) && (typeOfPenalty != 0)) {
    stop('The values of typeOfPenality is  0 ou 1!')
    return (FALSE)
  }
  if (!(approximation %in% c('rectangle', 'sphere'))) {
    stop('The value of approximation  is rectangle or sphere!')
    return (FALSE)
  }
  if (!(intersection %in% c('all', 'last', 'random', 'empty'))) {
    stop('The value of intersection is empty, all, last or random!')
    return (FALSE)
  }
  if (!(exclusion %in% c('all', 'last', 'random', 'empty', 'allmodif', 'allmodif2'))) {
    stop('The value of intersection is empty, all, last or random!')
    return (FALSE)
  }
  if (!(func %in% c('approxFPOP'))) {
    stop('The value of func is approxFPOP!')
    return (FALSE)
  }
  return (TRUE)
}
#-------------------------------------------------------------------------------
# Selecting parameters----------------------------------------------------------

# Number Of Simulations------------------
NumberOfSimulations <- 10                                                       # you can choose Number of simulation for calculation

# FPOP-method----------------------------
Approximation = 'rectangle'
Intersection = 'last'
Exclusion  = 'all'

Combination = c('(last,random)')
#Type Of Penalty-------------------------
TypeOfPenalty = 1                                                               # or 0

# Length Of Data-------------------------                                                            # you can choose
LengthOfData <-c(100000)

# Dimension------------------------------
Dim <- c(2)      # 2  3

#Noise-----------------------------------
Noise <- 1

#Penalty---------------------------------
Penalty <- ifelse( TypeOfPenalty == 0, 2*log(LengthOfData ), 2*Dim*log(LengthOfData ))

#Result Vector--------------------------- 
VectMeanNbCands <- c(1:LengthOfData)


#files for result------------------------
fileResNbCand <- NULL
fileResNbCand <- paste('Number of cadidates N = ', LengthOfData, 'Dim = ', Dim,  'A =', Approximation, 'I =', Intersection, 'E =', Exclusion,  ' TP = ', TypeOfPenalty,'.txt')

PlotResNbCand <- NULL
PlotResNbCand <- paste('Number of cadidates N = ', LengthOfData, 'Dim = ', Dim,  'A =', Approximation, 'I =', Intersection, 'E =', Exclusion,  ' TP = ', TypeOfPenalty,'Plot.png') 

#-------------------------------------------------------------------------------
MatrResNbCands <- matrix(0, nrow = LengthOfData, ncol = NumberOfSimulations)

#-------------------------------------------------------------------------------
for (i in 1 : NumberOfSimulations) {
  timeSeries <- rnormChanges(Dim, LengthOfData, changes = NULL, means = matrix(0, ncol = 1, nrow = Dim), Noise)
  approxFpop(timeSeries, Penalty, Approximation, Intersection, Exclusion, NbOfCands = TRUE)
  VectNbCands <- readLines(con = 'NbOfCands.txt', n = -1)
  VectNbCands <- strsplit(VectNbCands,split = ' ')
  VectNbCands <- sapply(VectNbCands, FUN = function(x) {as.integer(unlist(x))})
  MatrResNbCands[, i] <- VectNbCands
}
VectMeanNbCands <- apply(MatrResNbCands, 1, mean)
write.table(VectMeanNbCands, fileResNbCand,  row.names = FALSE, col.names = FALSE)
time <- c(1 : LengthOfData)
frameMeanNbCands <- data.frame (time, VectMeanNbCands)

PlotMean <- ggplot(frameMeanNbCands,  aes(time)) + geom_line(aes(y = frameMeanNbCands[, 2], col = Combination)) + labs( x = "Time", y = "Number of candidates being considered", title = paste('Number of change candidates (', NumberOfSimulations, 'iterations )'))+theme()
png(PlotResNbCand ,  width = 1500, height = 1000)
print(PlotMean)
dev.off()

#############################################################################################
#############################################################################################
##                                        End Test                                         ##
#############################################################################################
#############################################################################################
