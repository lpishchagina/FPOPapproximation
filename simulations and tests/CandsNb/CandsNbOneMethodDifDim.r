#############################################################################################
#############################################################################################
##                                Test Number of candidates                                ##
##             The number of candidates of change stored over time by FPOP                 ##
##                               one methode, dimensions 2,3,4,5,7,10, 12, 15, 20          ##
##                            (using parallel calculation)                                 ##
#############################################################################################
#############################################################################################

# Packages and libraries--------------------------------------------------------
# install.packages('')
#devtools::install_github("lpishchagina/FPOPapproximation")

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

set.seed(21)
#-------------------------------------------------------------------------------
# Number of candidates in each iteration in the FPOP-method----------------------
#'@title CandsNbOneSimu
#'
#' @description Number of candidates in each iteration in the FPOP-method(using approximation).
#' @param n is the length of data (positive).
#' @param dim is the dimension of data (positive).
#' @param typeOfPenalty is the type of the value of penalty: if 'typeOfPenalty = 0' then 'penalty = 2log(n)'; if 'typeOfPenalty = 1' then 'penalty = 2dim*log(n)' (by default, 1).
#' @param approximation is the type of approximation : 'rectangle' or 'sphere' (by default, 'rectangle') for the FPOP-method.
#' @param intersection is the type of intersection : 'empty', 'all', 'last' or 'random' (by default, 'last') for the FPOP-method.
#' @param exclusion is the type of intersection : 'empty', 'all'or 'random'(by default, 'all') for the FPOP-method.
#' @param noise standard deviation of an additional normal noise (by default 'noise = 1').
#' @return Number of candidates in each iteration (vector)  in the FPOP-method(using approximation).
#'
#' @example
#' CandsNbOneSimu(n = 100, dim = 2, typeOfPenalty = 1, approximation = 'rectangle', intersection = 'last', exclusion = 'all', noise = 1)
#' 
CandsNbOneSimu <- function(n, dim, typeOfPenalty = 1, approximation, intersection, exclusion, noise = 1) {
  penalty <- ifelse( typeOfPenalty == 0, 2*log(n), 2*dim*log(n))
  timeSeries <- rnormChanges(dim, n, changes = NULL, means = matrix(0, ncol = 1, nrow = dim), noise)
  res <-approxFpop(timeSeries, penalty, approximation, intersection, exclusion, NbOfCands = TRUE)$NumberOfCandidats
  return (res);
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Cluster reservation for parallel calculation----------------------------------
Cluster <- makeCluster(40, type = "SOCK", outfile="") # ATTENTION!!! 40 number of cores 
registerDoSNOW(Cluster)

# Сhoice of approximation for the FPOP-method-----------------------------------
# You can choose a combination of parameters from those offered (ATTENTION!!! REMOVE #) 

# Simple approximation---------------
#1.
Approximation = 'sphere'
Intersection = 'last'
Exclusion = 'all'

#2.
#Approximation = 'sphere'
#Intersection = 'random'
#Exclusion = 'all'

#3.
#Approximation = 'sphere'
#Intersection = 'random'
#Exclusion = 'random'

#4.
#Approximation = 'sphere'
#Intersection = 'last'
#Exclusion = 'random'

# Rectangle approximation---------------
# "one" Intersection----
#5.
#Approximation = 'rectangle'
#Intersection = 'last'
#Exclusion = 'random'

#6.
#Approximation = 'rectangle'
#Intersection = 'random'
#Exclusion = 'random'

#7.
# Approximation = 'rectangle'
# Intersection = 'random'
# Exclusion = 'all'  

#8.(using link vector)
# Approximation = 'rectangle'
#Intersection = 'last'
#Exclusion = 'all'

#9.(using link disk list and check 'inclusion')
# Approximation = 'rectangle'
#Intersection = 'last'
#Exclusion = 'allmodif2' 

#10.(using disk list)
# Approximation = 'rectangle'
#Intersection = 'last'
#Exclusion = 'allmodif' 

# "all" Intersections--
#11.
#Approximation = 'rectangle'
#Intersection = 'all'
#Exclusion = 'all'

#12.
#Approximation = 'rectangle'
#Intersection = 'all'
#Exclusion = 'random'

#13.
#Approximation = 'rectangle'
#Intersection = 'all'
#Exclusion = 'empty'

# PELT-----------------------------
#14.
# Approximation = 'rectangle'
# Intersection = 'empty'
# Exclusion = 'empty'    
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Сhoice of length Of data, dimension, noise, type of penalty and  number of simulations---------
Length <-c(10)
Dim <- c(2,3,4,5,7,10, 12, 15, 20)

Dimension <- NULL
for (k in 1 : length(Dim)) {
  Dimension <- c(Dimension, paste('dim =',Dim[k]))
} 

Noise <- 1
tPen = 1 # 0 # if '0' then 'penalty = 2log(n)'; if '1' then 'penalty = 2dim*log(n)'
NbSimus <- 200#0
#CandsNbDifDim------------------------------------------------------------------
Combination <- paste('(', Approximation ,',',Intersection ,',',Exclusion,')')
CandsNbDifDim <- matrix(0, Length, length(Dim))
CandsNbSimus <- matrix(0, Length, NbSimus)
for (k  in 1 : length(Dim)) {
  # result for one dimention-----------
  CandsNbSimus [,] <- foreach(i = 1 : NbSimus, .combine = c, .packages=c("rstream", "tidyverse","FPOPapproximation")) %dopar% { CandsNbOneSimu (Length, Dim[k], tPen, Approximation, Intersection, Exclusion, Noise) }
  CandsNbDifDim[,k] <- rowMeans(CandsNbSimus)
  #------------------------------------
  ##Recording results for one dimension
  ##you can save the results for each dimension otherwise comment out this section
  #CandsNb <- data.frame (time = c(1 : Length), CandsNbDifDim[,k])
  #PlotCandsNb <- ggplot(CandsNb, aes(time)) + geom_line(aes(y = CandsNb[, 2], col = Combination)) + labs( x = "Time", y = "Number of candidates being considered", title = paste('Number of change candidates (', NbSimus, 'iterations, dimension =', Dim[k],')',))+theme()
  #rFileRes <- paste('CandsNb_it', NbSimus, Dim[k],Length, tPen, Approximation, Intersection, Exclusion,'.RData', sep = '_')
  #txtFileRes <- paste('CandsNb_it', NbSimus, Dim[k],Length, tPen, Approximation, Intersection, Exclusion,'.txt', sep = '_')
  #pngFileRes <- paste('PlotCandsNb_it', NbSimus, Dim[k],Length, tPen, Approximation, Intersection, Exclusion,'.png', sep = '_')
  #write.table(CandsNb, txtFileRes, row.names = FALSE, col.names = FALSE)
  #save(CandsNb, file=rFileRes)
  #png(pngFileRes , width = 1500, height = 1000)
  #print(PlotCandsNb)
  #dev.off()
  #------------------------------------
}
CandsNbDifDim <- data.frame (time = c(1 : Length), CandsNbDifDim)
colnames(CandsNbDifDim) <- c('n', Dimension)
#Recording results 
rFileRes <- paste('CandsNbDifDim_it', NbSimus,Length, tPen, Approximation, Intersection, Exclusion,'.RData', sep = '_')
txtFileRes <- paste('CandsNbDifDim_it', NbSimus,Length, tPen, Approximation, Intersection, Exclusion,'.txt', sep = '_')
pngFileRes <- paste('PlotCandsNbDifDim_it', Dim[k],Length, tPen, Approximation, Intersection, Exclusion,'.png', sep = '_')

write.table(CandsNbDifDim, txtFileRes, row.names = FALSE, col.names = TRUE)
save(CandsNbDifDim, file=rFileRes)

PlotCandsNbDifDim <- ggplot(CandsNbDifDim, aes(n))+labs( x = "Time", y = "Number of candidates being considered", title = paste('Number of change candidates for', Combination, '(', NbSimus, 'iterations )'))+theme()+geom_line(aes(y = CandsNbDifDim[, 2], col =Dimension[1]))+geom_line(aes(y = CandsNbDifDim[, 3], col =Dimension[2]))+geom_line(aes(y = CandsNbDifDim[, 4], col =Dimension[3]))+geom_line(aes(y = CandsNbDifDim[, 5], col =Dimension[4]))+geom_line(aes(y = CandsNbDifDim[, 6], col =Dimension[5]))+geom_line(aes(y = CandsNbDifDim[, 7], col =Dimension[6]))+geom_line(aes(y = CandsNbDifDim[, 8], col =Dimension[7]))+geom_line(aes(y = CandsNbDifDim[, 9], col =Dimension[8]))+geom_line(aes(y = CandsNbDifDim[, 10], col =Dimension[9]))

print(PlotCandsNbDifDim)
png(pngFileRes , width = 1500, height = 1000)
print(PlotCandsNbDifDim)
dev.off()
#-------------------------------------------------------------------------------
stopCluster(Cluster)
################################################################################
########################### END ################################################
################################################################################
