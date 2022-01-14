#############################################################################################
#############################################################################################
##                                Test  Comparison of two Fpop-methods                     ##
#############################################################################################
#############################################################################################

# Packages and libraries--------------------------------------------------------

devtools::install_github("lpishchagina/FPOPapproximation")

library(FPOPapproximation)
library(base)
library(rstream)
library(tidyverse)
library(ggpubr)
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
#-------------------------------------------------------------------------------

# Cluster reservation for parallel execution of a method------------------------

Cluster <- makeCluster(7, type = "SOCK",outfile="")                             #  you can choose , for PC 7
registerDoSNOW(Cluster)
#-------------------------------------------------------------------------------

# Selecting parameters----------------------------------------------------------

# 1 FPOP-method----------------------------
Approximation1 = 'rectangle'
Intersection1 = 'random'
Exclusion1  = 'random'

# 2 FPOP-method----------------------------
Approximation2 = 'rectangle'
Intersection2 = 'last'
Exclusion2  = 'random'

# Type Of Penalty--------------------------
TypeOfPenalty = 1                                                               # or 0

# Length Of Data---------------------------                                                              # you can choose
LengthOfData <-100

# Dimension--------------------------------
Dim <- 2                                                                         # 2  3
#Noise-------------------------------------
Noise <- 1

# Comparison result------------------------
ResComparison <- TRUE;

# Number Of Simulations--------------------
NumberOfSimulations <- 100                         

# Changes and means------------------------
#Changes <- NULL
#Means <- matrix(0, ncol = 1, nrow = Dim)

Changes <- as.integer(LengthOfData / 2)
Means <-  matrix(c(0,1,1,10), nrow = 2)


#Penalty-----------------------------------
Penalty <- ifelse( TypeOfPenalty == 0, 2*log(LengthOfData ), 2*Dim*log(LengthOfData ))

#-------------------------------------------------------------------------------
res <- matrix(FALSE, ncol = 1, nrow = NumberOfSimulations)
res <- foreach(i = 1 : NumberOfSimulations, .combine = c, .packages=c("rstream", "tidyverse","FPOPapproximation")) %dopar% { TestTwoApproxFpop(data = rnormChanges(Dim, LengthOfData, Changes, Means, Noise), Penalty, Approximation1, Intersection1, Exclusion1, Approximation2, Intersection2, Exclusion2) }
ResComparison <- !(FALSE %in% res)

# Result------------------------------------------------------------------------  
ResComparison

# Stop cluster------------------------------------------------------------------
stopCluster(Cluster)
#-------------------------------------------------------------------------------
