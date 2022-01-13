#############################################################################################
#############################################################################################
##                                Test  Complexity                                         ##
##      (approximation = rectangle, intersection = last, exclusion = all-modif-modif2)     ##
#############################################################################################
#############################################################################################

# Packages and libraries--------------------------------------------------------

#devtools::install_github("lpishchagina/FPOPapproximation")

library(FPOPapproximation)
library(base)
library(rstream)
library(tidyverse)
library(ggpubr)
library("RColorBrewer")
library(iterators)
library(stats)
library(snow)
library(foreach)
library(doSNOW)

set.seed(13)
#-------------------------------------------------------------------------------

# Different values of FPOP-method parameters (using approximation)

# Rectangle-approximation---------
# Rectangle-approximation---------

# approximation = 'rectangle', intersection = 'all', exclusion = 'all'          (using left part of vector)
# approximation = 'rectangle', intersection = 'all', exclusion = 'allmodif'     (using all exclusions 1,.., tau-1)
# approximation = 'rectangle', intersection = 'all', exclusion = 'allmodif2'    (using list of disks)
#-------------------------------------------------------------------------------
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

#'@title RuntimeFpopOneSimulation
#'
#' @description Runtime of one implementation of FPOP-method for the multiple data of length N without changepoints (using approximation of the sets).
#' @param n is the length of data (positive)
#' @param dim is the dimension of data (positive)
#' @param typeOfPenalty is the type of the value of penalty: if 'typeOfPenalty = 0' then 'penality = 2log(n)'; if 'typeOfPenalty = 1' then 'penality = 2dim*log(n)' (by default, 1).
#' @param approximation is the type of approximation : 'rectangle' or 'sphere' (by default, 'rectangle') for the FPOP-method.
#' @param intersection is the type of intersection : 'empty', 'all', 'last' or 'random' (by default, 'last') for the FPOP-method.
#' @param exclusion is the type of intersection : 'empty', 'all'or 'random'(by default, 'all') for the FPOP-method.
#' @param func is the function of the implementation of FPOP-method (by  default, 'approxFPOP').
#' @param noise standard deviation of an additional normal noise (by default 'noise = 1').
#'
#' @return the runtime of one implementation of FPOP-method for the multiple data of length N without changepoints (using approximation of the sets).
#' @example
#' RuntimeFpopOneSimulation(10000, 2, typeOfPenalty = 1, 'rectangle', 'last', 'all')
#' RuntimeFpopOneSimulation(10000, 2, typeOfPenalty = 1, 'rectangle', 'last', 'random')
RuntimeFpopOneSimulation <- function(n, p, typeOfPenalty = 1, approximation, intersection, exclusion,  func = "approxFPOP", noise = 1) {
  #if (ParameterCheck(n, p, typeOfPenalty, approximation, intersection, exclusion,  func = "approxFPOP", noise)) {
  penalty <- ifelse( typeOfPenalty == 0, 2*log(n), 2*p*log(n))
  timeSeries <- rnormChanges(p, n, changes = NULL, means = matrix(0, ncol = 1, nrow = p), noise)
  res <- system.time(approxFpop(timeSeries, penalty, approximation, intersection,exclusion))[[1]]
  return (res);
  #} #end for if (ParameterCheck)

}
#-------------------------------------------------------------------------------

# Cluster reservation for parallel execution of a method------------------------
Cluster <- makeCluster(7, type = "SOCK", outfile="")
registerDoSNOW(Cluster)
#-------------------------------------------------------------------------------

# Selecting parameters----------------------------------------------------------

# Number Of Simulations------------------
NumberOfSimulations <- 35                                                        # you can choose Number of simulation for calculation of Runtime Mean

# FPOP-method----------------------------
Approximation = 'rectangle'
Intersection = 'last'
#Exclusion  = 'all'                                                             #MODIF1
#Exclusion  = 'allmodif'                                                        #MODIF2
#Exclusion  = 'allmodif2'                                                       #MODIF3
#Type Of Penalty-------------------------
TypeOfPenalty = 1                                                               # or 0

# Length Of Data-------------------------

LengthOfData <- c(1000, 3000, 5000,10000, 15000, 20000, 25000)

# Dimension------------------------------
Dim <- 2                                                                         # 2  3
#Noise-----------------------------------
Noise <- 1

#-------------------------------------------------------------------------------

# Table of Runtime for different lengths----------------------------------------
# (size: length(LimDegree)xNumberOfSimulations)
TableOfRuntime <- data.frame(matrix(0, length(LengthOfData),  NumberOfSimulations + 1))
colnames(TableOfRuntime) <- c("n", paste0("Rep", 1 : NumberOfSimulations ))
TableOfRuntime[, 1] <- LengthOfData

#table filling-------------------------
TableRes <- NULL
for(j in 1 : length(LengthOfData)) {
  N <- LengthOfData[j]
  TableOfRuntime[j,-1]  <- foreach(i = 1 : NumberOfSimulations, .combine = c, .packages=c("rstream", "tidyverse","FPOPapproximation")) %dopar% { RuntimeFpopOneSimulation (N, Dim, TypeOfPenalty, Approximation, Intersection, Exclusion,  func = "approxFPOP", noise = Noise) }
}

# Time complexity---------------------------------------------------------------
MeanOfRuntime <- rowMeans(TableOfRuntime[,-1])
#-------------------------------------------------------------------------------

#files for result------------------------
PlotOflgRes <- NULL
PlotOflgRes <- paste('RuntimeLg Dim =', Dim, 'A =', Approximation, ' I =', Intersection, ' E =', Exclusion,  ' TP = ', TypeOfPenalty, 'Plot .png' )
PlotOfRes <- NULL
PlotOfRes <-   paste('Runtime Dim =', Dim, 'A =', Approximation, ' I =', Intersection, ' E =', Exclusion,  ' TP = ', TypeOfPenalty, 'Plot .png' )
ResOfMeanRuntime <- NULL
ResOfMeanRuntime <- paste('Runtime Mean Dim =', Dim, 'A =', Approximation, 'I =', Intersection, 'E =', Exclusion,  ' TP = ', TypeOfPenalty,'.txt' )
ResTabRuntime <- NULL
ResTabRuntime <- paste('Runtime Dim =', Dim, 'A =', Approximation, 'I =', Intersection, 'E =', Exclusion,  ' TP = ', TypeOfPenalty,'.txt' )


# Save results
write.table(TableOfRuntime, ResTabRuntime  )
write.table(MeanOfRuntime, ResOfMeanRuntime, row.names = FALSE, col.names = FALSE  )

png(PlotOflgRes,  width = 1500, height = 1000)
plot(log10(LengthOfData), log10(MeanOfRuntime), xlab = "lg(data length)", ylab = "lg(mean Runtime) in second",  main = "Time complexity", col = "royalblue3")
lines(log10(LengthOfData), log10(MeanOfRuntime), col = "royalblue3", lwd = 3)
dev.off()

png(PlotOfRes,  width = 1500, height = 1000)
plot(LengthOfData, MeanOfRuntime, xlab = "data length", ylab = "Runtime mean in second",  main = "Time complexity", col = "royalblue3")
lines(LengthOfData, MeanOfRuntime, col = "royalblue3", lwd = 3)
dev.off()

# Stop cluster------------------------------------------------------------------
stopCluster(Cluster)
#---


#ATTENTION!!!!STOP!!!



#TabResModifAllExcl <- data.frame(matrix(0, length(LengthOfData),4))
#colnames(TabResModifAllExcl) <- c("n", "all(vector)",all, "all(list)")
#TabResModifAllExcl[, 1] <- LengthOfData
#TabResModifAllExcl[, 2] <- MeanOfRuntime
#TabResModifAllExcl[, 3] <- MeanOfRuntime
#TabResModifAllExcl[, 4] <- MeanOfRuntime



fileTabResModifAllExcl <- NULL
fileTabResModifAllExcl <- paste('Table_Runtime_mean_3_modif_Dim_=', Dim,  '(rectangle,last,all(vector)-all-all(list))', ' TP_=_', TypeOfPenalty,'.txt' )

plotTabResModifAllExcl <- NULL
plotTabResModifAllExcl <- paste('Plot_Runtime_mean_3_modif_Dim_=', Dim,  '(rectangle,last,all(vector)-all-all(list))', ' TP_=_',TypeOfPenalty,'.png',sep="_" )

# Save results
write.table(TabResModifAllExcl, fileTabResModifAllExcl )


Combination <-c("('rectangle', I ='last', E ='all(vector)')", "('rectangle', I ='last', E ='all')", "('rectangle', I ='last', E ='all(list)')")

PlotMeanRuntime <- ggplot(TabResModifAllExcl,  aes(TabResModifAllExcl[, 1])) + geom_line(aes(y = TabResModifAllExcl[, 2], col = Combination [1])) + geom_line(aes(y = TabResModifAllExcl[,3],  col = Combination [2])) + geom_line(aes(y = TabResModifAllExcl[, 4],  col = Combination [3]) ) + labs( x = "Time", y = "Runtime (in seconds)", title = "Time complexity")+theme()

png(filename = plotTabResModifAllExcl,  width = 1500, height = 1000)
print(PlotMeanRuntime)
dev.off()


plotTabResModif2AllExcl <- NULL
plotTabResModif2AllExcl <- paste('Plot_Runtime_mean_2_modif_Dim_=', Dim,  '(rectangle,last,all(vector)-all(list))', ' TP_=_',TypeOfPenalty,'.png',sep="_" )


Plot2MeanRuntime <- ggplot(TabResModifAllExcl,  aes(TabResModifAllExcl[, 1])) + geom_line(aes(y = TabResModifAllExcl[, 2], col = Combination [1]))  + geom_line(aes(y = TabResModifAllExcl[, 4],  col = Combination [3]) ) + labs( x = "Time", y = "Runtime (in seconds)", title = "Time complexity")+theme()

png(filename = plotTabResModif2AllExcl,  width = 1500, height = 1000)
print(Plot2MeanRuntime)
dev.off()
##############
#lg
plotlgTabResModifAllExcl <- NULL
plotlgTabResModifAllExcl <- paste('Plot_lgRuntime_mean_3_modif_Dim_=', Dim,  '(rectangle,last,all(vector)-all-all(list))', 'TP_=_',TypeOfPenalty,'.png',sep="_" )

PlotlgMeanRuntime <- ggplot(TabResModifAllExcl,  aes(log10(TabResModifAllExcl[, 1]))) + geom_line(aes(y = log10(TabResModifAllExcl[, 2]), col = Combination [1])) + geom_line(aes(y = log10(TabResModifAllExcl[,3]),  col = Combination [2])) + geom_line(aes(y = log10(TabResModifAllExcl[, 4]),  col = Combination [3]) ) + labs( x = "lg(Time)", y = "lg(Runtime) (in seconds)", title = "Time complexity")+theme()

png(filename = plotlgTabResModifAllExcl,  width = 1500, height = 1000)
print(PlotlgMeanRuntime)
dev.off()


plotTabResModiflg2AllExcl <- NULL
plotTabResModiflg2AllExcl <- paste('Plot_lgRuntime_mean_2_modif_Dim_=', Dim,  '(rectangle,last,all(vector)-all(list))', 'TP_=_',TypeOfPenalty,'.png',sep="_" )


Plotlg2MeanRuntime <- ggplot(TabResModifAllExcl,  aes(log10(TabResModifAllExcl[, 1]))) + geom_line(aes(y = log10(TabResModifAllExcl[, 2]), col = Combination [1]))  + geom_line(aes(y = log10(TabResModifAllExcl[, 4]),  col = Combination [3]) ) + labs( x = "lg(Time)", y = "lg(Runtime) (in seconds)", title = "Time complexity")+theme()



png(filename = plotTabResModiflg2AllExcl,  width = 1500, height = 1000)
print(Plotlg2MeanRuntime)
dev.off()



