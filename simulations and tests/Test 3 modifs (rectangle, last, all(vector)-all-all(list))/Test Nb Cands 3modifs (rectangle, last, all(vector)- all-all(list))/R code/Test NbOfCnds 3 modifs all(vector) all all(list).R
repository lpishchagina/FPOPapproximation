#############################################################################################
#############################################################################################
##                                Test Number of candidates                                ##
##      (approximation = rectangle, intersection = last, exclusion = all-modif-modif2)     ##
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

# approximation = 'rectangle', intersection = 'all', exclusion = 'all'          (using left part of vector)
# approximation = 'rectangle', intersection = 'all', exclusion = 'allmodif'     (using all exclusions 1,.., tau-1)
# approximation = 'rectangle', intersection = 'all', exclusion = 'allmodif2'    (using list of disks)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Selecting parameters----------------------------------------------------------

# Number Of Simulations------------------
NumberOfSimulations <- 100                                                      # you can choose Number of simulation for calculation

# FPOP-method----------------------------
Approximation = 'rectangle'
Intersection = 'last'
#Exclusion  = 'all'                                                             #MODIF1
#Exclusion  = 'allmodif'                                                        #MODIF2
#Exclusion  = 'allmodif2'                                                       #MODIF3
#Type Of Penalty-------------------------
#Combination = c('(rectangle,last,all (vector))')
#Combination = c('(rectangle,last,all)')
#Combination = c('(rectangle,last,all(list))')

#Type Of Penalty-------------------------
TypeOfPenalty = 1                                                               # or 0

# Length Of Data-------------------------                                                            # you can choose
LengthOfData <-c(10000)

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
fileResNbCand <- paste('Number_Of_Candidates_N_=', LengthOfData, 'Dim_=', Dim,  'A_=', Approximation, 'I_=', Intersection, 'E_=', Exclusion,  'TP_=', TypeOfPenalty,'.txt', sep="_")

PlotResNbCand <- NULL
PlotResNbCand <- paste('Number_Of_Candidates_N =', LengthOfData, 'Dim_=', Dim,  'A_=', Approximation, 'I_=', Intersection, 'E_=', Exclusion,  'TP_=', TypeOfPenalty,'Plot.png', sep="_")

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


#ATTENTION!!!!STOP!!!
#TabResNbCandsModifAllExcl <- data.frame(matrix(0, LengthOfData,4))
#colnames(TabResNbCandsModifAllExcl ) <- c("time", "all(vector)", "all", "all(list)")
#TabResNbCandsModifAllExcl [, 1] <- 1:LengthOfData

#TabResNbCandsModifAllExcl[, 2] <-  VectMeanNbCands
#TabResNbCandsModifAllExcl[, 3] <-  VectMeanNbCands
#TabResNbCandsModifAllExcl[, 4] <-  VectMeanNbCands

# Result files-----------------------
fileTabResNbCandsModifAllExcl <- NULL
fileTabResNbCandsModifAllExcl <- paste('Table_Mean_Number_Of_Candidates_3_modif_Dim_=', Dim, 'I_=', Intersection, 'E_=_all', 'TP_=', TypeOfPenalty,'.txt', sep = "_" )

plotTabResNbCandsModifAllExcl<- NULL
plotTabResNbCandsModifAllExcl<- paste('Plot_Mean_Number_Of_Candidates_3_modif_Dim_=', Dim, 'I_=', Intersection, 'E_=_all', 'TP_=', TypeOfPenalty,'.png' , sep = "_" )

# Save results----------------------
write.table(TabResNbCandsModifAllExcl,fileTabResNbCandsModifAllExcl )


Combination <-c("('rectangle', I ='last', E ='all(vector)')", "('rectangle', I ='last', E ='all')", "('rectangle', I ='last', E ='all(list)')")

PlotMeanNbOfCands <- ggplot(TabResNbCandsModifAllExcl,  aes(TabResNbCandsModifAllExcl[, 1])) + geom_line(aes(y = TabResNbCandsModifAllExcl[, 2], col = Combination [1])) + geom_line(aes(y = TabResNbCandsModifAllExcl[,3],  col = Combination [2])) + geom_line(aes(y = TabResNbCandsModifAllExcl[, 4],  col = Combination [3]) ) + labs(  x = "Time", y = "Number of candidates", title = "Number of candidates")+theme()

png(plotTabResNbCandsModifAllExcl,  width = 1500, height = 1000)
print(PlotMeanNbOfCands)
dev.off()


plotTabResNbCandsAllAllmodif2 <- paste('Plot_Mean_Number_Of_Candidates_2_modif(vector etlist)_Dim_=', Dim, 'I_=', Intersection, 'E_=_all', 'TP_=', TypeOfPenalty,'.png' , sep = "_" )

PlotMeanNbOfCandsAllAllmodif2 <- ggplot(TabResNbCandsModifAllExcl,  aes(TabResNbCandsModifAllExcl[, 1])) + geom_line(aes(y = TabResNbCandsModifAllExcl[, 2], col = Combination [1])) + geom_line(aes(y = TabResNbCandsModifAllExcl[, 4],  col = Combination [3]) ) + labs( x = "Time", y = "Number of candidates", title = "Number of candidates")+theme()

png(plotTabResNbCandsAllAllmodif2 ,  width = 1500, height = 1000)
print(PlotMeanNbOfCandsAllAllmodif2)
dev.off()




#############################################################################################
#############################################################################################
##                                        End Test                                         ##
#############################################################################################
#############################################################################################
