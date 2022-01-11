#############################################################################################
#############################################################################################
##                                Test  Candidates                                         ##
##             The number of candidates of change stored over time by FPOP                 ##
##                            The data WITHOUT changes                                     ##
##                              See Simulation and tests\                                  ##
#############################################################################################
#############################################################################################

#install.packages("RColorBrewer")

#devtools::install_github("lpishchagina/FPOPapproximation")

library(FPOPapproximation)
library(base)
library(rstream)
library(tidyverse)
#library(ggpubr)
library("RColorBrewer")

NbRep = 20

Combination <-c("(I ='sphere', E ='sphere') ","(I ='all', E ='all')"
                , "(I ='all', E ='empty')", "(I ='empty', E ='all')"
                , "(I ='last', E ='all')", "(I ='last', E ='random')"
                , "(I ='all', E ='random')", "(I ='random', E ='random')"
                , "(I ='empty', E ='empty')")
N <- c(1000)  # c(100, 1000, 10000)
Dim <- c(2, 3, 4, 5, 6, 7, 8, 9, 10)
Noise <- 1

Intersection = c('sphere','all', 'all', 'empty', 'last', 'last', 'all', 'random', "empty" )
Exclusion = c('sphere', 'all', 'empty', 'all', 'all', 'random', 'random', 'random', 'empty')

set.seed(100)
for (pN in N){
  for (pDim in Dim) {
    Penality <- 2*pDim*log(pN)
    MatrMeanNbCands = matrix(0, nrow = pN, ncol = length(Intersection))
    ResMeanFile = paste('MeanNbOfCands all combinations N =', N, 'Dim =', pDim, '.txt' )
    for (i in 1 : length((Intersection))) {
      ResFile = paste('NbOfCands N =', pN, 'I =', Intersection[i], 'E =', Exclusion[i], 'Dim =', pDim, '.txt' )
      MatrNbCands = matrix(0, nrow = pN, ncol = NbRep)
      for (pNbRep in 1 : NbRep) {
        TestTS <- changes_rnorm(p = pDim, n = pN, changes = NULL, means = matrix(0, ncol = 1, nrow = pDim), noise = Noise)
        approx_fpop(data = TestTS, penalty = Penality, intersection = Intersection[i], exclusion = Exclusion[i], NbOfCands = TRUE, NbOfExclus = FALSE)
        VectNbCands <- readLines(con = 'NbOfCands.txt', n = -1)
        VectNbCands <- strsplit(VectNbCands,split = ' ')
        VectNbCands <- sapply(VectNbCands, FUN = function(x) {as.integer(unlist(x))})
        MatrNbCands[, pNbRep] <- VectNbCands
      }
      MatrMeanNbCands[,i] <- apply(MatrNbCands, 1, mean)
      write.table(MatrMeanNbCands[,i], ResFile,  row.names = FALSE, col.names = FALSE)
    }
    write.table(MatrMeanNbCands, ResMeanFile,  row.names = FALSE, col.names = FALSE)
    ResPlotFile <- paste('Plot MeansNbOfCands all combinations N =', N, 'Dim =', pDim, '.png' )
    titlePlot <- paste('PloMeansNbOfCands all combinations N =', N, 'Dim =', pDim, '.png' )
    time <- c(1 : pN)
    frameMeanNbCands <- data.frame (time, MatrMeanNbCands)
    #print(frameMeanNbCands)
    PlotMean <- ggplot(frameMeanNbCands,  aes(time)) + geom_line(aes(y = frameMeanNbCands[, 2], col = Combination [1])) + geom_line(aes(y = frameMeanNbCands[,3],  col = Combination [2])) + geom_line(aes(y = frameMeanNbCands[, 4],  col = Combination [3]) ) + geom_line(aes(y = frameMeanNbCands[, 5], col = Combination [4])) + geom_line(aes(y = frameMeanNbCands[, 6], col = Combination  [5])) + geom_line(aes(y = frameMeanNbCands[, 7], col = Combination [6])) + geom_line(aes(y = frameMeanNbCands[, 8], col = Combination [7]) ) + geom_line(aes(y = frameMeanNbCands[, 9], col = Combination [8])) + geom_line(aes(y = frameMeanNbCands[, 10], col = Combination [9])) + labs( x = "Time", y = "Number of candidates being considered", title = paste('Number of change candidates (', NbRep, 'iterations )'))+theme()
    png(filename = ResPlotFile,  width = 1500, height = 1000)
    print(PlotMean)
    dev.off()
    #print(PlotMean)
  }
}

#############################################################################################
#############################################################################################
##                                        End Test                                         ##
#############################################################################################
#############################################################################################
