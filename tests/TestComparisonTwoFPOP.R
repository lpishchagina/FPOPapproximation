#############################################################################################
#install.packages("RColorBrewer")
#devtools::install_github("lpishchagina/FPOPapproximation")
#library------------------------------------------------------------------------
library(FPOPapproximation)
library(base)
library(rstream)
library(tidyverse)
#library(ggpubr)
library("RColorBrewer")
#Combinations-------------------------------------------------------------------
approximation = 'rectangle'

NbRep = 20

Combination <-c("(I ='sphere', E ='sphere') ","(I ='all', E ='all')"
                , "(I ='all', E ='empty')", "(I ='empty', E ='all')"
                , "(I ='last', E ='all')", "(I ='last', E ='random')"
                , "(I ='all', E ='random')", "(I ='random', E ='random')"
                , "(I ='empty', E ='empty')")
N <- c(1000)  # c(100, 1000, 10000)
Dim <- c(2, 3, 4, 5, 6, 7, 8, 9, 10)
Noise <- 1
