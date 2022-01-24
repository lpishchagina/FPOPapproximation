
N <- 10
Chpt <-5
Means <-  matrix(c(0,1,1,10), nrow = 2)
Noise <- 1
Dim <- 2
Penality <- 2*Dim*log(N)
time_series <- rnormChanges(p = Dim, n = N, changes = Chpt, means = Means, noise = Noise)

approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'all', exclusion = 'all', NbOfCands = TRUE, NbOfExclus = FALSE)

approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'last', exclusion = 'all', NbOfCands = TRUE, NbOfExclus = FALSE)


approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'last', exclusion = 'all', NbOfCands = TRUE, NbOfExclus = FALSE)

approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'last', exclusion = 'all', NbOfCands = TRUE, NbOfExclus = FALSE)

approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'last', exclusion = 'all', NbOfCands = TRUE, NbOfExclus = FALSE)


approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'last', exclusion = 'allmodif2', NbOfCands = TRUE, NbOfExclus = FALSE)



N <- 10
Chpt <-5
Means <-  matrix(c(0,1,1,10), nrow = 2)
Noise <- 1
Dim <- 2
Penality <- 2*Dim*log(N)
time_series <- rnormChanges(p = Dim, n = N, changes = Chpt, means = Means, noise = Noise)

approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'all', exclusion = 'all', NbOfCands = TRUE, NbOfExclus = FALSE)

approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'last', exclusion = 'all', NbOfCands = TRUE, NbOfExclus = FALSE)


approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'last', exclusion = 'all', NbOfCands = TRUE, NbOfExclus = FALSE)

approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'last', exclusion = 'all', NbOfCands = TRUE, NbOfExclus = FALSE)

approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'last', exclusion = 'all', NbOfCands = TRUE, NbOfExclus = FALSE)


approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'last', exclusion = 'allmodif2', NbOfCands = TRUE, NbOfExclus = FALSE)



N <- 1000
Chpt <-5
Means <-  matrix(c(0,1,1,10), nrow = 2)
Noise <- 1
Dim <- 2
time_series <- rnormChanges(p = 2, n = N, changes = NULL, means = matrix(0, ncol = 1, nrow = 2), noise = 1)

res<- list()
res[[1]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'last', exclusion = 'all', NbOfCands = TRUE, NbOfExclus = FALSE)

res[[2]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'last', exclusion = 'allmodif2', NbOfCands = TRUE, NbOfExclus = FALSE)

time_res <-list()
time_res[[1]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'last', exclusion = 'all', NbOfCands = TRUE, NbOfExclus = FALSE))

time_res[[2]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'last', exclusion = 'allmodif2', NbOfCands = TRUE, NbOfExclus = FALSE))

nb1<- res[[1]]$NumberOfCandidats
nb2<- res[[2]]$NumberOfCandidats
nb1-nb2
