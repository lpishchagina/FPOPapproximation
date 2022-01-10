N <- 100
Chpt <-NULL
Means <-  matrix(c(0,1,1,10), nrow = 2)
Noise <- 1
Dim <- 2
Penality <- 2*Dim*log(N)
time_series <- rnormChanges(p = Dim, n = N, changes = Chpt, means =  matrix(0, ncol = 1, nrow = Dim), noise = Noise)

Approx <- list()
Approx[[1]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'all', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE)
Approx[[2]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'all', exclusion = 'empty', NbOfCands = FALSE, NbOfExclus = FALSE)
Approx[[3]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'all', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE)

Approx[[4]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'empty', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE)

Approx[[5]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'last', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE)
Approx[[6]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'last', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE)

Approx[[7]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle',intersection = 'random', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE)
Approx[[8]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'random', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE)

Approx[[9]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'last', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE)
Approx[[10]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'random', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE)
Approx[[11]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'random', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE)
Approx[[12]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'last', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE)
Approx[[13]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'empty', exclusion = 'empty', NbOfCands = FALSE, NbOfExclus = FALSE)
Approx


ResultApprox <- function (ListApprox = list()) {
  for (i in 1:(length(ListApprox)-1)) {
    if (ListApprox[[i]]$changes != ListApprox[[i + 1]]$changes) return (FALSE)
    if (ListApprox[[i]]$UnpenalizedCost != ListApprox[[i + 1]]$UnpenalizedCost) return (FALSE)
  }
  return (TRUE)
}
res <-ResultApprox(Approx)

TimeApprox <- list()
TimeApprox[[1]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'all', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE))
TimeApprox[[2]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'all', exclusion = 'empty', NbOfCands = FALSE, NbOfExclus = FALSE))
TimeApprox[[3]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'all', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE))

TimeApprox[[4]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'empty', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE))

TimeApprox[[5]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'last', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE))
TimeApprox[[6]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'last', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE))

TimeApprox[[7]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle',intersection = 'random', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE))
TimeApprox[[8]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'random', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE))

TimeApprox[[9]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'last', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE))
TimeApprox[[10]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'random', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE))
TimeApprox[[11]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'random', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE))
TimeApprox[[12]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'last', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE))

TimeApprox[[13]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'empty', exclusion = 'empty', NbOfCands = FALSE, NbOfExclus = FALSE))
TimeApprox
