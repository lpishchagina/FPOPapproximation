<a id="top"></a>
#  FPOPapproximation Vignette
### Liudmila Pishchagina
### Novembre 15, 2021

## Quick Start

` FPOPapproximation ` is an R package written in Rcpp/C++ and developed to detection changes using the Functional Pruning Optimal Partitioning method (FPOP) in `p`-variate time series of length ` n `. 

The FPOP method consists in analysing the subsets ` Zit, i = 1,.., t ` at each iteration ` t `. The subset ` Zit ` is formed by the exclusion of the union of ` Sjt, j = 1,.., i-1 ` from the intersection of ` Sij, j = i,.., t `. In Gaussian case the sets ` S`  are `p`-balls. For dimension ` p > 1 ` the detection of a nonempty intersection and union between `p`-balls is difficult to solve. One way to solve this problem is to approximate the subsets ` Zit `.

To optimize this method we use the rectangle approximation and the different modifications of the number of intersections and exclusions.

Currently, the following parameter combinations (intersection, exclusion) are implemented:

` intersection = 'sphere', exclusion = 'sphere' ` (without rectangle approximation);

` intersection = 'all', exclusion = 'all' ` (all possible intersections and exclusions)

` intersection = 'all', exclusion = 'empty' ` (all possible intersections without exclusions)

` intersection = 'empty', exclusion = 'all' ` (rectangle approximation of the last ball and all possible exclusions)

` intersection = 'last', exclusion = 'all' ` (intersection of the last two balls and all exclusions)

` intersection = 'last', exclusion = 'random' ` (intersection of the last two balls and  exclusion only of random ball)

` intersection = 'all', exclusion = 'random' ` (all intersections and  one exclusion of random ball)

` intersection = 'random', exclusion = 'random' ` (one intersection and  one exclusion of random balls)

` intersection = 'empty', exclusion = 'empty' ` (PELT-method)

We present a basic use of the main functions of the `FPOPapproximation` package. 

We install the package from Github:

```r
#devtools::install_github("lpishchagina/FPOPapproximation")
library(FPOPapproximation)
```

## The function chpt_rnorm

The `changes_rnorm` is the generation of data (normal distribution) of dimension p with a given values of means and changes.

`n`  is the time series length.

`p`  is the time series dimension.

`chpts` is the changepoint vector that gives the last index of each segment.

The last element of `changes` is always less than to the length of time series.

By default, `changes = NULL` (for the data without changes). 

`means` is the matrix of successive means for the `p`-variate time series.

By default, `means = matrix(0, ncol = 1, nrow = p)` (for the data without changes). 

The length of each matrix row is equal to the length of `changes` plus one.

`noise` is a variance of the time series. By default, `noise =1`.

```r

#parameters

set.seed(13)
N <- 100
Chpt <-50
Means <-  matrix(c(0,1,1,10), nrow = 2)
Noise <- 1
Dim <- 2
Penality <- 2*Dim*log(N)

```

```r

#Data generation

##the time series with one change
time_series1 <- changes_rnorm(p = Dim, n = N, changes = Chpt, means = Means, noise = Noise)

##the time series without changes

time_series2 <- changes_rnorm(p = Dim, n = N, changes = NULL, means = matrix(0, ncol = 1, nrow = Dim), noise = Noise)

```
## The function approx_fpop

The ` approx_fpop ` function returns the result of the segmentation of FPOP-method using the rectangle approximation.

` data ` is the `p`-variate time series (matrix of real numbers with p-rows and n-columns).

` penalty ` is the value of penalty (a non-negative real number  equals to a classic `2*p*(noise^2)*log(n)` (`n` - data length). 

` intersection ` is the type of intersection : `'empty'`, `'all'`, `'last'`, `'random'` or `'sphere'`.

` exclusion ` is the type of intersection : `'empty'`, `'all'`, `'random'` or `'sphere'`.

The following parameter combinations are implemented:
` intersection ='sphere', exclusion ='sphere' `

` intersection ='all', exclusion ='all' ` 

` intersection ='all', exclusion ='empty' ` 

` intersection ='empty', exclusion ='all' ` 

` intersection ='last', exclusion ='all' ` 

` intersection ='last', exclusion ='random' ` 

` intersection ='all', exclusion ='random' ` 

` intersection ='random', exclusion ='random' ` 

` intersection ='empty', exclusion ='empty' ` 

` test_nb_cands `is the logical parameter (if ` test_nb_cands = TRUE `, than the file "Nb_cands.txt" contains the number of change candidates for each operation.

` test_nb_exclus` is the logical parameter (if ` test_nb_exclus = TRUE `, than the file "Nb_cands.txt" contains the number of exclusion for change candidates for each operation. (!!Change!!)

```r

Approx <- list()

Approx[[1]] <- approx_fpop(data = time_series1, penalty = Penality, intersection = 'sphere', exclusion = 'sphere', test_nb_cands = FALSE, test_nb_exclus = FALSE)
Approx[[2]] <- approx_fpop(data = time_series1, penalty = Penality, intersection = 'all', exclusion = 'all', test_nb_cands = FALSE, test_nb_exclus = FALSE)
Approx[[3]] <-approx_fpop(data = time_series1, penalty = Penality, intersection = 'all', exclusion = 'empty', test_nb_cands = FALSE, test_nb_exclus = FALSE)
Approx[[4]] <-approx_fpop(data = time_series1, penalty = Penality, intersection = 'empty', exclusion = 'all', test_nb_cands = FALSE, test_nb_exclus = FALSE)
Approx[[5]] <-approx_fpop(data = time_series1, penalty = Penality, intersection = 'last', exclusion = 'all', test_nb_cands = FALSE, test_nb_exclus = FALSE)
Approx[[6]] <-approx_fpop(data = time_series1, penalty = Penality, intersection = 'last', exclusion = 'random', test_nb_cands = FALSE, test_nb_exclus = FALSE)
Approx[[7]] <-approx_fpop(data = time_series1, penalty = Penality, intersection = 'all', exclusion = 'random', test_nb_cands = FALSE, test_nb_exclus = FALSE)
Approx[[8]] <-approx_fpop(data = time_series1, penalty = Penality, intersection = 'random', exclusion = 'random', test_nb_cands = FALSE, test_nb_exclus = FALSE)
Approx[[9]] <-approx_fpop(data = time_series1, penalty = Penality, intersection = 'empty', exclusion = 'empty', test_nb_cands = FALSE, test_nb_exclus = FALSE)

```

```r
#result for time_series1 the combination ('all','all')
Approx[[2]] 

#$changes
#[1] 50

#$means
#$means[[1]]
#[1] -0.02641869  0.94246956

#$means[[2]]
#[1]  0.902768 10.120444

#$globalCost
[1] 224.9253

```

`changes` is the changepoint vector that gives the last index of each segment.

The last element of `changes` always equals to the length of time series.

`means` is the list of successive means for the p-variate time series.

`globalCost` is the overall Gaussian cost of the segmented data. 

[Back to Top](#top)
