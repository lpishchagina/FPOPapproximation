#' @title chpt_rnorm
#'
#' @description Generation of data (normal distribution) of dimension p with a given values of means and changepoints
#'
#' @param n number of data point.
#' @param p parameter of dimension.
#' @param chpts a vector of increasing changepoint indices (last element is always less then n).
#' By default, 'chpts = NULL' (the data without changepoints).
#' @param means matrix of successive means for data, by default 'means = matrix(0, ncol = 1, nrow = p)'.
#' @param noise standard deviation of an additional normal noise, by default 'noise = 1'.
#'
#' @return matrix of data of dimension p x n with a given values of means by the segmentation.
#'
#' @examples
#' set.seed(13)
#'
#' size <- 1000
#'
#' data1 <- chpt_rnorm(p = 2, n = size, chpts = NULL, means = matrix(0, ncol = 1, nrow = 2), noise = 1)
#'
#' data2 =  chpt_rnorm(p = 2, n = size, chpts = 500, means = matrix(c(0,0,1,1), nrow = 2), noise = 1)
chpt_rnorm <- function(p, n, chpts = NULL, means = matrix(0, ncol = 1, nrow = p), noise = 1){
  #---stop---#
  if (!is.null(chpts) && n <= chpts[length(chpts)]){stop('last element of changepoints is always less than n')}

  if(!is.null(chpts) && !is.numeric(chpts)){stop('chpts are not all numeric')}
  if(is.unsorted(chpts)){stop('chpts should be an increasing vector')}
  if(!is.numeric(means)){stop('means are not all numeric')}

  if ((length(chpts)+1) !=  length(means[1,])){stop('The length of the means[1,] is always equal to the number of changepoints plus one')}
  if ((length(chpts)+1) !=  length(means[2,])){stop('The length of the means[2,] is always equal to the number of changepoints plus one')}

  if(!is.double(noise)){stop('noise is not a double')}
  if(noise < 0){stop('noise must be non-negative')}
  #---function---#
  res <- matrix(0,p,n)
  InttT<- diff(c(0,chpts, n))
  for (i in 1:p){res[i,] <- rep(means[i,], InttT) + rnorm(n, 0, noise)}
  return(res)
}
