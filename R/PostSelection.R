#' RemoveSmallJumps
#' @description Removes all jumps with change in mean less than \code{threshold*SD}
#' where \code{SD} is a measurement of the standard deviation using the median
#' absolut deviation estimator.
#' @param X data matrix or data vector
#' @param cps ordered list of change points the should be reconsidered
#' @param threshold jump height relative to standard deviation
#'  so that the change point is accepted
#' @return integer vector with change points with little change removed
#' @export
#' @examples set.seed(4)
#' X <- matrix(rep(c(1,0,1,0), c(80,30,50,40))+rnorm(200,0,0.6), 2, byrow = TRUE)
#' cps <- optimisticCBS(X, alpha = 0.2)
#' cps
#' RemoveSmallJumps(X, cps, 1)
RemoveSmallJumps <- function(X, cps, threshold = 1){
  if(length(cps)==0) return(cps)
  if(is.vector(X)) X <- matrix(X,1)
  T_ <- ncol(X)
  n <- nrow(X)
  for(i in 1:n){
    SD <- mad(diff(X[i,]))/sqrt(2)
    if(SD>1e-6)
      X[i,] <- (X[i,]-mean(X[i,]))/SD
  }
  if(cps[1]!=0)
    cps <- c(0, cps, T_)
  k <- length(cps)-2
  means <- matrix(0,nrow=k+1, ncol=n)
  for(i in 1:n){
    for(l in 1:(k+1)){
      means[l, i] <- mean(X[i,(cps[l]+1):cps[l+1]])
    }
  }
  # returns in which dimension the jumps are largest
  jumps <- abs(diff(means))
  max_i <- max.col(jumps)
  remove <- c(1,k+2)
  for(l in 1:k){
    if(jumps[l,max_i[l]] < threshold)
      remove <- c(remove, l+1)
  }
  cps <- cps[-remove]
  if(length(remove)>2)
    return(RemoveSmallJumps(X, cps))
  cps
}
