#' SmoothData
#'
#' @description Smoothing technique. See Olshen and Venkatraman "A faster circular binary segmentation algorithm for the analysis of array CGH data".
#'
#' @param X data matrix or data vector
#' @param R radius for smoothing
#' @param L threshold for correction
#' @param M correction towards median
#' @param times number of iterations of the smoothing procedure
#'
#' @return smoothed data in the same form as the input (matrix or vector)
#' @export
#'
#' @examples set.seed(1)
#' X <- rt(100,3)
#' X_smooth <- SmoothData(X)
#'
SmoothData <- function(X, R=3, L=4, M=2, times=1){

  if(is.matrix(X)){
    n <- nrow(X)
    for(i in 1:n)
      X[i,] <- SmoothData(X[i,], R=R, L=L, M=M, times=times)
    return(X)
  }
  sigma <- sd(X)
  T_ <- length(X)
  X_smooth <- X
  for(i in (R+1):(T_-R)){
    region <- X[(i-R):(i+R)]
    ord <- order(region)
    if(region[ord[2]] - region[ord[1]]> L*sigma)
      X_smooth[i-R-1+ord[1]] <- region[ord[R+1]] - M*sigma
    if(region[ord[2*R+1]] - region[ord[2*R]]> L*sigma)
      X_smooth[i-R-1+ord[2*R+1]] <- region[ord[R+1]] + M*sigma
  }
  if(times>1)
    return(SmoothData(X_smooth, R=R, L=L, M=M, times=times-1))
  X_smooth
}
