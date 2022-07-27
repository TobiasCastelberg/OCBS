## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib OCBS, .registration = TRUE
#' @exportPattern ^[[:alpha:]]+
## usethis namespace: end
NULL

#' Optimistic Circular Binary Segmentation
#' @description This is the main function of the package.
#' @param X data vector
#' @param optimization minimization of "L1" or "L2" loss
#' @param method "advanced" or "full" search
#' @param select Model selection approach "perm" (permutation test only), "mixed" (first Wilcoxon test, then permutation test), "wilc" (Wilcoxon test only)
#' @param alpha Significance level for testing change points
#' @param nr_perms Number of permutation for permutation tests
#' @param bd stopping boundary for early stopping
#' @param min_seg minimal segment length. Note that change points still can be closer together than \code{min_seg}
#'
#' @return Integer vector containing the change points found in the data
#' @export
#'
#' @examples set.seed(1)
#' X <- rep(c(0,1,0),c(100,20,80)) + rnorm(200,0,0.6)
#' optimisticCBS(X)
optimisticCBS <- function(X, optimization="L2", method="advanced", select="mixed",
                          alpha = 0.01, nr_perms = 10000, bd = NULL, min_seg = 2){
  if(is.vector(X)) X <- matrix(X,1)
  T_ <- ncol(X)
  n <- nrow(X)
  cps <- numeric()
  if(T_<3*min_seg)
    # Segment too short for (further) segments
    return(cps)

  # Standardize data to avoid overflow of cumulative sums
  for(i in 1:n){
    SD <- mad(diff(X[i,]))/sqrt(2)
    if(SD>1e-6)
      X[i,] <- X[i,]/SD

    # Add a bit of noise to avoid equal ranks (this randomizes Wilcoxon tests)
    X[i,] <- X[i,] + rnorm(T_, 0, 1e-8)
  }

  # Calculate stopping boundary
  if(is.null(bd))
    bd <- StoppingBoundary(nr_perms, alpha, eta_star = alpha * 0.05)

  # First check on CBS statistics
  circ_test <- CheckValid(X, bd, optimization, method, select, alpha, nr_perms, circular = T, min_seg)
  if(circ_test$cps_valid){
    s <- circ_test$s
    t <- circ_test$t

    # Second check on BS statistics
    s_valid = F
    t_valid = F
    if(s != 0 & t!=T_){
      X_left = matrix(X[,1:t], n)

      left_test <- CheckValid(X_left, bd, optimization, method, select, alpha, nr_perms, circular = F, min_seg=2)
      s_valid <- left_test$cps_valid
      X_right = matrix(X[,(s+1):T_], n)
      right_test <- CheckValid(X_right, bd, optimization, method, select, alpha, nr_perms, circular = F, min_seg=2)
      t_valid <- right_test$cps_valid
    }else{
      # If s=0 or t=T no boundary tests necessary
      if(s==0)
        t_valid = T
      if(t==T_)
        s_valid = T
    }

    # Add changepoints if valid and search on subsegments
    if(s_valid){
      # s valid
      cps_left <- optimisticCBS(X[,1:s], optimization, method, select, alpha, nr_perms, bd, min_seg)
      cps <- c(cps, s, cps_left)
      if(t_valid){
        #s and t valid
        cps_mid <- s + optimisticCBS(X[,(s+1):t], optimization, method, select, alpha, nr_perms, bd, min_seg)
        cps_right <- t + optimisticCBS(X[,(t+1):T_], optimization, method, select, alpha, nr_perms, bd, min_seg)
        cps <- c(cps, t, cps_mid, cps_right)
      }else{
        # s valid, t invalid
        cps_mid_right <- s + optimisticCBS(X[,(s+1):T_], optimization, method, select, alpha, nr_perms, bd, min_seg)
        cps <- c(cps, cps_mid_right)
      }
    }else if(t_valid){
      # s invalid, t valid
      cps_left_mid <- optimisticCBS(X[,1:t], optimization, method, select, alpha, nr_perms, bd, min_seg)
      cps_right <- t + optimisticCBS(X[,(t+1):T_], optimization, method, select, alpha, nr_perms, bd, min_seg)
      cps <- c(cps, t, cps_left_mid, cps_right)
    }
  }
  return(sort(cps))
}


#' Candidate testing
#'
#' @param X data matrix
#' @param bd stopping boundary for early stopping
#' @param optimization minimization of "L1" or "L2" loss
#' @param method "advanced" or "full" search
#' @param select Model selection approach "perm" (permutation test only), "mixed"
#' (first Wilcoxon test, then permutation test), "wilc" (Wilcoxon test only)
#' @param alpha Significance level for testing change points
#' @param nr_perms Number of permutation for permutation tests
#' @param circular \code{circular = T} does CBS, else the test will be for BS
#' @param min_seg minimal segment length. Note that change points still can be
#' closer together than \code{min_seg}
#'
#' @return list containing whether the change points are valid or not and the
#' corresponding candidates. For (non-circular) BS, there is only one potential
#' change point, which is stored in \code{t}. In that case \code{s=0.}
#'
#' @export
#'
#' @examples set.seed(1)
#' X <- matrix(rep(c(0,1,0),c(100,20,80)) + rnorm(200,0,0.6), nrow=1)
#' CheckValid(X, StoppingBoundary())
#' @examples set.seed(1)
#' X <- matrix(rnorm(200,0,0.6), nrow=1)
#' CheckValid(X, StoppingBoundary())
## usethis namespace: start
#' @importFrom stats mad rnorm sd wilcox.test
## usethis namespace: end
CheckValid <- function(X, bd = NULL, optimization = "L2", method = "advanced", select ="mixed",
                       alpha = 0.01, nr_perms = 10000, circular = T, min_seg = 2){
  T_ <- ncol(X)
  n <- nrow(X)
  if(T_>20 & select != "perm"){
    # WILCOXON APPROACH
    even <- seq(2,T_,2)
    odd <- seq(1,T_,2)
    X_batch1 <- matrix(X[,odd],n)
    X_batch2 <- matrix(X[,even],n)
    # Find candidates
    cand1 <- MaxStats(X_batch1, optimization, method, circular, max(2,floor(min_seg/2)))
    cand2 <- MaxStats(X_batch2, optimization, method, circular, max(2,floor(min_seg/2)))
    t1 <- cand1$ind
    t2 <- cand2$ind
    if(circular){
      s1 <- cand1$shift
      s2 <- cand2$shift
    }
    else{
      s1 <- 0
      s2 <- 0
    }

    # Wilcoxon test
    if(is.element(s1, c(0,1)+s2) & is.element(t1, c(0,1)+t2) & T_ > 200){
      # No sample splitting increases power in trade for significance
      # Only done if signal on both batches coincides and batch sizes are big enough
      # to reduce probability that they coincide under null hypothesis (i.e. no signal)
      pnWilc <- rep(0,n)
      for(i in 1:n){
        cand <- MaxStats(X, optimization, method, circular, min_seg)
        s <- cand$shift
        t <- cand$ind
        pnWilc[i] <- wilcox.test(x = X[(s+1):t], y = X[-c((s+1):t)],
                                 alternative = "two.sided")$p.value
      }
      # Geometric mean if time series has several dimensions (?)
      pWilc <- exp(mean(log(pnWilc)))
    }else{
      # Sample splitting (different candidates on the two batches)
      pnWilc1 <- rep(0,n)
      pnWilc2 <- rep(0,n)
      for(i in 1:n){
        x <- X_batch2[i,(s1+1):(t1-1)]
        y <- X_batch2[i,-c((s1+1):(t1-1))]
        if(mean(X_batch1[i,(s1+1):t1])> mean(X_batch1[i,-c((s1+1):t1)]))
          alternative <- "greater"
        else
          alternative <- "less"
        pnWilc1[i] <- wilcox.test(x, y, alternative = alternative)$p.value

        x <- X_batch1[i,(s2+2):t2]
        y <- X_batch1[i,-c((s2+2):t2)]
        if(mean(X_batch1[i,(s1+1):t1])> mean(X_batch1[i,-c((s1+1):t1)]))
          alternative <- "greater"
        else
          alternative <- "less"
        pnWilc2[i] <- wilcox.test(x, y, alternative = alternative)$p.value
      }
      # Geometric means
      pWilc1 <- exp(mean(log(pnWilc1)))
      pWilc2 <- exp(mean(log(pnWilc2)))
      # Check which batch performed better
      if(pWilc1 < pWilc2)
        cand <- correction(X, s1, t1, optimization, batch = 1)
      else
        cand <- correction(X, s2, t2, optimization, batch = 2)
      s <- cand$s
      t <- cand$t
      # Correction to multiple independent testing
      pWilc <- 1- (1-min(pWilc1,pWilc2))^2
    }

    # If Wilcoxon test rejects the Null, we return
    if(pWilc < alpha)
      return(list("cps_valid"=T, "s"=s, "t"=t))
    if(select == "wilc")
      return(list("cps_valid"=F, "s"=s, "t"=t))
  }


  # PERMUTATION APPROACH
  cand <- MaxStats(X, optimization, method, circular, min_seg)
  cps_valid <- PermTest(X, bd, cand$stat, optimization, method, alpha, nr_perms, circular, min_seg)
  return(list("cps_valid"=cps_valid, "s"=cand$shift, "t"=cand$ind))
}


.correction <- function(X, s, t, optimization, batch){
  # Optimizes locations after using sample splitting
  T_ <- ncol(X)
  if(batch==1){
    t <- 2*t-1
    if(s>0)
      s <- 2*s-1
  }
  else{
    s <- 2*s
    t <- 2*t
  }
  stat <- EvalStatSlow(X, s, t, optimization)
  if(s>0){
    stat_s <- EvalStatSlow(X, s+1, t, optimization)
    if(stat_s>stat){
      stat <- stat_s
      s <- s+1
    }
  }
  if(t<T_){
    stat_t <- EvalStatSlow(X, s, t+1, optimization)
    if(stat_t > stat){
      stat <- stat_t
      t <- t+1
    }
  }
  list("s"=s, "t"=t, "stat"=stat)
}

