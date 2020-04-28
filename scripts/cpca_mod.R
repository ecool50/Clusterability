#-------------------------------------------------------------------------------
# `cpca_stepwise_eigen` is an updated version of `cpc_stepwise_base`.
# - the input covariance matrices are of one of `Matrix` types;
# - the property, that the covariance matrices are symmetric, is used for speed up.
# 
# Note:
# - `%*%` operator is replaced by more efficient ones, e.g. `crossprod`;
#-------------------------------------------------------------------------------

#' @importFrom Matrix Matrix isSymmetric crossprod tcrossprod
#' @importFrom plyr laply
.cpca_stepwise_eigen <- function(cov, ng, ncomp = 0, 
                                tol = 1e-6, maxit = 1e3,
                                start = c("eigen", "random"), symmetric = TRUE,
                                verbose = 0, ...)
{
  ### par
  stopifnot(length(cov) == length(ng))
  
  cl <- plyr::laply(cov, function(x) as.character(class(x)))
  stopifnot(all(grepl("*matrix$", cl)))
  
  stopifnot(symmetric)
  covSymm <- plyr::laply(cov, Matrix::isSymmetric)
  stopifnot(all(covSymm))
  
  start <- match.arg(start)
  
  ng <- as.numeric(ng) # needed for further multiplication with double precision like `ng / n` 
  
  ### var
  k <- length(cov) # number of groups `k`
  p <- nrow(cov[[1]]) # number of variables `p` (the covariance matrices p x p)
  n <- sum(ng) # the number of samples
  
  # If ncomp = 0 retrieve all components
  if(ncomp == 0) {
    ncomp <- p
  }
  
  ### output variables
  D <- matrix(0, nrow = ncomp, ncol = k)
  CPC <- matrix(0, nrow = p, ncol = ncomp)
  
  Qw <- Matrix::Matrix(diag(1, p))
  
  convergedComp <- rep(FALSE, ncomp)
  itComp <- rep(0, ncomp)
  
  ### step 1: compute the staring estimation
  if(start == "eigen") {
    S <- Matrix::Matrix(0, nrow = p, ncol = p)
    stopifnot(Matrix::isSymmetric(S))
    
    # val <- ng[1] / n
    # for(i in 1:k) {
    #   S <- S +  (val * cov[[i]])
    #   gc()
    # }
    
    S <- Reduce('+', cov)*(1/3)
    
    res <- RSpectra::eigs_sym(as.matrix(S), k = ncomp)
    gc()
      # eigen(S, symmetric = symmetric) # `?eigen`: a vector containing the p eigenvalues of ‘x’, sorted in _decreasing_ order
    # all(order(res$values[1:ncomp], decreasing = TRUE) == seq(1, ncomp))
    
    q0 <- Matrix::Matrix(res$vectors[, 1:ncomp, drop = FALSE])
  } else {
    q0 <- Matrix::Matrix(runif(p * ncomp), nrow = p, ncol = ncomp)
  }
  
  ### step 2: estimation of `p` components in a loop
  for(comp in 1:ncomp) {
    if(verbose > 1) {
      cat(" * component:", comp, "/", ncomp, "\n")
    }

    q <- q0[, comp]

    d <- rep(0, k)
    for(i in 1:k) {
      d[i] <- as.numeric(Matrix::crossprod(q, Matrix::crossprod(Matrix(cov[[i]], sparse = T), Matrix(q, sparse = T))))
      gc()
    }

    # loop along `it`
    cost0 <- 0
    for(it in 1:maxit) {
      if(verbose > 1) {
        cat(" * it:", it, "/", maxit, "\n")
      }

      S <- Matrix::Matrix(0, nrow = p, ncol = p)
      for(i in 1:k) {
        S <- S + (ng[i] / d[i]) * cov[[i]]
        gc()
      }

      w <- Matrix::crossprod(Matrix(S, sparse = T), Matrix(q, sparse = T)) # (1) S %*% q; (2) tcrossprod(S, t(q)), as S is symmetric
      gc()
      if(comp != 1) {
        w <- Matrix::crossprod(Matrix(Qw, sparse = T), Matrix(w, sparse = T)) # (1) Qw %*% w; (2) tcrossprod(Qw, t(w)), as Qw is symmetric
        gc()
      }

      q <- as.numeric(w) / sqrt(as.numeric(Matrix::crossprod(Matrix(w, sparse = T)))) # normalize
      gc()

      # compute `cost` & `d`
      for(i in 1:k) {
        d[i] <- as.numeric(Matrix::crossprod(Matrix(q,sparse = T), Matrix::crossprod(Matrix(cov[[i]],sparse = T), Matrix(q, sparse = T))))
        gc()
      }

      cost <- sum(log(d) * ng)

      delta <- abs((cost - cost0) / cost)
      if(verbose > 1) {
        cat("  --  delta:", delta, "\n")
      }
      if(delta <= tol) {
        break
      }

      cost0 <- cost
    }
    # end of loop along `it`
    itComp[comp] <- it
    convergedComp[comp] <- (it < maxit)

    D[comp, ] <- d
    CPC[, comp] <- q

    Qw <- Qw - Matrix::tcrossprod(Matrix(q, sparse = T))
    gc()
  }
  # end of loop along `comp`

  ### return
  out <- list(D = D, CPC = CPC, ncomp = ncomp,
              convergedComp = convergedComp, converged = all(convergedComp),
              itComp = itComp, maxit = maxit)

  oldClass(out) <- c("CPCAPowerEigen", "CPCAPower")
  
  return(out)
}
