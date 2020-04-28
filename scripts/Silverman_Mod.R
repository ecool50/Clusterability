
#' Number of modes
#'
#' Calculates the number of modes for given y-values of a density function.
#'
#' @param y vector of y-values of a density function
#'
#' @return The number of modes
#'
#' @export
nr.modes <- function(y) {
  d1 <- diff(y)
  signs <- diff(d1 / abs(d1))
  length(signs[signs == -2])
}

#' Bandwidth calculation for density estimation
#'
#' Calculates the smallest value so that the gaussian kernel density estimate of the given data \code{x} has \code{k} modes.
#' The smaller you choose the bandwidth for a kernel density estimate, the larger the number of modes becomes. This function calculates the smallest value leading to a kernel density estimate with \code{k} number of modes.
#'
#' @param x vector of data
#' @param k number of modes
#' @param prec number of digits for precision of calculation
#' @param density.fun A function that returns a vector of density estimates
#'
#' @return the smallest value so that the gaussian kernel density estimate of the given data \code{x} has \code{k} modes.
#'
#' @examples h.crit(rnorm(10), k = 1)
#'
#' @export
h.crit <- function(x, k, prec = 6, density.fun = NULL) {
  if (is.null(density.fun)) {
    density.fun <- function(x, h) {
      density(x, bw = h, kernel = "gaussian")$y
      # fkde::fkde(x, h = h)$y
    }
  }
  
  digits <- prec
  prec <- 10 ^ (-prec)
  x <- sort(x)
  minh <- min(diff(x))		#minimal possible h
  maxh <- diff(range(x)) / 2	#maximal possible h
  a <- maxh
  b <- minh
  zaehler <- 0
  
  while (abs(b - a) > prec) {
    message(abs(b - a))
    m <- nr.modes(density.fun(x, a))
    
    b <- a
    if (m > k) {
      minh <- a
      a <- (a + maxh) / 2
    }
    else {
      maxh <- a
      a <- (a - minh) / 2
    }
  }
  a <- round(a, digits)
  
  if (nr.modes(density.fun(x, a)) <= k) {
    #subtract until more than k modes
    while (nr.modes(density.fun(x, a)) <= k) {
      a <- a - prec
    }
    a <- a + prec
  }
  if (nr.modes(density.fun(x, a)) > k) {
    #add until nr. of moodes correct
    while (nr.modes(density.fun(x, a)) > k) {
      a <- a + prec
    }
  }
  a
}


#' Silvermantest
#'
#' The silvermantest tests the null hypothesis that an underlying density has at most \code{k} modes.
#'
#' @param x vector of data
#' @param k number of modes for the null hypothesis
#' @param R number of bootstrap replications
#' @param adjust boolean to activate the adjusting of the p-value (valid if k=1) (see Hall and York)
#' @param digits number of digits of the p-value
#' @param density.fun A function that returns a vector of density estimates
#'
#' @return An object of the class Silvermantest (see: \code{\link{Silvermantest-class}}).
#'
#' @examples
#' x <- c(rpois(n = 50, lambda = 1), rnorm(n = 100, mean = 4))
#' silverman.test(x, k = 1)
#'
#' @importFrom stats sd
#' @importFrom stats rnorm
#' @importFrom stats density
#' @importFrom stats predict
#' @importFrom methods new
#' @export
silverman.test <- function(x, k = 1, R=100, adjust=FALSE, digits=6, density.fun=NULL){
  
  library(boot)
  library(parallel)
  # x: data
  # k: number of modes to be tested
  # M: number of bootstrap replications
  
  #check if seed is available (as done in boot package)
  #if so save it
  seedAvailable = exists(x=".Random.seed", envir=.GlobalEnv, inherits=FALSE)
  if(seedAvailable)
    saved_seed <- .Random.seed
  else{
    rnorm(1)
    saved_seed <- .Random.seed
  }
  
  #temp function for bootstrapping
  y.obs <- function(x,h,sig=sd(x)){
    mean(x) + (x-mean(x)+h*rnorm(length(x),0,1))/((1+h^2/sig^2)^(1/2))
    #(x+h*rnorm(length(x),0,1))/((1+h^2/sig^2)^(1/2))
  }
  
  #temp function for density calculation
  if(is.null(density.fun)) {
    density.fun <- function(x,h){density(x, bw = h, kernel = "gaussian")$y}
  }
  
  #start of the test
  message("Computing critical bandwidth")
  # h0 <- h.crit(x, k, density.fun=density.fun)
  h0 <- multimode::bw.crit(x, mod0 = k)
  
  # statistic function
  mode.fun <- function(d, i, h0) {
    x.boot <- sort(y.obs(d[i], h0))
    nr.modes(density.fun(x.boot, h0))
  }
  
  cl <- makePSOCKcluster(detectCores() - 1)
  clusterExport(cl, 'mode.fun')
  message("Running bootstrap")
  mod.boot <- boot::boot(x, statistic = mode.fun, R = R, h0 = h0, parallel = "multicore", ncpus = 6, cl = cl)
  
  n <- sum(as.vector(mod.boot$t) > k)
  p <- n/R
  
  if (adjust) {
    if (k==1) {
      #asymptotic levels of silvermantest by Hall/York
      x=c(0,0.005,0.010,0.020,0.030,0.040,0.050,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.25,0.30,0.35,0.40,0.50)
      y=c(0,0,0,0.002,0.004,0.006,0.010,0.012,0.016,0.021,0.025,0.032,0.038,0.043,0.050,0.057,0.062,0.07,0.079,0.088,0.094,0.102,0.149,0.202,0.252,0.308,0.423)
      sp = splines::interpSpline(x,y)
      #adjusting the p-value
      if (p<0.005)
        p <- 0
      else{
        p <- predict(sp,p)$y
        p <- round(p,digits)
      }
    } else{
      warning("The option to adjust the p-value is valid only for k=1")
    }
  }
  test_obj <- list(data=x, p_value = p, saved_seed=saved_seed, k=k)
  return(test_obj)
}