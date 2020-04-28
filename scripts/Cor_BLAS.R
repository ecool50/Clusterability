#' Correlation calculation. use BLAS and data.table to speed up.
#'
#' @importFrom data.table frank
#' @importFrom RhpcBLASctl omp_get_num_procs omp_set_num_threads
#' @param x matrix; input data, rows for variable (genes), columns for observations (cells).
#' @param y matrix; input data, rows for variable (genes), columns for observations (cells) (default: NULL)
#' @param method character; method used. (default: "pearson")
#' @param nthreads integer; number of threads to use. if NULL, automatically detect the number. (default: NULL)
#' @param na.rm logical; remove missing values. (default: T)
#' @details calcualte the correlation among variables(rows)
#' @return correlation coefficient matrix among rows
#' @export
cor.BLAS <- function(x,y=NULL,method="pearson",nthreads=NULL,na.rm=T)
{
  if(is.null(nthreads))
  {
    nprocs <- RhpcBLASctl::omp_get_num_procs()
    RhpcBLASctl::omp_set_num_threads(max(nprocs-1,1))
  }else{
    RhpcBLASctl::omp_set_num_threads(nthreads)
  }
  cor.pearson <- function(x,y=NULL)
  {
    if(is.null(y)){
      ### x = x - rowMeans(t(na.omit(t(x))))
      x = x - rowMeans(x,na.rm=na.rm)
      x = x / sqrt(rowSums(x^2,na.rm=na.rm))
      ### cause 'memory not mapped' :( ; and slower in my evaluation: 38 sec .vs. 12 sec.
      #x.cor = tcrossprod(x)
      if(na.rm){ x[is.na(x)] <- 0 }
      x.cor = x %*% t(x)
      return(x.cor)
    }else{
      x = x - rowMeans(x,na.rm=na.rm)
      x = x / sqrt(rowSums(x^2,na.rm=na.rm))
      y = y - rowMeans(y,na.rm=na.rm)
      y = y / sqrt(rowSums(y^2,na.rm=na.rm))
      #xy.cor <- tcrossprod(x,y)
      if(na.rm){ x[is.na(x)] <- 0 }
      if(na.rm){ y[is.na(y)] <- 0 }
      xy.cor <- x %*% t(y)
      return(xy.cor)
    }
  }
  x <- as.matrix(x)
  if(!is.matrix(x)){
    warning("x is not like a matrix")
    return(NULL)
  }
  if(!is.null(y)){
    y <- as.matrix(y)
    if(!is.matrix(y)){
      warning("y is not like a matrix")
      return(NULL)
    }
  }
  if(method=="pearson"){
    return(cor.pearson(x,y))
  }else if(method=="spearman"){
    if(is.null(y)){
      return(cor.pearson(t(apply(x, 1, data.table::frank, na.last="keep"))))
    }else{
      return(cor.pearson(t(apply(x, 1, data.table::frank, na.last="keep")),
                         t(apply(y, 1, data.table::frank, na.last="keep"))))
    }
  }else{
    warning("method must be pearson or spearman")
    return(NULL)
  }
}