CI<-function(x,cluster){
  wcss.feature <- numeric(ncol(x))
  for (k in unique(cluster)) {
    indices <- (cluster == k)
    if (sum(indices) > 1)
      wcss.feature <- wcss.feature + apply(scale(x[indices,], center = TRUE, scale = FALSE)^2, 2, sum)
  }
  wcss = sum(wcss.feature)
  Tss=sum(apply(scale(x, center = TRUE, scale = FALSE)^2, 2, sum))
  wcss/Tss
}


UNPaC_null_sig<- function(x, k,cluster.fun,nsim=100, rho=0.02) {
  ## This version generates a new U every simulation
  ## Also does not center features prior to clustering
  out <- NA
  # if (cov=="glasso"){cur.scov=glasso::glasso(var(x), rho=rho)$w} else {cur.scov=coop::covar(t(x))}
  cur.scov <- coop::covar(t(x))
  cur.iscov <- chol(cur.scov)
  cur.h1 <- apply(x, 2, function(x) multimode::bw.critx)
  x.var <- apply(x, 2, var)
  for (i in 1:nsim) {
    future::plan(multiprocess)
    cur.ystar <- future.apply::future_apply(x, 2, sample, replace=TRUE)
    
    xn <- scale(cur.ystar, scale=FALSE) +
      t(cur.h1*matrix(rnorm(nrow(x)*ncol(x)), ncol=nrow(x), nrow=ncol(x)))
    
    xn <- t((1/sqrt(1+cur.h1^2/x.var))*t(xn))
    xn <- t(colMeans(cur.ystar)+t(xn))
    
    U <- pnorm(matrix(rnorm(ncol(x)*row(x)), ncol = ncol(x)) %*% cur.iscov)
    S <- future.apply::future_sapply(1:ncol(x),function(i) quantile(xn[,i],U[,i]))
    xn.km <- cluster.fun(S, k)$cluster
    out[i] <- CI(S,xn.km)
  }
  return(out)
}



UNPaC_Copula <- function(x,cluster,cluster.fun, nsim=100) {
    test.km.ci <- CI(x,cluster)
    
    out=UNPaC_null_sig(x, k,cluster.fun,nsim=nsim, rho=rho, cov=cov,center=center,scale=scale)
    
    pvalue_emp=sum(test.km.ci>out) / length(out)
    pvalue_norm=pnorm((test.km.ci-mean(out))/sd(out))
    return(list(sim_CI=out, pvalue_emp=pvalue_emp, pvalue_norm=pvalue_norm))
  }