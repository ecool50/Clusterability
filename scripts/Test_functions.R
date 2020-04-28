clusters <- hclust(dist(iris[, 1:4]))
plot(clusters)


idx_hc <- function(hc, n) {
  ## list array of cluster indices at each of the n-1 merges
  idx_hc <- array(list(), c(2*n-1, 2))
  idx_hc[1:n, 1] <- as.list(n:1)
  idx_hc[(n+1):(2*n-1), ] <- hc$merge + n + (hc$merge<0)
  
  ## complete idx_hc
  for (k in 1:(n-1)) {
    idx_hc[[n+k, 1]] <- unlist(idx_hc[idx_hc[[n+k, 1]], ])
    idx_hc[[n+k, 2]] <- unlist(idx_hc[idx_hc[[n+k, 2]], ])
  }
  
  ## flip index, hclust and shc use revered ordering
  idx_hc[(2*n-1):(n+1), ]
}

n <- nrow(iris)

temp <- idx_hc(clusters, n)

groups <- cutree(clusters, k = 3)

temp1 <- rbind(iris[temp[[1]], 1:4], iris[temp[[n]], 1:4])

dip.test(parallelDist(as.matrix(temp1)))$p.value

for(i in 1:10){
  full_dat <- rbind(iris[temp[[i]], 1:4], iris[temp[[i+n-1]], 1:4])
  print(dip.test(parallelDist(as.matrix(full_dat)))$p.value)
}
