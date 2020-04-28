 DipScale <- function(data){
  # vector to store the dip values
  dips <- vector(mode = "numeric", length = ncol(data))
  
  # compute the maximum dip
  for(i in 1:ncol(data)){
    cur_dip <- dip.test(data[, i])$statistic[[1]]
    dips[i] <- cur_dip
  }
  
  # scale each dimension with it's dip value
  for (i in 1:ncol(data)){
    temp_dat <- data[,i]
    min_val <- min(temp_dat)
    max_val <- max(temp_dat)
    
    # scale the data
    data[, i] <- (data[, i] - min_val)/(max_val - min_val)*dips[i]
  }
  
  return(data)
}

.ComputeDips <- function(data){
  # vector to store the dip values
  dips <- vector(mode = "numeric", length = ncol(data))
  
  # compute the maximum dip
  for(i in 1:ncol(data)){
    cur_dip <- dip.test(data[, i])$statistic[[1]]
    dips[i] <- cur_dip
  }
  return(dips)
}

# function to rotate the data
RotateData <- function(theta, x, y){
  theta <- 0.01745*theta
  sinalpha <- sin(theta);
  cosalpha <- cos(theta);
  
  # rotate the x
  x <- x*cosalpha-y*sinalpha;
  y <- y*sinalpha+y*cosalpha;

  return(data.frame(x = x, y = y))
}

DipTransform <- function(data, c = 5){
  # compute dip values
  dip.max <- max(.ComputeDips(data))
  deg = 0
  k <- ncol(data)
  while(deg < 180*k){
    for(i in seq(1,k-1,1)){
      for(j in seq(i+1,k,1)){
        dip_i <- dip.test(data[, i])$statistic[[1]]
        dip_j <- dip.test(data[, j])$statistic[[1]]
        a <- max(dip_i/dip_j, dip_j/dip_i)
        x <- data[, i]
        y <- data[, j]
        res <- RotateData(c/a, x,y)
        data[, i] <- res$x
        data[, j] <- res$y
        deg <- deg + c/a
        dip.vals <- .ComputeDips(data)
        if(max(dip.vals) > dip.max){
          data <- DipScale(data)
          dip.max <- max(dip.vals)
        }
      }
    }
  }
  return(data)
}

