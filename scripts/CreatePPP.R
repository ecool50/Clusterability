CreatePPP <- function(dat){
  library(sp)
  library(maptools)
  library(spatstat)
  
  pts <- dat
  dimnames(pts)[[1]] = rownames(dat)
  df = data.frame(a = 1:nrow(dat))
  row.names(df) = rownames(pts)
  k <- SpatialPointsDataFrame(pts, df)
  ppp <- as.ppp(k)
  return(ppp)
}