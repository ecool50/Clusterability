scatter <- function(v, un, up){
  ( (v - un) %*% t(v - un) ) + ( (v - up) %*% t(v - up) )}

action <- function(uv, m) {
  abs( uv %*% m %*% matrix(uv) ) }

t <- seq(- pi , pi, 0.05)
uv <- cbind(cos(t), sin(t))

distance.from.plane = function(x,w,b) {
 b - sum(x*w)
  }

# classify.fisher = function(x,w,b) {
#  distances = apply(x, 1, distance.from.plane, w, b)
# return(ifelse(distances < 0, 1, 2))}

# data(iris)
ComputeFisher <- function(x_comb, assignments){
  x <- x_comb
  Y <- assignments
  u <- apply(x,2,mean)
  up <- apply(subset(x,Y==1),2,mean)
  un <- apply(subset(x,Y==2),2,mean)
  np <- sum(Y==1)
  nn <- sum(Y==2)
  SB <- nn * (un - u) %*% t(up - u)
  SW <- matrix( apply( apply(x,1,scatter,un=un,up=up), 1, sum ), nrow=2 )
  ratios <- apply(uv, 1, action, SB) / apply(uv, 1, action, SW)
  mr <- which.max(ratios)
  muv <- uv[mr,]
  mv <- 40*ratios[mr]*muv
  xp <- as.vector(x %*% muv)
  rxp <- round(range(xp),0)+c(-1,1)
  xpn <- subset(xp,Y==1)
  xpp <- subset(xp,Y==2)
  b = (mean(xpp) * sd(xpn) + mean(xpn) * sd(xpp)) / (sd(xpp) + sd(xpn))
  # plot(x,col=Y,asp=1)
  # par(lwd=2)
  return(list(intercept = b/muv[2], slope = -muv[1]/muv[2], proj = xp, dat = x_comb, labels = Y))
  # abline(,)
}

