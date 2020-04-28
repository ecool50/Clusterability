# create


# Minimum Density Hyperplane

X <- t(obj_sub)
sol <- mdh(X)
plot(sol)
sol_dist <- parallelDist(sol$fitted)
dip_stat <- dip.test(sol_dist)$p.value

plot(sol$fitted[, 1], sol$fitted[, 2], xlab = "Plane 1", ylab = "Plane 2", main = "Minimum density hyperplane data projection")
hist(sol_dist, breaks = 300, main = paste("Hash noCycle noMyeloid noTfh noTh2 \ndensity plot on Min Density Hyperplane", "\nDip Test p-value =", round(dip_stat,4)), xlab = "Distance", probability = T, col = "lightblue")
lines(density(sol_dist), col="red", lwd=2)
