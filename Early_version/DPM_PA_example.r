library(RgeoProfile)
sim <- rDPM(n_traps,alpha=5,sigma=1,tau=3)

trap_locations <- cbind(sim$longitude,sim$latitude)
hits <- sample(0:1, n_traps, replace = TRUE)
trap_data <-  cbind(trap_locations, hits)
colnames(trap_data) <- c("x", "y", "hits")
trap_data


range(trap_data[,1])
range(trap_data[,2])
xminMax <- c(-0.2, -0.07)
yminMax <- c(51.4, 51.6)

anchor_points_x <- seq(xminMax[1], xminMax[2], 0.01)
anchor_points_y <- seq(yminMax[1], yminMax[2], 0.01)

