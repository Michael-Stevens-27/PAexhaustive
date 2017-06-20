############################################################
grid_x <- seq(-1,1,1)
grid_y <- seq(-1,1,1)
gridded <- expand.grid(grid_x,grid_y)
r <- 5
theta <- seq(0,2*pi, 1)
circle_x <- r*cos(theta)
circle_y <- r*sin(theta)
circle <- cbind(circle_x, circle_y)
colnames(circle) <- c("Var1", "Var2")
k <- rbind(gridded, circle)
dummy_hits <- rep(0, length(gridded[,1]))
dummy_miss <- rep(1, length(circle[,1]))
dummy_hm <- c(dummy_hits, dummy_miss)
outer_dummy_data <- cbind(k, dummy_hm)
############################################################
line_hx <- seq(-2, 2, 0.5)
line_hy <- rep(2.456, length(line_hx))
line_mx <- seq(-2,2, 0.5)
line_my <- rep(1.456, length(line_mx))

line_hits <- cbind(line_hx, line_hy)
line_miss <- cbind(line_mx, line_my)
colnames(line_hits) <- c("x", "y")
colnames(line_miss) <- c("x", "y")

lines <- rbind(line_hits, line_miss)
line_one <- rep(1, length(line_hits[,1]))
line_zero <- rep(0, length(line_miss[,1]))

line_onezero <- c(line_one, line_zero)
lines_dummy_data <- cbind(lines, line_onezero)
############################################################
cross_hx <- seq(-2,2,0.5)
cross_hy <- seq(-2,2,0.5)
cross_mx <- seq(-2,2,0.5)
cross_my <- seq(2,-2,-0.5)
cross_hxy <- cbind(cross_hx, cross_hy)
cross_mxy <- cbind(cross_mx, cross_my)
colnames(cross_mxy) <- c("x","y")
colnames(cross_hxy) <- c("x","y")
cross <- rbind(cross_hxy, cross_mxy)
cross_one <- rep(1, length(cross_hxy[,1]))
cross_zero <- rep(0, length(cross_mxy[,1]))
cross_01 <- c(cross_one, cross_zero)
cross_data <- cbind(cross, cross_01)
############################################################



## Extract ALL Parameters from your data
Data_params <- Extract_Params(cross_data, x_grid_cells =15, y_grid_cells = 15, Guard_Rail = 1, Trap_Radius = 0.9)
Trap_Po_Params <- Trap_Po_Parameters(Data_params)

##_____________COMPUTE PROBABILITY MATRICES______________##

## Return the probability of seeing your data from a SINGLE SOURCE model
Single_Source_Prob <- Single_Source_Probability(Data_params, Trap_Po_Params)

## Single source normalising
### NORMALISE ??????
Single_Source_Prob$SS_Hits <- Single_Source_Prob$SS_Hits/sum(Single_Source_Prob$SS_Hits)
Single_Source_Prob$SS_Miss <- Single_Source_Prob$SS_Miss/sum(Single_Source_Prob$SS_Miss)
Single_Source_Prob$SS_Both <- Single_Source_Prob$SS_Both/sum(Single_Source_Prob$SS_Both)

png("cross.png")
## Return the probability of seeing your data from a Two Source Model
Two_Source_Prob <-DS_Po_Params(Data_params, Trap_Po_Params)
#___________________________________
##_____________PLOTS______________##

### Single Source

#plot1source(Data_params, Single_Source_Prob)

### Two Source

plot2source(Data_params, Two_Source_Prob)
dev.off()

getwd()
