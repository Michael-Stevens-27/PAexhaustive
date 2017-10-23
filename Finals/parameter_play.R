################################################################################
############################## HIT MISS RATIO ##################################
################################################################################

n_offenders <- seq(10,50, 1)
n_traps_x <-  seq(5, 15, 1)
n_traps_y <- seq(5, 15, 1)
repeats <- 200
average_ratio_mat <- matrix(NA, ncol = 3, nrow = length(n_offenders)*length(n_traps_x))
colnames(average_ratio_mat) <- c("n_offenders", "n_traps", "average ratio of hits to misses")
k <- 1
while(k <= length(n_traps_x))
  {
    j <- 1
    while(j <= length(n_offenders))
    {
      i <- 1
      ratio_matrix <- matrix(NA, ncol = 3, nrow = repeats)
      colnames(ratio_matrix) <- c("hits", "misses", "ratio")
      while(i <= repeats)
      {
        abc <- rpois(1,n_offenders[j])
        if(abc > 1)
        {
        sim <- rDPM(abc, alpha = 0, sigma = 1)
        point_data <- geoData(sim$longitude, sim$latitude)

        master_params <- geoParams(data = point_data, sigma_mean = 1, sigma_squared_shape = 2, samples= 50000, chains = 200, burnin = 2000, priorMean_longitude = mean(point_data$longitude), priorMean_latitude = mean(point_data$latitude), guardRail = 0.05)
        Trap_Data <- trap_assignment_data(sim, master_params, n_traps_x =  n_traps_x[k], n_traps_y = n_traps_y[k], "title")
        detection_TR <- Trap_Data$detection_TR
        Trap_Data <- Trap_Data$Trap_data

        ### DPM
        #dpm_hits <- subset(Trap_Data, Trap_Data[,3]>0)
        #dpm_hits <- dpm_hits[rep(1:nrow(dpm_hits), dpm_hits[,3]),]
        #dpm_misses <- subset(Trap_Data, Trap_Data[,3]==0)
        #trap_loc_data <- geoData(Trap_Data[,1], Trap_Data[,2])
        #s <- geoDataSource(sim$source_lon, sim$source_lat)
        #hit_data <- geoData(dpm_hits[,1], dpm_hits[,2])
        #hit_params <- geoParams(data = hit_data, sigma_mean = 1, sigma_squared_shape = 2, samples= 50000, chains = 200, burnin = 2000, priorMean_longitude = mean(hit_data$longitude), priorMean_latitude = mean(hit_data$latitude), guardRail = 0.05)
        #hit_params$output$longitude_minMax <- master_params$output$longitude_minMax
        #hit_params$output$latitude_minMax <- master_params$output$latitude_minMax
        #m <- geoMCMC(data=hit_data, params= hit_params)

        ### PA
        #Trap_Data <- as.data.frame(Trap_Data)
        #Data_parameters <- Extract_Params(sim, Trap_Data, PA_x_grid_cells = 10, PA_y_grid_cells = 10, Guard_Rail = 0.05, Trap_Radius = detection_TR, n_sources = 1, n_cores = 1, n_offenders = abc)
        #Trap_Poisson_Params <- Trap_Po_Parameters(Data_parameters)
        #Source_Probabilities <- Multisource_probs(Data_parameters, Trap_Poisson_Params)
        misses <- subset(Trap_Data, Trap_Data[,3] ==0)
        hits <- subset(Trap_Data, Trap_Data[,3] >0)
        ratio_matrix[i, 1] <- sum(hits[,3])
        ratio_matrix[i, 2] <- length(misses[,3])
        if(length(misses[,3]) == 0)
        {ratio_matrix[i, 3] <- 0}
        else
        {ratio_matrix[i, 3] <- sum(hits[,3])/length(misses[,3])}
        #ratio_matrix[i, 4] <- Michael_geoReportHitscores(Data_parameters, s, Source_Probabilities$Source_Both)[,3] - Michael_geoReportHitscores(Data_parameters, s, Source_Probabilities$Source_Hits)[,3]
        #message("repeats = ",dQuote(i))
        i <- i + 1
      }
      else{}
      }
      #message("n_offenders = ",dQuote(index))
      index <- j + (k - 1)*length(n_offenders)
      message("index= ",dQuote(index))

      average_ratio_mat[index, 1] <- n_offenders[j]
      average_ratio_mat[index, 2] <- n_traps_x[k]^2
      average_ratio_mat[index, 3] <- mean(ratio_matrix[,3])
      j <- j + 1
    }
      #message("n_traps_x = ",dQuote(k))
      k <- k + 1
  }
}
save(average_ratio_mat, n_traps_x, n_offenders, file= "TrapToOffender.rdata")
which(average_ratio_mat[,3] ==0)
max(average_ratio_mat[,3])

load("TrapToOffender.rdata")
data_matrix <- matrix(average_ratio_mat[,3], ncol = length(n_traps_x), nrow = length(n_offenders))
traps <- n_traps_x*n_traps_x
hit_to_miss <- data_matrix
miss_to_hit <- 1/data_matrix
x11()
contour(n_offenders, traps, hit_to_miss, nlevels = 20, xlab = "Traps", ylab = "Number of Offenders", col= "green") #, theta = 260, zlim = c(0,100))

################################################################################
########################## FINDING SENSIBLE PARAMS #############################
################################################################################
par(mfrow=c(1,2))
n_sources <- 2
n_offenders <- 25
r_n_offenders <- rpois(1,n_offenders)
sigma <- runif(1, 1, 3)
tau <- runif(1, 1, 1)

  sim <- rDPM(r_n_offenders, alpha = 1, sigma = sigma, tau)
  point_data <- geoData(sim$longitude, sim$latitude)
  ##### only need this for the lon/lat min max value with combine misses and hits - that's all, the rest of the params don't matter
  master_params <- geoParams(data = point_data, sigma_mean = 1, sigma_squared_shape = 2,  samples= 50000, chains = 200, burnin = 2000, priorMean_longitude = mean(point_data$longitude), priorMean_latitude = mean(point_data$latitude), guardRail = 0.05)
  s <- geoDataSource(sim$source_lon, sim$source_lat)
  Trap_Data <- trap_assignment_data(simulation = sim, sim_params = master_params, sources = s, trap_range = c(10,20),"rDPM data")

  Standard_dev <- Trap_Data$Sd
  detection_TR <- Trap_Data$detection_TR
  Trap_Data <- Trap_Data$Trap_data

  ### DPM ###
  dpm_hits <- subset(Trap_Data, Trap_Data[,3]>0)
  dpm_hits <- dpm_hits[rep(1:nrow(dpm_hits), dpm_hits[,3]),]
  dpm_misses <- subset(Trap_Data, Trap_Data[,3]==0)
  trap_loc_data <- geoData(Trap_Data[,1], Trap_Data[,2])
  hit_data <- geoData(dpm_hits[,1], dpm_hits[,2])

  hit_params <- geoParams(data = hit_data, sigma_mean = Standard_dev, sigma_squared_shape = 2, samples= 50000, chains = 200, burnin = 2000, priorMean_longitude = mean(hit_data$longitude), priorMean_latitude = mean(hit_data$latitude), guardRail = 0.05)
  hit_params$output$longitude_minMax <- master_params$output$longitude_minMax
  hit_params$output$latitude_minMax <- master_params$output$latitude_minMax
  m <- geoMCMC(data=hit_data, params= hit_params)

  ### PA ###
  Trap_Data <- as.data.frame(Trap_Data)
  Data_parameters <- Extract_Params(sim, Trap_Data, PA_x_grid_cells = 10, PA_y_grid_cells = 10, Guard_Rail = 0.05, Trap_Radius = detection_TR, n_sources = n_sources, n_cores = 1, n_offenders = n_offenders, Time = 1, Sd_x = Standard_dev, Sd_y = Standard_dev)
  Trap_Poisson_Params <- Trap_Po_Parameters(Data_parameters)
  Source_Probabilities <- Multisource_probs(Data_parameters, Trap_Poisson_Params)
  Michael_geoReportHitscores(Data_parameters, s, Source_Probabilities$Source_Hits)[,3]
  Michael_geoReportHitscores(Data_parameters, s, Source_Probabilities$Source_Both)[,3]
  geoReportHitscores(hit_params, s, m$posteriorSurface)[,3]

  #x11()
  plot_sources(Data_parameters, Source_Probabilities)
  points(s$source_longitude, s$source_latitude, pch = 15, col = "blue")

Source_Probabilities$LSource_Both <- expandMatrix(Source_Probabilities$Source_Both, 500, 500)
Source_Probabilities$LSource_Hits <- expandMatrix(Source_Probabilities$Source_Hits, 500, 500)

### PLOT ###

x11()
geoPlotMap(data = hit_data, source = s, params = master_params, breakPercent = seq(0, 100, 10), mapType = "roadmap", contourCols =c("darkred", "red", "orange", "yellow"),
  					crimeCol = "darkgreen", crimeCex = 2, sourceCol = "red", sourceCex = 2, surface = m$geoProfile)

x11()
geoPlotMap(data = hit_data, source = s, params = master_params, breakPercent = seq(0, 100, 10), mapType = "roadmap", contourCols =c("darkred", "red", "orange", "yellow"),
					crimeCol = "darkgreen", crimeCex = 2, sourceCol = "red", sourceCex = 2, surface = rank(-Source_Probabilities$LSource_Both))

x11()
geoPlotMap(data = hit_data, params = master_params, source = s, breakPercent = seq(0, 10,  1), mapType = "roadmap", contourCols =c("darkred","red", "orange", "yellow"),
           crimeCol = "darkgreen", crimeCex = 2, sourceCol = "red", sourceCex = 2, surface = rank(-Source_Probabilities$LSource_Hits))

################################################################################

######  sim some data, and by eye find the approriate looking maps where the misses have a decent impact
######  list off each value for the parameters, then fiddle with each indivdually,
######  sigma, time, TR rho, Number of offenders, number of traps
###### then sim up another set of data and set the parameters etc to find ratio change t

################################################################################

################################################################################
###################### PLOTLY ########################
################################################################################

library(reshape2)
library(plotly)
time.seq <- seq(0,4, 0.1)
tr.seq <- seq(0,4, 0.1)

test_params <- function(time_sequence, trap.R.sequence)
{
	hs_matrix <- matrix(NA, nrow = length(time_sequence), ncol = length(trap.R.sequence))
for(ti in 1:length(time_sequence))
	{

		for(tr in 1:length(trap.R.sequence))
		{
		Data_parameters <- Extract_Params(sim, Trap_Data, x_grid_cells = 50, y_grid_cells = 50, Guard_Rail = 0.05, Trap_Radius = trap.R.sequence[tr], n_sources = 1, n_cores = 1, Time = time_sequence[ti])
		Trap_Poisson_Params <- Trap_Po_Parameters(Data_parameters)
		Source_Probabilities <- Multisource_probs(Data_parameters, Trap_Poisson_Params)
		hs_matrix[ti,tr] <- Source_Probabilities$Source_Miss[25,25]
		}
		print(ti)
	}
	return(list(hs_matrix = hs_matrix, time_sequence = time_sequence, trap.R.sequence=trap.R.sequence))
}
param_test <- test_params(time.seq, tr.seq)

grid_points <- expand.grid(param_test$time_sequence, param_test$trap.R.sequence)
hs_vector <- c(param_test$hs_matrix)
z_val_mat <- cbind(grid_points, hs_vector)
head(z_val_mat)

hs_z <- acast(z_val_mat, Var1~Var2, value.var = "hs_vector")
plot_ly(z = hs_z, type = "surface")


################################################################################
################################################################################
################################################################################
################################################################################

n_offenders <- 25
n_sources <- 1

sim <- rDPM(n_offenders, alpha = 0, sigma = 1)
point_data <- geoData(sim$longitude, sim$latitude)
master_params <- geoParams(data = point_data, sigma_mean = 1, sigma_squared_shape = 2,  samples= 50000, chains = 200, burnin = 2000, priorMean_longitude = mean(point_data$longitude), priorMean_latitude = mean(point_data$latitude), guardRail = 0.05)
par(mfrow=c(1,2))
s <- geoDataSource(sim$source_lon, sim$source_lat)
Trap_Data <- trap_assignment_data(sim, master_params, sources =s, c(5,20),"Title")
detection_TR <- Trap_Data$detection_TR
Trap_Data <- Trap_Data$Trap_data

### DPM ###

dpm_hits <- subset(Trap_Data, Trap_Data[,3]>0)
dpm_hits <- dpm_hits[rep(1:nrow(dpm_hits), dpm_hits[,3]),]
dpm_misses <- subset(Trap_Data, Trap_Data[,3]==0)
trap_loc_data <- geoData(Trap_Data[,1], Trap_Data[,2])
s <- geoDataSource(sim$source_lon, sim$source_lat)
hit_data <- geoData(dpm_hits[,1], dpm_hits[,2])
points <- cbind(hit_data$longitude, hit_data$latitude)

################################################
### FIX params ####
################################################

fixed_sigma <- c()
fixed_sigma_HS <- c()
fixed_TR <- c()
fixed_TR_HS <- c()
fixed_NO <- c()
fixed_NO_HS <- c()

pw_dist <- pairwise_distance(points)
Sd <- 0.5*mean(pw_dist$distance, na.rm = TRUE)

hit_params <- geoParams(data = hit_data, sigma_mean = Trap_Data$Sd, sigma_squared_shape = 2, samples= 50000, chains = 200, burnin = 2000, priorMean_longitude = mean(hit_data$longitude), priorMean_latitude = mean(hit_data$latitude), guardRail = 0.05)
hit_params$output$longitude_minMax <- master_params$output$longitude_minMax
hit_params$output$latitude_minMax <- master_params$output$latitude_minMax
m <- geoMCMC(data=hit_data, params= hit_params)
Trap_Data <- as.data.frame(Trap_Data)

geoReportHitscores(master_params, s, m$posteriorSurface)

Sigma <- seq(0.5, 10, 0.05)
TRS <- seq(0.5, 10, 0.05)
NO_off <- seq(1,191, 1)

for(i in 1:length(NO_off))

Data_parameters <- Extract_Params(sim, Trap_Data, PA_x_grid_cells = 50, PA_y_grid_cells = 50, Guard_Rail = 0.05, Trap_Radius = detection_TR, n_sources = n_sources, n_cores = 1, n_offenders = NO_off[i], Time = 1)
Trap_Poisson_Params <- Trap_Po_Parameters(Data_parameters)
Source_Probabilities <- Multisource_probs(Data_parameters, Trap_Poisson_Params)

#fixed_sigma[i] <- (Data_parameters$Trap_Radius^2)*Data_parameters$n_offenders/(Sigma[i]^2)
#fixed_sigma_HS[i] <- Michael_geoReportHitscores(Data_parameters, s, Source_Probabilities$Source_Both)[,3]

#fixed_TR[i] <- (TRS[i]^2)*(Data_parameters$n_offenders)/(Data_parameters$Sd_x*Data_parameters$Sd_y)
#fixed_TR_HS[i] <- Michael_geoReportHitscores(Data_parameters, s, Source_Probabilities$Source_Both)[,3]

#fixed_NO[i] <-  (Data_parameters$Trap_Radius^2)*(i)/(Data_parameters$Sd_x*Data_parameters$Sd_y)
#fixed_NO_HS[i] <- Michael_geoReportHitscores(Data_parameters, s, Source_Probabilities$Source_Both)[,3]


par(mfrow=c(2,2))

Trap_Data <- trap_assignment_data(sim, master_params, 5,5,"Title")

plot(fixed_sigma, fixed_sigma_HS, type = 'l', main = "Altering Sigma", ylab = "Hitscore", xlab = expression(rho^2*lambda/sigma^2))
abline(0.001216, 0)

plot(fixed_TR, fixed_TR_HS, type = 'l',  main = "Altering Trap Radius", ylab = "Hitscore", xlab = expression(rho^2*lambda/sigma^2))
abline(0.001216, 0)

plot(fixed_NO, fixed_NO_HS, type = 'l',  main = "Altering Number  of offenders", ylab = "Hitscore", xlab = expression(rho^2*lambda/sigma^2))
abline(0.001216, 0)


### save ###

setwd("/home/mstevens/Desktop/Main\ work/Presence\ Absence/Presence_Absence-master/Finals")

save(Sigma, TRS, NO_off, fixed_sigma, fixed_TR, fixed_TR_HS, fixed_sigma_HS, fixed_NO, fixed_NO_HS, file = "params.rdata")
