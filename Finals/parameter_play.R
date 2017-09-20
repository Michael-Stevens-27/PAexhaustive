################################################################################
############################## HIT MISS RATIO ##################################
################################################################################
n_offenders <- seq(10,100, 1)
n_traps_x <-  seq(3, 15, 1)
n_traps_y <- seq(3, 15, 1)
repeats <- 50
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
        sim <- rDPM(abc, priorMean_longitude = -0.04217481, priorMean_latitude = 51.5235505, alpha = 0, sigma = 1)
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
        ratio_matrix[i, 3] <- sum(hits[,3])/length(misses[,3])
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

load("TrapToOffender.rdata")
data_matrix <- matrix(average_ratio_mat[,3], ncol = length(n_traps_x), nrow = length(n_offenders))
traps <- n_traps_x*n_traps_y
hit_to_miss <- data_matrix
miss_to_hit <- 1/data_matrix
x11()
contour(n_offenders, traps, miss_to_hit, nlevels = 50, xlab = "Traps", ylab = "Number of Offenders") #, theta = 260, zlim = c(0,100))

################################################################################
########################## FINDING SENSIBLE PARAMS #############################
################################################################################

n_offenders <- 49
n_soucres <- 2
sim <- rDPM(n_offenders, alpha = 1, sigma = 1)

point_data <- geoData(sim$longitude, sim$latitude)
master_params <- geoParams(data = point_data, sigma_mean = 1, sigma_squared_shape = 2, samples= 50000, chains = 200, burnin = 2000, priorMean_longitude = mean(point_data$longitude), priorMean_latitude = mean(point_data$latitude), guardRail = 0.05)
Trap_Data <- trap_assignment_data(sim, master_params, 7, 7, 1)
detection_TR <- Trap_Data$detection_TR
Trap_Data <- Trap_Data$Trap_data

### DPM ###
dpm_hits <- subset(Trap_Data, Trap_Data[,3]>0)
dpm_hits <- dpm_hits[rep(1:nrow(dpm_hits), dpm_hits[,3]),]
dpm_misses <- subset(Trap_Data, Trap_Data[,3]==0)
trap_loc_data <- geoData(Trap_Data[,1], Trap_Data[,2])
s <- geoDataSource(sim$source_lon, sim$source_lat)
hit_data <- geoData(dpm_hits[,1], dpm_hits[,2])
hit_params <- geoParams(data = hit_data, sigma_mean = 1, sigma_squared_shape = 2, samples= 50000, chains = 200, burnin = 2000, priorMean_longitude = mean(hit_data$longitude), priorMean_latitude = mean(hit_data$latitude), guardRail = 0.05)
hit_params$output$longitude_minMax <- master_params$output$longitude_minMax
hit_params$output$latitude_minMax <- master_params$output$latitude_minMax
#m <- geoMCMC(data=hit_data, params= hit_params)

### PA ###
Trap_Data <- as.data.frame(Trap_Data)
Data_parameters <- Extract_Params(sim, Trap_Data, PA_x_grid_cells = 30, PA_y_grid_cells = 30, Guard_Rail = 0.05, Trap_Radius = detection_TR, n_sources = 2, n_cores = 1, n_offenders = n_offenders, Sd_x = 1, Sd_y = 1, Time = 5)
Trap_Poisson_Params <- Trap_Po_Parameters(Data_parameters)
Source_Probabilities <- Multisource_probs(Data_parameters, Trap_Poisson_Params)

#x11(display = "Original")
#geoPlotMap(data = misc_sim$hit_data, source = misc_sim$s, params = misc_sim$hit_params, breakPercent = seq(0, 10, 1), mapType = "roadmap", contourCols =c("darkred", "red", "orange", "yellow"),#
#          crimeCol = "darkgreen", crimeCex = 2, sourceCol = "blue", sourceCex = 2, surface = misc_sim$m$geoProfile)

Source_Probabilities$Source_Both <- expandMatrix(Source_Probabilities$Source_Both, 500, 500)
Source_Probabilities$Source_Hits <- expandMatrix(Source_Probabilities$Source_Hits, 500, 500)

x11()
geoPlotMap(data = hit_data, source = s, params = master_params, breakPercent = seq(0, 100, 10), mapType = "roadmap", contourCols =c("darkred", "red", "orange", "yellow"),
						crimeCol = "darkgreen", crimeCex = 2, sourceCol = "red", sourceCex = 2, surface = rank(-Source_Probabilities$Source_Both))

x11()
geoPlotMap(data = hit_data, params = master_params, source = s, breakPercent = seq(0, 100, 10), mapType = "roadmap", contourCols =c("darkred","red", "orange", "yellow"),
           crimeCol = "darkgreen", crimeCex = 2, sourceCol = "red", sourceCex = 2, surface = rank(-Source_Probabilities$Source_Hits))

################################################################################

######  sim some data, and by eye find the approriate looking maps where the misses have a decent impact
######  list off each value for the parameters, then fiddle with each indivdually,
######  sigma, time, TR rho, Number of offenders, number of traps
###### then sim up another set of data and set the parameters etc to find ratio change t
