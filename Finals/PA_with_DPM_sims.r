rm(list = ls()) #remove all objects
library(flux)

library(gRbase) ## For the "combnPrim" function, efficiently listing combinations
library(parallel)
# _______________________________________________________
# SIM DATA AND SET PARAMS
# _______________________________________________________
library(RgeoProfile)
# generate simulated data
#sim <- rDPM(50, priorMean_longitude = -0.04217481, priorMean_latitude = 51.5235505, alpha = 1, sigma = 1, tau = 12)
#sim
###########################################################################################################################################

# allocate hits and misses to crimes
#hits <- sample(0:1,length(sim$longitude),replace=TRUE)

#trap_data <-  cbind(sim$longitude,sim$latitude, hits)
#colnames(trap_data) <- c("x", "y", "hits")
#head(trap_data)

setwd("/home/mstevens/Desktop/Presence Absence/Presence_Absence-master/Finals/data")
trap_data <- read.table("FootballToyExample.txt", header = FALSE)

dummy_source <- read.table("Joe_Sources.txt", header = FALSE)


# _______________________________________________________
# RUN AS DPM TO EXTRACT PARAMS
# _______________________________________________________
# convert
d <- geoData(trap_data$V1,trap_data$V2)
s <- geoDataSource(dummy_source$V1,dummy_source$V2)

# params
params <- geoParams(data = d, sigma_mean = 1, sigma_squared_shape = 2, samples= 100000, chains = 200, burnin = 10000, priorMean_longitude = mean(d$longitude), priorMean_latitude = mean(d$latitude), guardRail = 0.05)
params$output$longitude_cells <-300
params$output$latitude_cells <- 300

params$output$longitude_minMax
all_params$Long_Max_Bound

m <- geoMCMC(data=d,params= params)

##_____________FUNCTIONS ______________##
###########################################################################################################################################
#
# The "pairwise_distance" function returns half the average of the minimum distances between traps. This value takes the place of the
# x and y standard deviations in our model.
#
###########################################################################################################################################

pairwise_distance <- function(points){
	distance <- matrix(NA, ncol = length(points[,1]), nrow = length(points[,1]))

	for(i in 1:length(points[,1]))
	{
	for(j in 1:i)
		{
		x1 <- points[i, 1]
		y1 <- points[i, 2]
		x2 <- points[j, 1]
		y2 <- points[j, 2]
		dist <- sqrt(abs(x1 - x2)^2 + abs(y1 - y2)^2)
		distance[i,j] <- dist
		}
	}
	diag(distance)  <- NA
	distance[distance==0] <- NA
	distance_min <- apply(distance, 1, min, na.rm = TRUE)
	distance_min[distance_min =="Inf"] <- NA
	TR <- 0.5*mean(distance_min, na.rm = TRUE)
	SD_TR <- list(distance = distance, TR = TR)
	return(SD_TR)
}

###########################################################################################################################################
#
# The "Extract_Params" function returns a list of parameters for an "n_sources" source model. Should "n_sources" be set to anything greater
# than one, "Extract_Params" will produce an additional parameter "AP_allocation," a list ilustrating the number of ways of choosing
# "n_sources" sources from the total number of grid cells.
#
###########################################################################################################################################

Extract_Params <- function(Trap_Data, x_grid_cells = 10, y_grid_cells = 10, Time = 1, Guard_Rail = 0.05, Trap_Radius = 0.6, n_sources = 1, n_cores = 1)
			{
			colnames(Trap_Data) <- c("Longitude", "Latitude", "Hits")
			Trap_Data <- data.frame(Trap_Data)
			long_diff <- (max(Trap_Data$Longitude) - min(Trap_Data$Longitude))
			lat_diff <- (max(Trap_Data$Latitude) - min(Trap_Data$Latitude))

			MnMx_Long <- c(min(Trap_Data$Longitude), max(Trap_Data$Longitude))
			MnMx_Lat <- c(min(Trap_Data$Latitude), max(Trap_Data$Latitude))

			Long_Max_Bound <- MnMx_Long[2] + (Guard_Rail*long_diff)
			Long_Min_Bound <- MnMx_Long[1] - (Guard_Rail*long_diff)
			Lat_Max_Bound <-  MnMx_Lat[2] + (Guard_Rail*lat_diff)
			Lat_Min_Bound <- MnMx_Lat[1] - (Guard_Rail*lat_diff)
			Anchor_Points_Long <- seq(min(Long_Min_Bound,Long_Max_Bound), max(Long_Min_Bound,Long_Max_Bound), abs(Long_Max_Bound-Long_Min_Bound)/(x_grid_cells - 1))
			Anchor_Points_Lat <- seq(min(Lat_Min_Bound, Lat_Max_Bound), max(Lat_Min_Bound, Lat_Max_Bound), abs(Lat_Max_Bound-Lat_Min_Bound)/(y_grid_cells - 1))
			if(n_sources > 1)
      {
        AP_allocation <- t(combnPrim(x_grid_cells*y_grid_cells, n_sources))
      } else {
        AP_allocation <- "NA"
        }
      Hits_Only <- subset(Trap_Data, Hits != 0)
			Miss_Only <- subset(Trap_Data, Hits == 0)
			pairwise <- pairwise_distance(Trap_Data)
			Sd_x <- 0.2*mean(pairwise$distance, na.rm = TRUE)
			Sd_y <- 0.2*mean(pairwise$distance, na.rm = TRUE)

      params <- list(n_sources = n_sources, x_grid_cells = x_grid_cells, y_grid_cells = y_grid_cells, Sd_x = Sd_x, Sd_y = Sd_y,
				    Trap_Radius=Trap_Radius, Time=Time, Anchor_Points_Long=Anchor_Points_Long, Anchor_Points_Lat=Anchor_Points_Lat,
				    AP_allocation=AP_allocation, Hits_Only=Hits_Only, Miss_Only=Miss_Only, Long_Max_Bound = Long_Max_Bound, n_cores = n_cores)
			return(params)
			}

###########################################################################################################################################
#
# The "Poisson_Parameter" function returns the average number of hits we expect to see within a trap located at "(x, y)." This value
# is an approximation of the integral of a bivariate normal distribution bounded by the trap radius. We state expected number of hits
# is a value observed over some time interval, thus the multiplication by "t."
#
###########################################################################################################################################

Poisson_Parameter <- function(x, y, mu_x, mu_y, trap_radius, t, Sd_x, Sd_y)
                     	{
				co_efficient <- t*4*(trap_radius^2)*(2*pi*Sd_x*Sd_y)^(-1)
				param <- co_efficient*exp( -0.5*( (x - mu_x)^2/(Sd_x*Sd_x) + (y - mu_y)^2/(Sd_y*Sd_y) ) )
				return(param)
                        }

########################################################################################################################################
#
# The "Trap_Po_Parameters" function returns the poisson parameter for each trap in the form of an array. Where the individual entries of
# matrix i represent the Poisson parameters of event i.
#
#######################################################################################################################################

Trap_Po_Parameters <- function(Params)
			    {
					Po_Array_Hits <- array(NA, c(length(Params$Anchor_Points_Long), length(Params$Anchor_Points_Lat), length(Params$Hits_Only$Longitude) ) )
					Po_Array_Miss <- array(NA, c(length(Params$Anchor_Points_Long), length(Params$Anchor_Points_Lat), length(Params$Miss_Only$Longitude) ) )

		for(i in 1:length(Params$Anchor_Points_Long))
			{
			for(j in 1:length(Params$Anchor_Points_Lat))
				{
				Po_hits <- mapply(Poisson_Parameter, Params$Hits_Only$Longitude, Params$Hits_Only$Latitude, mu_x = Params$Anchor_Points_Long[i],
					      mu_y = Params$Anchor_Points_Lat[j], trap_radius = Params$Trap_Radius, t = Params$Time,
						Sd_x = Params$Sd_x, Sd_y = Params$Sd_y)
				Po_miss <- mapply(Poisson_Parameter, Params$Miss_Only$Longitude, Params$Miss_Only$Latitude, mu_x = Params$Anchor_Points_Long[i],
					      mu_y = Params$Anchor_Points_Lat[j], trap_radius = Params$Trap_Radius, t = Params$Time,
						Sd_x = Params$Sd_x, Sd_y = Params$Sd_y)
				Po_Array_Hits[i, j,] <-  Po_hits
				Po_Array_Miss[i, j,] <-  Po_miss
				}
			}
			Po_Param <- list(Po_Array_Hits = Po_Array_Hits, Po_Array_Miss = Po_Array_Miss)
			return(Po_Param)
}

###########################################################################################################################################
#
# The "Multisource_probs" function returns the porbability for source locations. Given the user specifies the number of sources expected.
#
###########################################################################################################################################

Multisource_probs <- function(Data_Params, Po_Params)
		{
    ####################################################
    # SINGLE -  SOURCE
		####################################################
		if(Data_Params$n_sources == 1)
		{
		SS_Array_Hits <- array(NA, c(length(Po_Params$Po_Array_Hits[,1,1]), length(Po_Params$Po_Array_Hits[1,,1]), length(Po_Params$Po_Array_Hits[1,1,]) ) )
		SS_Array_Miss <- array(NA, c(length(Po_Params$Po_Array_Miss[,1,1]), length(Po_Params$Po_Array_Miss[1,,1]), length(Po_Params$Po_Array_Miss[1,1,]) ) )

			for(i in 1:length(Po_Params$Po_Array_Hits[,1,1]))
				{
				for(j in 1:length(Po_Params$Po_Array_Hits[1,,1]))
					{
					Hits_Prob <- mcmapply(dpois, Data_Params$Hits_Only$Hits, Po_Params$Po_Array_Hits[i,j,], mc.cores = Data_Params$n_cores)
					Miss_Prob <- mcmapply(dpois, Data_Params$Miss_Only$Hits, Po_Params$Po_Array_Miss[i,j,], mc.cores = Data_Params$n_cores)
					SS_Array_Hits[i,j,] <- Hits_Prob
					SS_Array_Miss[i,j,] <- Miss_Prob
					print(i)
					}
				}
				SS_Hits <- exp(apply(log(SS_Array_Hits), c(1,2), sum))
				SS_Miss <- exp(apply(log(SS_Array_Miss), c(1,2), sum))
				SS_Both <- SS_Hits * SS_Miss
				Source_Hits <- SS_Hits/sum(SS_Hits)
				Source_Miss <- SS_Miss/sum(SS_Miss)
				Source_Both <- SS_Both/sum(SS_Both)
				Source_Prob <- list(Source_Hits = Source_Hits, Source_Miss = Source_Miss, Source_Both = Source_Both)

				return(Source_Prob)

			}	 else{
					Sum_hit_po <- matrix(NA, ncol = length(Data_Params$Hits_Only$Longitude), nrow =length(Data_Params$AP_allocation[,1]))
  				Sum_miss_po <- matrix(NA, ncol = length(Data_Params$Miss_Only$Longitude), nrow =length(Data_Params$AP_allocation[,1]))

  				################################################################################################
					print("Calculate Poisson Parameters by summing hazard surfaces")
					Sys.sleep(5)
					print("Poisson Parameters for the Hits")
					Sys.sleep(3)
					for(j in 1:length(Data_Params$Hits_Only[,1]))
						{
							matrix_h <- Po_Params$Po_Array_Hits[,,j]
							Sum_hit_po[, j] <- parRapply(cluster, Data_Params$AP_allocation, FUN = function(x) sum(matrix_h[x]))
							print(j)
						}
						print("Poisson Parameters for the Miss")
						Sys.sleep(3)
						for(j in 1:length(Data_Params$Miss_Only[,1]))
						{
							matrix_m <- Po_Params$Po_Array_Miss[,,j]
							Sum_miss_po[, j] <- parRapply(cluster, Data_Params$AP_allocation, FUN = function(x) sum(matrix_m[x]))
							print(j)
						}

					S_hit_prob <- matrix(NA, ncol = length(Data_Params$Hits_Only$Longitude), nrow =length(Data_Params$AP_allocation[,1]))
					S_miss_prob <- matrix(NA, ncol = length(Data_Params$Miss_Only$Longitude), nrow =length(Data_Params$AP_allocation[,1]))
					print("Calculate probabilities for observing data given our fixed sources")
					Sys.sleep(5)
					print("For the Hits")
					Sys.sleep(3)
					for(i in 1:length(Data_Params$Hits_Only[,1]))
					{
						S_hit_prob[,i] <- mcmapply(dpois, Data_Params$Hits_Only$Hits[i], Sum_hit_po[,i], mc.cores = Data_Params$n_cores)
						print(i)
					}
					S_hit_prob <- log(S_hit_prob)
					S_hit_prob <- exp(apply(S_hit_prob, FUN = sum, MARGIN = 1))
					print("For the misses")
					Sys.sleep(3)
					for(i in 1:length(Data_Params$Miss_Only[,1]))
					{
						S_miss_prob[,i] <- mcmapply(dpois, Data_Params$Miss_Only$Hits[i], Sum_miss_po[,i], mc.cores = Data_Params$n_cores)
						print(i)
					}

					S_miss_prob <- log(S_miss_prob)
					S_miss_prob <- exp(apply(S_miss_prob, FUN = sum, MARGIN = 1))

					S_both_prob <- S_hit_prob*S_miss_prob
  				S_probs <- cbind(Data_Params$AP_allocation, S_hit_prob, S_miss_prob, S_both_prob)
				################################################################################################
				final_hits <- c()
				final_miss <- c()
				final_both <- c()
				print("Sum over probabilities to obtain likelihood for individual sources")
				Sys.sleep(5)
				for(i in 1:((Data_Params$x_grid_cells)*(Data_Params$y_grid_cells)))
				{
					print(i)
					a <- parRapply(cluster, S_probs[,1:Data_Params$n_sources], FUN=function (x) any(x == i))
					final_hits[i] <- sum(S_probs[a, Data_Params$n_sources + 1])
					final_miss[i] <- sum(S_probs[a, Data_Params$n_sources + 2])
					final_both[i] <- sum(S_probs[a, Data_Params$n_sources + 3])
				}
				Source_Hits <- matrix(final_hits, ncol = length(Data_Params$Anchor_Points_Lat), nrow = length(Data_Params$Anchor_Points_Long))
				Source_Miss <- matrix(final_miss, ncol = length(Data_Params$Anchor_Points_Lat), nrow = length(Data_Params$Anchor_Points_Long))
			  Source_Both <- matrix(final_both, ncol = length(Data_Params$Anchor_Points_Lat), nrow = length(Data_Params$Anchor_Points_Long))

				Source_Hits <-  Source_Hits/sum(Source_Hits)
				Source_Miss <-  Source_Miss/sum(Source_Miss)
				Source_Both <-  Source_Both/sum(Source_Both)
				Source_Prob <- list( Source_Hits = Source_Hits, Source_Miss = Source_Miss, Source_Both = Source_Both)
				return(Source_Prob)
				}
			}

###########################################################################################################################################
# function to resize this sub-matrix to the original resolution
# mat - matrix to be resized
# output_long/lat - new number of long/let cells
###########################################################################################################################################

expandMatrix <- function(mat,output_long,output_lat)
  {
	# define function expanding vector
	expandVector <- function(input_vec,output_length)
		{
			my_vec <- input_vec
			desired_length <- output_length
			new_vec <- rep(NA, desired_length)

			vec_ID <- seq(1,length(my_vec),length=desired_length)

			for(i in 1:length(new_vec))
				{
					ifelse(vec_ID[i] %% 1 == 0,
		new_vec[i] <- my_vec[floor(vec_ID[i])],
		new_vec[i] <- mean((1-vec_ID[i] %% 1) * my_vec[floor(vec_ID[i])] + (vec_ID[i] %% 1) * my_vec[ceiling(vec_ID[i])])
	)
				}
  return(new_vec)
		}
  mat1 <- apply(mat,2, function(x) expandVector(x, output_long))
  mat2 <- apply(mat1,1, function(x) expandVector(x, output_lat))
  return(t(mat2))
  }


###########################################################################################################################################
#
# The "plot1source" function plots a contour map of the single source probabilities
#
###########################################################################################################################################

plot_sources <- function(Data_Params, Probs)
		   {
				 if(Data_Params$n_sources == 1)
				 {
				 #x11()
				 contour(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat, Probs$SS_Hits, col = "darkgreen", nlevels = 5)
				 contour(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat, Probs$SS_Miss,col = "red",add=TRUE, nlevels = 5)
				 contour(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat, Probs$SS_Both,col = "blue", add=TRUE, nlevels = 5)
				 points(Data_Params$Hits_Only$Longitude, Data_Params$Hits_Only$Latitude , pch = 16, col = "green")
				 points(Data_Params$Miss_Only$Longitude, Data_Params$Miss_Only$Latitude , pch = 16, col = "red")

				 #x11()
				 #par(mfrow=c(1,3))
				 #persp(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat,  SS_Prob$SS_Hits, col = "green")
				 #persp(Data_Params$Anchor_Points_Long,Data_Params$Anchor_Points_Lat,  SS_Prob$SS_Miss, col = "red")
				 #persp(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat,  SS_Prob$SS_Both, col = "blue")

				 #par(mfrow=c(2,2))
			 } else {
				 contour(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat, Probs$Source_Hits, col = "darkgreen", nlevels = 5)
				 contour(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat, Probs$Source_Miss,col = "red",add=TRUE, nlevels = 3)
				 contour(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat, Probs$Source_Both,col = "blue",add=TRUE, nlevels = 5)
				 points(Data_Params$Hits_Only$Longitude, Data_Params$Hits_Only$Latitude , pch = 16, col = "green")
				 points(Data_Params$Miss_Only$Longitude, Data_Params$Miss_Only$Latitude , pch = 16, col = "red")
				#symbols(Data_Params$Hits_Only$Longitude, Data_Params$Hits_Only$Latitude, circles = rep(Data_Params$Trap_Radius, length(Data_Params$Hits_Only$Longitude)), add = T, inches = F)

				#contour(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat, Two_Source$TwoD_Hits, col = "darkgreen", nlevels = 5)
				#	points(Data_Params$Hits_Only$Longitude, Data_Params$Hits_Only$Latitude , pch = 16, col = "green")

				#contour(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat, Two_Source$TwoD_Miss,col = "red", nlevels = 5)
				#points(Data_Params$Miss_Only$Longitude, Data_Params$Miss_Only$Latitude , pch = 16, col = "red")

				#contour(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat, Two_Source$TwoD_Both,col = "blue", nlevels = 5)
				#points(Data_Params$Hits_Only$Longitude, Data_Params$Hits_Only$Latitude , pch = 16, col = "green")
				#points(Data_Params$Miss_Only$Longitude, Data_Params$Miss_Only$Latitude , pch = 16, col = "red")
				#points(DPM_SIM$source_lon, DPM_SIM$source_lat, pch = 16, col = "blue")
				#
				#x11()
				#par(mfrow=c(1,3))
				#persp(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat, Two_Source$TwoD_Hits,  col = "green")
				#persp(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat, Two_Source$TwoD_Miss,  col = "red")
				#persp(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat, Two_Source$TwoD_Both,  col = "blue")
			}
			}

###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################

# _______________________________________________________
# RUN PRESENCE/ABSENCE WITH THESE PARAMS
# _______________________________________________________
trap_data <- as.data.frame(trap_data)
start.time <- Sys.time()
## Extract ALL Parameters from your data
time.taken <- end.time - start.time
time.taken

all_params <- Extract_Params(My_trap_data, x_grid_cells = 20, y_grid_cells = 20, Guard_Rail = 0.5, Trap_Radius = 0.28, n_sources = 3, n_cores = 3)

## Compute Poisson Parameters
po_params <- Trap_Po_Parameters(all_params)

## Compute probability matrices in parallel
cluster <- makeCluster(3)
clusterExport(cluster, c("Data_parameters", "Trap_Poisson_Params"))
source_probs <- Multisource_probs(Data_parameters, Trap_Poisson_Params)
#stopCluster(cluster)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

###########################################################################################################################################

hits <- geoData(all_params$Hits_Only$Longitude, all_params$Hits_Only$Latitude)

misses <- geoDataSource(all_params$Miss_Only$Longitude,all_params$Miss_Only$Latitude)


######

# the map alone

a <- expandMatrix(source_probs$Source_Both, 300, 300)

x11()
geoPlotMap(data = d, source = s, params = params, breakPercent = seq(0, 20, 2), mapType = "roadmap", contourCols =c("red", "orange", "yellow", "white"),
           crimeCol = "black", crimeCex = 2, sourceCol = "red", sourceCex = 2, surface = m$geoProfile)

######

# both
x11()

geoPlotMap(data = hits, source = misses, params = params, breakPercent = seq(0, 5, 0.5), mapType = "roadmap", contourCols =c("red", "orange", "yellow", "white"),
           crimeCol = "darkgreen", crimeCex = 5, sourceCol = "red", sourceCex = 5, surface = rank(-a))

# hits
x11()
geoPlotMap(data = hits, params = params, breakPercent = seq(0, 100, 10), mapType = "roadmap", contourCols =c("red", "orange", "yellow", "white"),
           crimeCol = "darkgreen", crimeCex = 5, sourceCol = "red", sourceCex = 5, surface = rank(-source_probs$Source_Hits))


# misses
x11()
geoPlotMap(data = hits, source = misses, params = params, breakPercent = seq(0, 100, 10), mapType = "roadmap", contourCols =c("red", "orange", "yellow", "white"),
           crimeCol = "darkgreen", crimeCex = 5, sourceCol = "red", sourceCex = 5, surface = rank(-source_probs$Source_Miss))

#x11()
#perspGP(surface=t(probs$TwoD_miss),aggregate_size=5,surface_type="prob",phiGP=70,thetaGP=-10)

#x11()
#perspGP(surface=t(probs$TwoD_Both),aggregate_size=3,surface_type="prob")
#############################################################################################################################################
