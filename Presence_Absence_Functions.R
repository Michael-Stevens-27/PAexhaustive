rm(list = ls()) #remove all objects

#.libPaths(.libPaths()[2:4]) for R to access packages from a different location - USE ON CONATINER ONLY

library(flux)
library(gRbase)      ## For the "combnPrim" function, efficiently listing combinations
                     ## Initial install, run commands:
								     ## source("http://bioconductor.org/biocLite.R")
								     ## biocLite("gRbase")
library(parallel)    ## Package for paralellisation
library(RgeoProfile) ## Package to run the dirichlet process mixture model

##_____________FUNCTIONS ______________##
###### CLUSTER HITSCORE #################
geoReportHitscores <- function(params,source_data,surface)
  {
    sources <- cbind(source_data$source_longitude,source_data$source_latitude)
    ordermat = matrix(0, params$output$latitude_cells, params$output$longitude_cells)
    profile_order = order(surface)
    for (i in 1:(params$output$latitude_cells*params$output$longitude_cells))
    {
      ordermat[profile_order[i]] = i
    }
    hitscoremat <<- 1 - ordermat/(params$output$latitude_cells*params$output$longitude_cells)
    hitscoremat2 <- hitscoremat[nrow(hitscoremat):1, ]
    xvec = seq(params$output$longitude_minMax[1], params$output$longitude_minMax[2],
    length = params$output$longitude_cells)
    yvec = seq(params$output$latitude_minMax[1], params$output$latitude_minMax[2],
    length = params$output$latitude_cells)
    xdiff = abs(outer(rep(1, nrow(sources)), xvec) - outer(sources[,
    1], rep(1, params$output$longitude_cells)))
    ydiff = abs(outer(rep(1, nrow(sources)), yvec) - outer(sources[,
    2], rep(1, params$output$latitude_cells)))
    msourcex = mapply(which.min, x = split(xdiff, row(xdiff)))
    msourcey = params$output$longitude_cells - (mapply(which.min,
    x = split(ydiff, row(ydiff)))) + 1
    if(nrow(sources) > 1) {
    hitscores = diag(hitscoremat2[msourcey, msourcex])
    }
    else {
    hitscores = hitscoremat2[msourcey, msourcex]
    }
    hit_output <<- cbind(sources, hitscores)
    colnames(hit_output) <- c("lon", "lat", "hs")
    return(hit_output)
  }

################################################################################
# The "Poisson_Parameter" function returns the average number of hits we expect to see within a sentinel site located at "(x, y)." This value
# is an approximation of the integral of a bivariate normal distribution bounded by the site radius. We state expected number of hits
# is a value observed over some time interval, thus the multiplication by "t." We also multiply by the value "hazard_param," the expected number of
# sentinels we expect there to be within the entire search area.
# example: Poisson_Parameter(x = 0.25, y = 0.25, mu_x = 0, mu_y = 0, trap_radius = 15, t = 60, Sd_x = 7, Sd_y = 7, hazard_param = 10)
################################################################################

Poisson_Parameter <- function(x, y, mu_x, mu_y, trap_radius, t, Sd_x, Sd_y, hazard_param = 1)
                     	{
				co_efficient <- hazard_param*t*(pi*trap_radius^2)*(2*pi*Sd_x*Sd_y)^(-1)
				param <- co_efficient*exp( -0.5*( (latlon_to_bearing(y, x, y, mu_x)$gc_dist)^2/(Sd_x*Sd_x) + (latlon_to_bearing(y, mu_x, mu_y, mu_x)$gc_dist)^2/(Sd_y*Sd_y) ) )
				return(param)
                     }

################################################################################
# Function to resize this sub-matrix to the original resolution, mat-matrix to be resized, output_long/lat - new number of long/let cells
# Example: my_mat <- matrix(1:9, 3, 3)
# expandMatrix(my_mat, 20, 20)
################################################################################

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

################################################################################
# Function to calculate euclidean distance between two points specified in latitude and longitude
################################################################################

latlon_to_bearing <- function(origin_lat, origin_lon, dest_lat, dest_lon)

	{
	# convert input arguments to radians
	origin_lat <- origin_lat*2*pi/360
	dest_lat <- dest_lat*2*pi/360
	origin_lon <- origin_lon*2*pi/360
	dest_lon <- dest_lon*2*pi/360

	delta_lon <- dest_lon-origin_lon

	# calculate bearing and great circle distance
	bearing <- atan2(sin(delta_lon)*cos(dest_lat), cos(origin_lat)*sin(dest_lat)-sin(origin_lat)*cos(dest_lat)*cos(delta_lon))
	gc_angle <- acos(sin(origin_lat)*sin(dest_lat) + cos(origin_lat)*cos(dest_lat)*cos(delta_lon))

	# convert bearing from radians to degrees measured clockwise from due north, and convert gc_angle to great circle distance via radius of earth (km)
	bearing <- bearing*360/(2*pi)
	bearing <- (bearing+360)%%360
	earthRad <- 6371
	gc_dist <- earthRad*gc_angle

	return(list(bearing=bearing, gc_dist=gc_dist))
}

################################################################################
# The "pairwise_distance" function returns half the average of the minimum distances between sentinel sites.
# We use half the mean nearest neighbour distance (NND) between sites for site radius and half the mean NND between hits for the standard deviation
################################################################################

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
		dist <- latlon_to_bearing(y1, x1, y2, x2)$gc_dist
		distance[i,j] <- dist
		}
	}
	diag(distance)  <- NA
	distance[distance==0] <- NA
	distance_min <- apply(distance, 1, min, na.rm = TRUE)
	distance_min[distance_min =="Inf"] <- NA
	SD_TR <- list(distance = distance, distance_min = distance_min)
	return(SD_TR)
}

################################################################################
# The "Extract_Params" function returns a list of parameters for an "n_sources" source model. Should "n_sources" be set to anything greater
# than one, "Extract_Params" will produce an additional parameter "AP_allocation," a list illustrating the number of ways of choosing
# "n_sources" sources from the total number of grid cells, the most computationally expensive action.
# example: Extract_Params(cbind(1:10,1:10, sample(0:10, 10, replace=T)), PA_x_grid_cells = 10, PA_y_grid_cells = 10, Time = 1, Trap_Radius = 1, n_offenders = 1, Guard_Rail = 0.05, n_sources = 1, n_cores = 1, Sd_x = 1, Sd_y = 1)
################################################################################

Extract_Params <- function(Trap_Data, PA_x_grid_cells = 10, PA_y_grid_cells = 10, Time = 1, Trap_Radius = 1, n_offenders = 1, Guard_Rail = 0.05, n_sources = 1, n_cores = 1, Sd_x = 1, Sd_y = 1)
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
			Source_Points_Long <- seq(min(Long_Min_Bound,Long_Max_Bound), max(Long_Min_Bound,Long_Max_Bound), abs(Long_Max_Bound-Long_Min_Bound)/(PA_x_grid_cells - 1))
			Source_Points_Lat <- seq(min(Lat_Min_Bound, Lat_Max_Bound), max(Lat_Min_Bound, Lat_Max_Bound), abs(Lat_Max_Bound-Lat_Min_Bound)/(PA_y_grid_cells - 1))
      # Generate a list of locations of n sources given the number of potential sources
      if(n_sources > 1) {AP_allocation <- t(combnPrim(PA_x_grid_cells*PA_y_grid_cells, n_sources))}
      else{AP_allocation <- "NA"}
      Hits_Only <- subset(Trap_Data, Hits != 0)
			Miss_Only <- subset(Trap_Data, Hits == 0)
			trapping_locations <- cbind(Hits_Only$Latitude, Hits_Only$Longitude)

      params <- list(n_sources = n_sources, PA_x_grid_cells = PA_x_grid_cells, PA_y_grid_cells = PA_y_grid_cells, Sd_x = Sd_x, Sd_y = Sd_y,
      Trap_Radius=Trap_Radius, Time=Time, Source_Points_Long=Source_Points_Long, Source_Points_Lat=Source_Points_Lat,
      AP_allocation=AP_allocation, Hits_Only=Hits_Only, Miss_Only=Miss_Only, Long_Max_Bound = Long_Max_Bound, n_cores = n_cores,
      MnMx_Long= MnMx_Long, MnMx_Lat= MnMx_Lat, n_offenders = n_offenders)
      return(params)
			}

################################################################################
# The "Trap_Po_Parameters" function returns the poisson parameter for each site in the form of an array. Where the individual entries of
# matrix i represent the Poisson parameters of event i.
# example:
# misc_params <- Extract_Params(cbind(1:10,1:10, sample(0:10, 10, replace=T)), PA_x_grid_cells = 10, PA_y_grid_cells = 10, Time = 25, Trap_Radius = 10, n_offenders = 1, Guard_Rail = 0.05, n_sources = 1, n_cores = 1, Sd_x = 2.5, Sd_y = 2.5)
# Trap_Po_Parameters(misc_params)
# rm(misc_params)
################################################################################

Trap_Po_Parameters <- function(Params)
		{
		Po_Array_Hits <- array(NA, c(length(Params$Source_Points_Long), length(Params$Source_Points_Lat), length(Params$Hits_Only$Longitude) ) )
		long_lat_grid <- expand.grid(Params$Source_Points_Long, Params$Source_Points_Lat)

		for(i in 1:length(Params$Hits_Only$Longitude))
		  {
			     Po_hits <- mapply(Poisson_Parameter, Params$Hits_Only$Longitude[i], Params$Hits_Only$Latitude[i], mu_x = long_lat_grid$Var1, mu_y = long_lat_grid$Var2,
										       trap_radius = Params$Trap_Radius, t = Params$Time, Sd_x = Params$Sd_x, Sd_y = Params$Sd_y, hazard_param = Params$n_offenders)
			     Hit_matrix <- matrix(Po_hits, nrow = length(Params$Source_Points_Long), ncol = length(Params$Source_Points_Lat))
			     Po_Array_Hits[, ,i] <-  Hit_matrix
		  }

		if(length(Params$Miss_Only$Longitude)>0)
		  {
		      Po_Array_Miss <- array(NA, c(length(Params$Source_Points_Long), length(Params$Source_Points_Lat), length(Params$Miss_Only$Longitude) ) )

		        for(j in 1:length(Params$Miss_Only$Longitude))
	 		          {
		            Po_miss <- mapply(Poisson_Parameter, Params$Miss_Only$Longitude[j], Params$Miss_Only$Latitude[j], mu_x = long_lat_grid$Var1, mu_y = long_lat_grid$Var2,
											             trap_radius = Params$Trap_Radius, t = Params$Time, Sd_x = Params$Sd_x, Sd_y = Params$Sd_y, hazard_param = Params$n_offenders)
		            Miss_matrix <- matrix(Po_miss, nrow = length(Params$Source_Points_Long), ncol = length(Params$Source_Points_Lat))
		            Po_Array_Miss[, ,j] <-  Miss_matrix
	 		          }
	   	          Po_Param <- list(Po_Array_Hits = Po_Array_Hits, Po_Array_Miss = Po_Array_Miss)
			          return(Po_Param)
		  }
		  else{ Po_Param <- list(Po_Array_Hits = Po_Array_Hits)
		  return(Po_Param)
		      }
		}

###############################################################################
# The "Multisource_probs" function returns the probability surfaces for the hits, misses and both combine.
# Example:
# misc_params <- Extract_Params(cbind(1:10,1:10, sample(0:3, 10, replace=T)), PA_x_grid_cells = 5, PA_y_grid_cells = 5, Time = 5, Trap_Radius = 20, n_offenders = 5, Guard_Rail = 0.05, n_sources = 2, n_cores = 1, Sd_x = 2, Sd_y = 2)
# misc_po_params <- Trap_Po_Parameters(misc_params)
# Multisource_probs(Data_Params = misc_params, Po_Params = misc_po_params)
###############################################################################

Multisource_probs <- function(Data_Params, Po_Params)
		{
    if(length(Data_Params$Miss_Only$Longitude) > 0)
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
					Hits_Prob <- mapply(dpois, Data_Params$Hits_Only$Hits, Po_Params$Po_Array_Hits[i,j,])
					Miss_Prob <- mapply(dpois, Data_Params$Miss_Only$Hits, Po_Params$Po_Array_Miss[i,j,])
					SS_Array_Hits[i,j,] <- Hits_Prob
					SS_Array_Miss[i,j,] <- Miss_Prob
					}
				}
				SS_Hits <- exp(apply(log(SS_Array_Hits), c(1,2), sum))
				SS_Miss <- exp(apply(log(SS_Array_Miss), c(1,2), sum))
				SS_Both <- SS_Hits * SS_Miss
        print(SS_Both)
				Source_Hits <- SS_Hits/sum(SS_Hits)
				Source_Miss <- SS_Miss/sum(SS_Miss)
				Source_Both <- SS_Both/sum(SS_Both)
				Source_Prob <- list(Source_Hits = Source_Hits, Source_Miss = Source_Miss, Source_Both = Source_Both)
			  return(Source_Prob)

			}	 else{
				cluster <- makeCluster(Data_Params$n_cores)
					Sum_hit_po <- matrix(NA, ncol = length(Data_Params$Hits_Only$Longitude), nrow =length(Data_Params$AP_allocation[,1]))
  				Sum_miss_po <- matrix(NA, ncol = length(Data_Params$Miss_Only$Longitude), nrow =length(Data_Params$AP_allocation[,1]))

  				################################################################################################
					print("Calculate Poisson Parameters by summing hazard surfaces")
					Sys.sleep(2)
					print("Poisson Parameters for the Hits")
					Sys.sleep(2)
					for(j in 1:length(Data_Params$Hits_Only[,1]))
						{
							matrix_h <- Po_Params$Po_Array_Hits[,,j]
							Sum_hit_po[, j] <- parRapply(cluster, Data_Params$AP_allocation, FUN = function(x) sum(matrix_h[x]))
							print(j)
						}
						print("Poisson Parameters for the Miss")
						Sys.sleep(2)
					for(j in 1:length(Data_Params$Miss_Only[,1]))
						{
							matrix_m <- Po_Params$Po_Array_Miss[,,j]
							Sum_miss_po[, j] <- parRapply(cluster, Data_Params$AP_allocation, FUN = function(x) sum(matrix_m[x]))
							print(j)
						}

					S_hit_prob <- matrix(NA, ncol = length(Data_Params$Hits_Only$Longitude), nrow =length(Data_Params$AP_allocation[,1]))
					S_miss_prob <- matrix(NA, ncol = length(Data_Params$Miss_Only$Longitude), nrow =length(Data_Params$AP_allocation[,1]))
					print("Calculate probabilities for observing data given our fixed sources")
					Sys.sleep(2)
					print("For the Hits")
					Sys.sleep(2)
					for(i in 1:length(Data_Params$Hits_Only[,1]))
					{
						S_hit_prob[,i] <- mcmapply(dpois, Data_Params$Hits_Only$Hits[i], Sum_hit_po[,i], mc.cores = Data_Params$n_cores)
						print(i)
					}
					S_hit_prob <- log(S_hit_prob)
					S_hit_prob <- exp(apply(S_hit_prob, FUN = sum, MARGIN = 1))
					print("For the misses")
					Sys.sleep(2)
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
				Sys.sleep(2)
				for(i in 1:((Data_Params$PA_x_grid_cells)*(Data_Params$PA_y_grid_cells)))
				{
					print(i)
					a <- parRapply(cluster, S_probs[,1:Data_Params$n_sources], FUN=function (x) any(x == i))
					final_hits[i] <- sum(S_probs[a, Data_Params$n_sources + 1])
					final_miss[i] <- sum(S_probs[a, Data_Params$n_sources + 2])
					final_both[i] <- sum(S_probs[a, Data_Params$n_sources + 3])
				}
				Source_Hits <- matrix(final_hits, ncol = length(Data_Params$Source_Points_Lat), nrow = length(Data_Params$Source_Points_Long))
				Source_Miss <- matrix(final_miss, ncol = length(Data_Params$Source_Points_Lat), nrow = length(Data_Params$Source_Points_Long))
			  Source_Both <- matrix(final_both, ncol = length(Data_Params$Source_Points_Lat), nrow = length(Data_Params$Source_Points_Long))

				Source_Hits <-  Source_Hits/sum(Source_Hits)
				Source_Miss <-  Source_Miss/sum(Source_Miss)
				Source_Both <-  Source_Both/sum(Source_Both)
				Source_Prob <- list( Source_Hits = Source_Hits, Source_Miss = Source_Miss, Source_Both = Source_Both)
				return(Source_Prob)
				stopCluster(cluster)
				}
      }
      else{
        if(Data_Params$n_sources == 1)
    		{
    		SS_Array_Hits <- array(NA, c(length(Po_Params$Po_Array_Hits[,1,1]), length(Po_Params$Po_Array_Hits[1,,1]), length(Po_Params$Po_Array_Hits[1,1,]) ) )

    			for(i in 1:length(Po_Params$Po_Array_Hits[,1,1]))
    				{
    				for(j in 1:length(Po_Params$Po_Array_Hits[1,,1]))
    					{
    					Hits_Prob <- mapply(dpois, Data_Params$Hits_Only$Hits, Po_Params$Po_Array_Hits[i,j,])
    					SS_Array_Hits[i,j,] <- Hits_Prob
    					}
    				}
    				SS_Hits <- exp(apply(log(SS_Array_Hits), c(1,2), sum))
    				Source_Hits <- SS_Hits/sum(SS_Hits)
    			  Source_Prob <- list(Source_Hits = Source_Hits)

    				return(Source_Prob)

    			}	 else{
						cluster <- makeCluster(Data_Params$n_cores)
    					Sum_hit_po <- matrix(NA, ncol = length(Data_Params$Hits_Only$Longitude), nrow =length(Data_Params$AP_allocation[,1]))

      				################################################################################################
    					print("Calculate Poisson Parameters by summing hazard surfaces")
    					Sys.sleep(2)
    					print("Poisson Parameters for the Hits")
    					Sys.sleep(2)
    					for(j in 1:length(Data_Params$Hits_Only[,1]))
    						{
    							matrix_h <- Po_Params$Po_Array_Hits[,,j]
    							Sum_hit_po[, j] <- parRapply(cluster, Data_Params$AP_allocation, FUN = function(x) sum(matrix_h[x]))
    							print(j)
    						}

    					S_hit_prob <- matrix(NA, ncol = length(Data_Params$Hits_Only$Longitude), nrow =length(Data_Params$AP_allocation[,1]))

    					print("Calculate probabilities for observing data given our fixed sources")
    					Sys.sleep(2)
    					print("For the Hits")
    					Sys.sleep(2)
    					for(i in 1:length(Data_Params$Hits_Only[,1]))
    					{
    						S_hit_prob[,i] <- mcmapply(dpois, Data_Params$Hits_Only$Hits[i], Sum_hit_po[,i], mc.cores = Data_Params$n_cores)
    						print(i)
    					}
    					S_hit_prob <- log(S_hit_prob)
    					S_hit_prob <- exp(apply(S_hit_prob, FUN = sum, MARGIN = 1))

              S_probs <- cbind(Data_Params$AP_allocation, S_hit_prob)
    				################################################################################################
    				final_hits <- c()
    				print("Sum over probabilities to obtain likelihood for individual sources")
    				Sys.sleep(2)
    				for(i in 1:((Data_Params$PA_x_grid_cells)*(Data_Params$PA_y_grid_cells)))
    				{
    					print(i)
    					a <- parRapply(cluster, S_probs[,1:Data_Params$n_sources], FUN=function (x) any(x == i))
    					final_hits[i] <- sum(S_probs[a, Data_Params$n_sources + 1])
    				}
    				Source_Hits <- matrix(final_hits, ncol = length(Data_Params$Source_Points_Lat), nrow = length(Data_Params$Source_Points_Long))
    				Source_Hits <-  Source_Hits/sum(Source_Hits)
    				Source_Prob <- list( Source_Hits = Source_Hits)
    				return(Source_Prob)
						stopCluster(cluster)
    				}
      }
			}
###############################################################################
# The "PA_Hitscores" function reports the hitscores of sources taken from the Presence absence profile, much like the geoReportHitscores function
# Example: PA_Hitscores(misc_params, misc_sources, PA_profile)
###############################################################################

PA_Hitscores <- function(params, sources, probability_matrix)
	 {
	 hitscore <- c()
	 for(k in 1:length(sources$source_longitude))
	 {
	 source_matrix <- matrix(NA, ncol = params$PA_y_grid_cells, nrow = params$PA_x_grid_cells)
	 S_x <- rep(sources$source_longitude[k], length(params$Source_Points_Long))
	 diff_x <- which.min(abs(S_x - params$Source_Points_Long))
	 s_grid_long_index <- diff_x

	 S_y <- rep(sources$source_latitude[k], length(params$Source_Points_Lat))
	 diff_y <- which.min(abs(S_y - params$Source_Points_Lat))
	 s_grid_lat_index <- diff_y

	 source_matrix[s_grid_long_index, s_grid_lat_index] = T
	 source_matrix[is.na(source_matrix)] <- F
	 hit_score_index <- order(probability_matrix, decreasing = T)
	 for(t in 1:length(hit_score_index))
	  {
	  if(source_matrix[hit_score_index[t]]==T)
	  {
	  single_hitscore <- t/(params$PA_x_grid_cells*params$PA_y_grid_cells)
	  }
	  else{}
	  }
		hitscore[k] <- single_hitscore
		}
		hitscores <- cbind(sources$source_longitude, sources$source_latitude, hitscore)
	  return(hitscores = hitscores)
	 }

###############################################################################
# The "trap_assignment_data" function takes data generated by the "rDPM" function and converts it into
# trap density data. The main input is the number of sentinels you expect there to be within your search area.
# The function creates most of its parameters using this (e.g. the number of traps)
# Given the rDPM output generate a collection of traps, generate a trap radius for the traps (0.5mean(nnd))
# and construct densities for each trap given the rDPM data falls within the trap radius.
###############################################################################

trap_assignment_data <- function(simulation, sim_params, sources, exp_population = 15, title, n_traps_y = 5, n_traps_x = 5, spacing = NULL)
	{
	#### RANDOM ARRAY OF TRAPS OVER DATA
  if(spacing == "random")
  {
  lower_traps <- floor(0.9*exp_population)
  upper_traps <- ceiling(1.1*exp_population)
  n_traps <- sample(lower_traps:upper_traps, 1)
  longs <- runif(n_traps, sim_params$output$longitude_minMax[1], sim_params$output$longitude_minMax[2])
  lats <- runif(n_traps, sim_params$output$latitude_minMax[1], sim_params$output$latitude_minMax[2])
  trap_loc <- cbind(longs, lats)
  }
  else if(spacing == "uniform")
  {
  #### UNIFORM ARRAY OF TRAPS OVER DATA

  lon_dist <- latlon_to_bearing(sim_params$output$latitude_minMax[1], sim_params$output$longitude_minMax[1], sim_params$output$latitude_minMax[1], sim_params$output$longitude_minMax[2])$gc_dist
  lat_dist <- latlon_to_bearing(sim_params$output$latitude_minMax[1], sim_params$output$longitude_minMax[1], sim_params$output$latitude_minMax[2], sim_params$output$longitude_minMax[1])$gc_dist

  lats <- seq(sim_params$output$latitude_minMax[1], sim_params$output$latitude_minMax[1] + ((lon_dist/lat_dist)*(sim_params$output$latitude_minMax[2] - sim_params$output$latitude_minMax[1])),
             (lon_dist/lat_dist)*(sim_params$output$latitude_minMax[2] - sim_params$output$latitude_minMax[1])/(n_traps_y-1))
	longs <- seq(sim_params$output$longitude_minMax[1], sim_params$output$longitude_minMax[2], (sim_params$output$longitude_minMax[2] - sim_params$output$longitude_minMax[1])/(n_traps_x-1))

  trap_loc <- expand.grid(longs, lats)
  }
  pairwise_traps <- pairwise_distance(trap_loc)
  #### TRAP RADIUS
  detection_TR <- 0.5*mean(pairwise_traps$distance_min, na.rm = TRUE)
  n_traps <- n_traps_x*n_traps_y

  #### PLOT RDPM DATA
  plot(trap_loc, cex = 1.5, ylab = "Latitude", xlab = "Longitude", pch = 0)
  title(main=title)
  points(simulation$longitude, simulation$latitude, pch = 16, col= "green", cex = 1)
  points(simulation$source_lon, simulation$source_lat, col ="blue", pch = 18, cex = 2.5)

	#### DATA WITHIN TRAP RADIUS
	lat_long <- cbind(simulation$longitude, simulation$latitude)
	all_distances <- matrix(NA, nrow = length(simulation$longitude), ncol = length(trap_loc[,1]))

	for(a in 1:length(trap_loc[,1]))
	 {
		   for(b in 1:length(simulation$longitude))
		     {
			        all_distances[b,a] <- latlon_to_bearing(trap_loc[a,2], trap_loc[a,1], simulation$latitude[b], simulation$longitude[b])$gc_dist
		     }
	 }
	Trap_cap_density <- cbind(trap_loc, hit_miss = rep(0, length(trap_loc[,1])))

	for(crime in 1:length(simulation$longitude))
	 {
    row <- all_distances[crime,]
    trap_i <- which(row == min(row))

    if(all_distances[crime, trap_i] < detection_TR)
    {
    Trap_cap_density[trap_i,3] <- Trap_cap_density[trap_i,3] + 1
    }
    else{}
   }
  ##
  Trap_cap_density <- as.data.frame(Trap_cap_density)
  the_hits <- subset(Trap_cap_density, Trap_cap_density$hit_miss >0)
  the_miss <- subset(Trap_cap_density, Trap_cap_density$hit_miss  == 0)
  colnames(the_hits) <- c("longs", "lats", "hit_miss")
  colnames(the_miss) <- c("longs", "lats", "hit_miss")

  if(length(the_hits$longs) > 1 & length(the_miss$longs)>1)
    {
      lon_max <- max(the_hits$longs, the_miss$longs)
      lon_min <- min(the_hits$longs, the_miss$longs)
      lat_max <- max(the_hits$lats, the_miss$lats)
      lat_min <- min(the_hits$lats, the_miss$lats)
      lon_minmax <- c(lon_min, lon_max)
      lat_minmax <- c(lat_min, lat_max)

      sd_hits <- cbind(the_hits$longs, the_hits$lats)
      pw_dist <- pairwise_distance(sd_hits)
      #### STANDARD DEVIATION
      Sd <- 0.5*mean(pw_dist$distance_min, na.rm = TRUE)

      source_hits_dist <- matrix(NA, nrow = length(sources$source_longitude), ncol = length(the_hits[,1]))
      source_miss_dist <- matrix(NA, nrow = length(sources$source_longitude), ncol = length(the_miss[,1]))

      for(j in 1:length(sources$source_longitude))
 	     {
 		       for(k in 1:length(the_hits[,1]))
 		         {
 			            source_hits_dist[j,k] <- latlon_to_bearing(the_hits[k,2], the_hits[k,1], sources$source_latitude[j], sources$source_longitude[j])$gc_dist
 		         }
 	     }
       for(x in 1:length(sources$source_longitude))
        {
            for(y in 1:length(the_miss[,1]))
              {
                   source_miss_dist[x,y] <- latlon_to_bearing(the_miss[y,2], the_miss[y,1], sources$source_latitude[x], sources$source_longitude[x])$gc_dist
              }
        }
      hits_near_sources <- matrix(NA, nrow = length(sources$source_longitude), ncol = 3)
      miss_near_sources <- matrix(NA, nrow = length(sources$source_longitude), ncol = 3)

	     for(Sources in 1:length(sources$source_longitude))
	      {
          hits_one_SD <- length(which(source_hits_dist[Sources ,] <= Sd))
          hits_two_SD <- length(which(source_hits_dist[Sources ,] <= 2*Sd))
          hits_three_SD <- length(which(source_hits_dist[Sources,] <= 3*Sd))
          miss_one_SD <- length(which(source_miss_dist[Sources ,] <= Sd))
          miss_two_SD <- length(which(source_miss_dist[Sources ,] <= 2*Sd))
          miss_three_SD <- length(which(source_miss_dist[Sources,] <= 3*Sd))
          hits_near_sources[Sources,] <- cbind(hits_one_SD, hits_two_SD, hits_three_SD)
          miss_near_sources[Sources,] <- cbind(miss_one_SD, miss_two_SD, miss_three_SD)
        }
        #### DISTANCE BETWEEN 2 SOURCES
        Source_Distances <- NA
        if(length(sources$source_longitude) == 2)
        {
          Source_Distances <- latlon_to_bearing(sources$source_latitude[1], sources$source_longitude[1], sources$source_latitude[2], sources$source_longitude[2])$gc_dist
        }
        else{}
        ##### plot 2 - trap data
        plot(the_hits$longs, the_hits$lats, cex = the_hits$hit_miss, ylab = "Latitude", xlab = "Longitude", col = "green", pch = 19, xlim = lon_minmax, ylim = lat_minmax, main = "The Model's Data")
        points(the_miss$longs, the_miss$lats, cex = 1, col = "red", pch = 4)
        points(simulation$source_lon, simulation$source_lat, col ="blue", pch = 18, cex = 2.5)
	      return(list(Trap_cap_density = Trap_cap_density, Sd = Sd, detection_TR= detection_TR, n_traps = n_traps, hits_near_sources = hits_near_sources, miss_near_sources = miss_near_sources, Source_Distances = Source_Distances))
       }
       else{
       return(list(Trap_cap_density =Trap_cap_density))
       }
    }

################################################################################
# The "PA_simulation" function generates a random data set consisting of trap locations and densities,
# along with a number of sources. This data is run on the PA and DPM and hitscores are extracted amongst all
# parameters and values that are comparable. Note the explicit profiles are not extracted for memory purposes
# The "surface" argument in the geoReportHitscore() function must be changed when using version 2.0.0 of
# RgeopProfile (namely from "posteriorSurface" to "surface")
################################################################################

PA_simulation <- function(trap_spacing = "random", replications = 5, n_offenders = 50, n_sources = n_sources, n_cores = 1, PA_x_grid_cells = 50, PA_y_grid_cells= 50, alpha = 3, sigma= 1, priorMean_longitude = -0.04217481, priorMean_latitude = 51.5235505, guardRail = 0.05)
	 {
   #PA_Hitscores, DPM_hitscores,	Difference, number of traps within 1:3 SD's
		Hitscore_Output <- matrix(NA, nrow = (n_sources*replications), ncol = 13)
    colnames(Hitscore_Output) <- c("DPM_HS", "PA_HS", "DPM-PA", "TRAPS-1_SD", "TRAPS-2_SD", "TRAPS-3_SD", "HITS_1_SD", "HITS_2_SD", "HITS_3_SD", "MISS_1_SD", "MISS_2_SD", "MISS_3_SD", "ALLOC_PTS")
    # Pop_density, actual_hits, actual_misses, potential_offenders, sigma and tau or rDPM, sigma and trap_radiusgenerated from data
    Param_Output <- matrix(NA, nrow = replications, ncol = 10)
    colnames(Param_Output) <- c("Pop_density", "#hits", "#misses", "#traps", "Uniform_sigma_draw","half_meanNND(traps)", "fitted_sigma", "rDPM_tau", "Trap_Radius", "Source_distances")
		Rep <- 1
		while(Rep <= replications)
		{
      random_offenders <- rpois(1, n_offenders)
      print(random_offenders)
      if(random_offenders > 1)
        {
          tau <- runif(1, 1, 2*sigma)
      if(n_sources == 1)
        {
			    sim <- rDPM(random_offenders, priorMean_longitude = priorMean_longitude, priorMean_latitude = priorMean_latitude, alpha = alpha, sigma = sigma, tau = tau)
        }
      else
        {
          sim <- rDPM(5000, priorMean_longitude = priorMean_longitude, priorMean_latitude = priorMean_latitude, alpha = 10, sigma = sigma, tau = tau)
          crime_per_source <- c()
          alloc_leng <- split(sim$group, sim$group)
          for(i in 1:length(unique(sim$group)))
          {
            crime_per_source[i] <- length(alloc_leng[[i]])
          }
          a <- which(crime_per_source > 2*random_offenders)
          pois_one <- rpois(1, floor(random_offenders*0.5))
          pois_two <- random_offenders - pois_one
          one <- sample(1:crime_per_source[a[1]], pois_one)
          two <- sample(1:crime_per_source[a[2]], pois_two)
          sim$longitude <- c(sim$longitude[one], sim$longitude[two])
          sim$latitude <- c(sim$latitude[one], sim$latitude[two])
          sim$source_lon <- c(sim$source_lon[a[1]], sim$source_lon[a[2]])
          sim$source_lat <- c(sim$source_lat[a[1]], sim$source_lat[a[2]])
          sim$group <- c(rep(1, length(one)), rep(2, length(two)))
        }
      point_data <- geoData(sim$longitude, sim$latitude)
      master_params <- geoParams(data = point_data, sigma_mean = 1, sigma_squared_shape = 2, samples= 50000, chains = 200, burnin = 2000, priorMean_longitude = mean(point_data$longitude), priorMean_latitude = mean(point_data$latitude), guardRail = guardRail)
      s <- geoDataSource(sim$source_lon, sim$source_lat)
      Trap_Data <- trap_assignment_data(simulation = sim, sim_params = master_params, sources = s, exp_population = n_offenders, title = "Data Generated Via rDPM", spacing = trap_spacing)
      #############################################################################################################################################################
      if(sum(Trap_Data$hits_near_sources[,3]>4)==n_sources)
      {
      if(length(which(Trap_Data$Trap_cap_density$hit_miss > 1) > 1))
        {
      				### DPM ###
							dpm_hits <- subset(Trap_Data$Trap_cap_density, Trap_Data$Trap_cap_density[,3]>0)
							##### For the dpm repeat the points that already exist given the number of trapped animals. dpm_hits <- dpm_hits[rep(1:nrow(dpm_hits), dpm_hits[,3]),]
							dpm_misses <- subset(Trap_Data$Trap_cap_density, Trap_Data$Trap_cap_density[,3]==0)
							trap_loc_data <- geoData(Trap_Data$Trap_cap_density[,1], Trap_Data$Trap_cap_density[,2])
							hit_data <- geoData(dpm_hits[,1], dpm_hits[,2])
              hit_params <- geoParams(data = hit_data, sigma_mean = 0.9, sigma_squared_shape = 2, samples= 50000, chains = 200, burnin = 2000, priorMean_longitude = mean(hit_data$longitude), priorMean_latitude = mean(hit_data$latitude), guardRail = guardRail)
							hit_params$output$longitude_minMax <- master_params$output$longitude_minMax
							hit_params$output$latitude_minMax <- master_params$output$latitude_minMax
							m <- geoMCMC(data=hit_data, params= hit_params)
              fitted_sigma <- mean(m$sigma)
              dpm_hits <- dpm_hits[rep(1:nrow(dpm_hits), dpm_hits[,3]),]
              hit_data <- geoData(dpm_hits[,1], dpm_hits[,2])
              hit_params <- geoParams(data = hit_data, sigma_mean = fitted_sigma, sigma_var = 0, samples= 50000, chains = 200, burnin = 2000, priorMean_longitude = mean(hit_data$longitude), priorMean_latitude = mean(hit_data$latitude), guardRail = guardRail)
							hit_params$output$longitude_minMax <- master_params$output$longitude_minMax
							hit_params$output$latitude_minMax <- master_params$output$latitude_minMax
							m <- geoMCMC(data=hit_data, params= hit_params)
							### PA ###
							Trap_cap_density <- as.data.frame(Trap_Data$Trap_cap_density)
							Data_parameters <- Extract_Params(Trap_Data = Trap_Data$Trap_cap_density, PA_x_grid_cells = PA_x_grid_cells, PA_y_grid_cells = PA_y_grid_cells, Guard_Rail = guardRail, Trap_Radius = Trap_Data$detection_TR, n_sources = n_sources, n_cores = n_cores, n_offenders = n_offenders, Sd_x = Trap_Data$Sd, Sd_y = Trap_Data$Sd)
							Trap_Poisson_Params <- Trap_Po_Parameters(Data_parameters)
							Source_Probabilities <- Multisource_probs(Data_parameters, Trap_Poisson_Params)
              ### Extract data
              lower_index <- (n_sources*(Rep-1) + 1)
              upper_index <- (n_sources*Rep)
              Hitscore_Output[lower_index:upper_index, 1] <- geoReportHitscores(hit_params, s, m$posteriorSurface)[,3]
              Hitscore_Output[lower_index:upper_index, 2] <- PA_Hitscores(Data_parameters, s, Source_Probabilities$Source_Both)[,3]
              Hitscore_Output[lower_index:upper_index, 3] <- geoReportHitscores(hit_params, s, m$posteriorSurface)[,3] - PA_Hitscores(Data_parameters, s, Source_Probabilities$Source_Both)[,3]
              Hitscore_Output[lower_index:upper_index, 4] <- Trap_Data$hits_near_sources[1:n_sources, 1] + Trap_Data$miss_near_sources[1:n_sources, 1]
              Hitscore_Output[lower_index:upper_index, 5] <- Trap_Data$hits_near_sources[1:n_sources, 2] + Trap_Data$miss_near_sources[1:n_sources, 2]
              Hitscore_Output[lower_index:upper_index, 6] <- Trap_Data$hits_near_sources[1:n_sources, 3] + Trap_Data$miss_near_sources[1:n_sources, 3]
              Hitscore_Output[lower_index:upper_index, 7:9] <- Trap_Data$hits_near_sources[1:n_sources, 1:3]
              Hitscore_Output[lower_index:upper_index, 10:12] <- Trap_Data$miss_near_sources[1:n_sources, 1:3]
              allocations <- c()
              for(i in 1:n_sources)
                {
                  allocations[i] <- length(which(sim$group == i))
                }
              Hitscore_Output[lower_index:upper_index, 13] <- allocations

              Param_Output[Rep, 1] <- random_offenders
              Param_Output[Rep, 2] <- sum(Data_parameters$Hits_Only[,3])
              Param_Output[Rep, 3] <- length(Data_parameters$Miss_Only[,3])
              Param_Output[Rep, 4] <- Trap_Data$n_traps
              Param_Output[Rep, 5] <- sigma
              Param_Output[Rep, 6] <- Trap_Data$Sd
              Param_Output[Rep, 7] <- fitted_sigma
              Param_Output[Rep, 8] <- tau
              Param_Output[Rep, 9] <- Trap_Data$detection_TR
              if(n_sources == 2)
                {
                  Param_Output[Rep, 10] <- Trap_Data$Source_Distances
                }
              else
                {
                  Param_Output[Rep, 10] <- NA
                }
              all_objects <- list(dpm_hits = dpm_hits,dpm_misses =dpm_misses, trap_loc_data = trap_loc_data, hit_data =hit_data, hit_params =  hit_params,m = m, fitted_sigma = fitted_sigma,Trap_cap_density = Trap_cap_density, Data_parameters =  Data_parameters,
                                Trap_Poisson_Params = Trap_Poisson_Params, Source_Probabilities = Source_Probabilities, random_offenders = random_offenders, sim = sim, sigma= sigma, tau = tau)
                                #crime_per_source = crime_per_source, pois_one  =pois_one, pois_two = pois_two, a = a, alloc_leng = alloc_leng)
              Rep <- Rep + 1
					}
      else{}
      }
      else{}
      }
      else{}
      }
		return(list(Hitscores_Output= Hitscore_Output, Param_Output = Param_Output, all_objects = all_objects))
	 }

############################################################# RUN AND FILE SAVES
# LOCAL
# start <-  Sys.time()
#
# par(mfrow =c(1,2))
# simulation <- PA_simulation(replications = 1, n_offenders = 25, n_sources = 2, n_cores = 1, PA_x_grid_cells = 10, PA_y_grid_cells = 10, trap_spacing = "random", alpha = 1, sigma = 1,priorMean_longitude = -0.04217481, priorMean_latitude = 51.5235505, guardRail = 0.05)
# simulation$Hitscores_Output

# end <- Sys.time()
# end - start
#save(simulation, file= "testing_twoS.rdata")

# CLUSTER
# args=commandArgs(trailingOnly=TRUE)
# TwoS_ <- PA_simulation(replications = 1, n_offenders = 25, n_sources = 2, n_cores = 2, PA_x_grid_cells = 30, PA_y_grid_cells = 30,  trap_spacing = "random", priorMean_longitude = -0.04217481, priorMean_latitude = 51.5235505, alpha = 0.5, sigma = 1, guardRail = 0.05)
#
# save(TwoS_, file=paste("TwoS_",args,sep=""))

################################################################################
