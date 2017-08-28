rm(list = ls()) #remove all objects
library(flux)
library(gRbase) ## For the "combnPrim" function, efficiently listing combinations
                ## Initial install, run commands:
								## source("http://bioconductor.org/biocLite.R")
								## biocLite("gRbase")
library(parallel)
library(RgeoProfile)
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
		dist <- latlon_to_bearing(y1, x1, y2, x2)$gc_dist
		distance[i,j] <- dist
		}
	}
	diag(distance)  <- NA
	distance[distance==0] <- NA
	distance_min <- apply(distance, 1, min, na.rm = TRUE)
	distance_min[distance_min =="Inf"] <- NA
	SD_TR <- list(distance = distance)
	return(SD_TR)
}

###########################################################################################################################################
#
# The "Extract_Params" function returns a list of parameters for an "n_sources" source model. Should "n_sources" be set to anything greater
# than one, "Extract_Params" will produce an additional parameter "AP_allocation," a list ilustrating the number of ways of choosing
# "n_sources" sources from the total number of grid cells.
#
###########################################################################################################################################

Extract_Params <- function(Simulated_Data, Trap_Data, x_grid_cells = 10, y_grid_cells = 10, Time = 1, Guard_Rail = 0.05, Trap_Radius = 0.6, n_sources = 1, n_cores = 1)
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
			Simulated_Data <- cbind(Simulated_Data$latitude, Simulated_Data$longitude)
			pairwise <- pairwise_distance(Simulated_Data)
			Sd_x <- 0.5*mean(pairwise$distance, na.rm = TRUE)
			Sd_y <- 0.5*mean(pairwise$distance, na.rm = TRUE)

      params <- list(n_sources = n_sources, x_grid_cells = x_grid_cells, y_grid_cells = y_grid_cells, Sd_x = Sd_x, Sd_y = Sd_y,
				    Trap_Radius=Trap_Radius, Time=Time, Anchor_Points_Long=Anchor_Points_Long, Anchor_Points_Lat=Anchor_Points_Lat,
				    AP_allocation=AP_allocation, Hits_Only=Hits_Only, Miss_Only=Miss_Only, Long_Max_Bound = Long_Max_Bound, n_cores = n_cores,
					  MnMx_Long= MnMx_Long, MnMx_Lat= MnMx_Lat)
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
				param <- co_efficient*exp( -0.5*( (latlon_to_bearing(y, x, y, mu_x)$gc_dist)^2/(Sd_x*Sd_x) + (latlon_to_bearing(y, mu_x, mu_y, mu_x)$gc_dist)^2/(Sd_y*Sd_y) ) )
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
		long_lat_grid <- expand.grid(Params$Anchor_Points_Long, Params$Anchor_Points_Lat)

		for(i in 1:length(Params$Hits_Only$Longitude))
		{
			Po_hits <- mapply(Poisson_Parameter, Params$Hits_Only$Longitude[i], Params$Hits_Only$Latitude[i], mu_x = long_lat_grid$Var1, mu_y = long_lat_grid$Var2,
												trap_radius = Params$Trap_Radius, t = Params$Time, Sd_x = Params$Sd_x, Sd_y = Params$Sd_y)
			Hit_matrix <- matrix(Po_hits, nrow = length(Params$Anchor_Points_Long), ncol = length(Params$Anchor_Points_Lat))
			Po_Array_Hits[, ,i] <-  Hit_matrix

		}
		if(length(Params$Miss_Only$Longitude)>0)
		{
		Po_Array_Miss <- array(NA, c(length(Params$Anchor_Points_Long), length(Params$Anchor_Points_Lat), length(Params$Miss_Only$Longitude) ) )

		for(j in 1:length(Params$Miss_Only$Longitude))
	 		{
		  Po_miss <- mapply(Poisson_Parameter, Params$Miss_Only$Longitude[j], Params$Miss_Only$Latitude[j], mu_x = long_lat_grid$Var1, mu_y = long_lat_grid$Var2,
											 trap_radius = Params$Trap_Radius, t = Params$Time, Sd_x = Params$Sd_x, Sd_y = Params$Sd_y)
		  Miss_matrix <- matrix(Po_miss, nrow = length(Params$Anchor_Points_Long), ncol = length(Params$Anchor_Points_Lat))
		  Po_Array_Miss[, ,j] <-  Miss_matrix
	 		}
	   	Po_Param <- list(Po_Array_Hits = Po_Array_Hits, Po_Array_Miss = Po_Array_Miss)
			return(Po_Param)
		}
		else{
		Po_Param <- list(Po_Array_Hits = Po_Array_Hits)
		return(Po_Param)
		}
		}

###########################################################################################################################################
#
# The "Multisource_probs" function returns the porbability for source locations. Given the user specifies the number of sources expected.
#
###########################################################################################################################################

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
					Sys.sleep(3)
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
					Sys.sleep(3)
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
				Sys.sleep(3)
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
    					Sum_hit_po <- matrix(NA, ncol = length(Data_Params$Hits_Only$Longitude), nrow =length(Data_Params$AP_allocation[,1]))

      				################################################################################################
    					print("Calculate Poisson Parameters by summing hazard surfaces")
    					Sys.sleep(3)
    					print("Poisson Parameters for the Hits")
    					Sys.sleep(3)
    					for(j in 1:length(Data_Params$Hits_Only[,1]))
    						{
    							matrix_h <- Po_Params$Po_Array_Hits[,,j]
    							Sum_hit_po[, j] <- parRapply(cluster, Data_Params$AP_allocation, FUN = function(x) sum(matrix_h[x]))
    							print(j)
    						}

    					S_hit_prob <- matrix(NA, ncol = length(Data_Params$Hits_Only$Longitude), nrow =length(Data_Params$AP_allocation[,1]))

    					print("Calculate probabilities for observing data given our fixed sources")
    					Sys.sleep(3)
    					print("For the Hits")
    					Sys.sleep(3)
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
    				Sys.sleep(3)
    				for(i in 1:((Data_Params$x_grid_cells)*(Data_Params$y_grid_cells)))
    				{
    					print(i)
    					a <- parRapply(cluster, S_probs[,1:Data_Params$n_sources], FUN=function (x) any(x == i))
    					final_hits[i] <- sum(S_probs[a, Data_Params$n_sources + 1])
    				}
    				Source_Hits <- matrix(final_hits, ncol = length(Data_Params$Anchor_Points_Lat), nrow = length(Data_Params$Anchor_Points_Long))
    				Source_Hits <-  Source_Hits/sum(Source_Hits)
    				Source_Prob <- list( Source_Hits = Source_Hits)
    				return(Source_Prob)
    				}
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

##########################################################################################################################################
##########################################################################################################################################

Michael_geoReportHitscores <- function(params, sources, probability_matrix)
	 {
	 source_matrix <- matrix(NA, ncol = params$y_grid_cells, nrow = params$x_grid_cells)
	    for(i in 1:(length(params$Anchor_Points_Long)-1))
	    {
	      a <- abs(sources$source_longitude[1] - params$Anchor_Points_Long[i])
	      b <- abs(sources$source_longitude[1] - params$Anchor_Points_Long[i+1])
	      if(a < abs(params$Anchor_Points_Long[i+1] - params$Anchor_Points_Long[i]) & b < abs(params$Anchor_Points_Long[i+1] - params$Anchor_Points_Long[i]))
	      {
	        s_grid_long <- params$Anchor_Points_Long[i]
	        s_grid_long_index <- i
	      }
	      else{
	      }
	    }
	    for(j in 1:(length(params$Anchor_Points_Lat)-1))
	    {
	      a <- abs(sources$source_latitude[1] - params$Anchor_Points_Lat[j])
	      b <- abs(sources$source_latitude[1] - params$Anchor_Points_Lat[j+1])
	      if(a < abs(params$Anchor_Points_Lat[j+1] - params$Anchor_Points_Lat[j]) & b < abs(params$Anchor_Points_Lat[j+1] - params$Anchor_Points_Lat[j]))
	      {
	        s_grid_lat <- params$Anchor_Points_Long[j]
	        s_grid_lat_index <- j
	      }
	      else if(sources$source_latitude[1] > params$Anchor_Points_Lat[j+1])
        {
        s_grid_lat <- params$Anchor_Points_Long[j+1]
        s_grid_lat_index <- j
	      }
        else {}
	    }
	 source_matrix[s_grid_long_index, s_grid_lat_index] = T
	 source_matrix[is.na(source_matrix)] <- F
	 hit_score_index <- order(-probability_matrix)
	 for(i in 1:length(hit_score_index))
	  {
	  if(source_matrix[hit_score_index[i]]==T)
	  {
	  hitscore <- i/(params$x_grid_cells*params$y_grid_cells)
	  }
	  else{}
	  }
	  return(hitscore)
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

trap_assignment_data <- function(simulation, sim_params, n_traps_x, n_traps_y, p)
	{
	#### makee a grid of traps
	lon_dist <- latlon_to_bearing(sim_params$output$latitude_minMax[1], sim_params$output$longitude_minMax[1], sim_params$output$latitude_minMax[1], sim_params$output$longitude_minMax[2])$gc_dist
	lat_dist <- latlon_to_bearing(sim_params$output$latitude_minMax[1], sim_params$output$longitude_minMax[1], sim_params$output$latitude_minMax[2], sim_params$output$longitude_minMax[1])$gc_dist

	#lats <- seq(sim_params$output$latitude_minMax[1], sim_params$output$latitude_minMax[1] + ((lon_dist/lat_dist)*(sim_params$output$latitude_minMax[2] - sim_params$output$latitude_minMax[1])), (lon_dist/lat_dist)*(sim_params$output$latitude_minMax[2] - sim_params$output$latitude_minMax[1])/(n_traps_y-1))
	#lats <- seq(sim_params$output$latitude_minMax[1], sim_params$output$latitude_minMax[2], (lon_dist/lat_dist)*(sim_params$output$latitude_minMax[2] - sim_params$output$latitude_minMax[1])/(n_traps_y-1))

	longs <- seq(sim_params$output$longitude_minMax[1], sim_params$output$longitude_minMax[2], (sim_params$output$longitude_minMax[2] - sim_params$output$longitude_minMax[1])/(n_traps_x-1))
	lats <- seq(sim_params$output$latitude_minMax[1], sim_params$output$latitude_minMax[2], (sim_params$output$latitude_minMax[2] - sim_params$output$latitude_minMax[1])/(n_traps_y-1))
	trap_loc <- expand.grid(longs, lats)
	detection_TR <- min(0.4999*(latlon_to_bearing(trap_loc[1,2], trap_loc[1,1], trap_loc[2,2], trap_loc[2,1])$gc_dist), 0.4999*(latlon_to_bearing(trap_loc[1,2], trap_loc[1,1], trap_loc[n_traps_x+1,2], trap_loc[n_traps_x+1,1])$gc_dist))
	plot(trap_loc, xlim = c(sim_params$output$longitude_minMax[1],sim_params$output$longitude_minMax[2]), ylim = c(sim_params$output$latitude_minMax[1], sim_params$output$latitude_minMax[2]), cex = 2, ylab = "Latitude", xlab = "Longitude")
	title(main=p)
	points(simulation$source_lon, simulation$source_lat, col ="blue", pch = 15, cex = 2)
  points(simulation$longitude, simulation$latitude, pch = 18, col= "green", cex = 2)

	##### distance metric
	lat_long <- cbind(simulation$longitude, simulation$latitude)
	all_distances <- matrix(NA, nrow = length(simulation$longitude), ncol = length(trap_loc[,1]))

	for(a in 1:length(trap_loc[,1]))
	{
		for(b in 1:length(simulation$longitude))
		{
			all_distances[b,a] <- latlon_to_bearing(trap_loc[a,2], trap_loc[a,1], simulation$latitude[b], simulation$longitude[b])$gc_dist
		}
	}
	Trap_data <- cbind(trap_loc, hit_miss = rep(NA, length(trap_loc[,1])))

	for(crime in 1:length(trap_loc[,1]))
	{
	detected_by <- length(which(all_distances[,crime] <= detection_TR))
  Trap_data[crime,3] <- detected_by
	}
	return(list(Trap_data = Trap_data, detection_TR= detection_TR))
	}

###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
# _______________________________________________________
# SIM DATA AND SET PARAMS
# _______________________________________________________
# generate simulated data
# sim data for single source using rDPM

################################################################################
############################## SIMULATION ######################################
################################################################################

DPM_hitscores3x3 <- c()
PA_Hitscores3x3 <- c()
Difference3x3 <- c()

PA_simulation(replications = 100, sim_params, no_x_traps = 5, no_y_traps = 5, DPM_params)
{
	for(Rep in 1:replications)
	{
	sim <- rDPM(15, priorMean_longitude = -0.04217481, priorMean_latitude = 51.5235505, alpha = 0, sigma =1, tau = 3)
	point_data <- geoData(sim$longitude, sim$latitude)
	master_params <- geoParams(data = point_data, sigma_mean = 1, sigma_squared_shape = 2, samples= 50000, chains = 200, burnin = 2000, priorMean_longitude = mean(point_data$longitude), priorMean_latitude = mean(point_data$latitude), guardRail = 0.05)
	Trap_Data <- trap_assignment_data(sim, master_params, no_x_traps, no_y_traps, Rep)

	if(length(which(Trap_Data$Trap_data$hit_miss>0))>1)
	{
	### DPM ###
	dpm_hits <- subset(Trap_Data$Trap_data, Trap_Data$Trap_data[,3]>0)
	dpm_hits <- dpm_hits[rep(1:nrow(dpm_hits), dpm_hits[,3]),]
	dpm_misses <- subset(Trap_Data$Trap_data, Trap_Data$Trap_data[,3]==0)
	trap_loc_data <- geoData(Trap_Data$Trap_data[,1], Trap_Data$Trap_data[,2])
	s <- geoDataSource(sim$source_lon, sim$source_lat)
	hit_data <- geoData(dpm_hits[,1], dpm_hits[,2])
	hit_params <- geoParams(data = hit_data, sigma_mean = 1, sigma_squared_shape = 2, samples= 50000, chains = 200, burnin = 2000, priorMean_longitude = mean(hit_data$longitude), priorMean_latitude = mean(hit_data$latitude), guardRail = 0.05)
	hit_params$output$longitude_minMax <- master_params$output$longitude_minMax
	hit_params$output$latitude_minMax <- master_params$output$latitude_minMax
	m <- geoMCMC(data=hit_data, params= hit_params)

	### PA ###
	Trap_Data <- as.data.frame(Trap_Data)
	Data_parameters <- Extract_Params(sim, Trap_Data, x_grid_cells = 50, y_grid_cells = 50, Guard_Rail = 0.05, Trap_Radius = Trap_Data$detection_TR, n_sources = 1, n_cores = 1)
	Trap_Poisson_Params <- Trap_Po_Parameters(Data_parameters)
	Source_Probabilities <- Multisource_probs(Data_parameters, Trap_Poisson_Params)

	DPM_hitscores[p] <- geoReportHitscores(hit_params, s, m$posteriorSurface)[,3]
	PA_Hitscores[p] <- Michael_geoReportHitscores(Data_parameters, s, Source_Probabilities$Source_Both)
	Difference[p] <- geoReportHitscores(hit_params, s, m$posteriorSurface)[,3] - Michael_geoReportHitscores(Data_parameters, s, Source_Probabilities$Source_Both)
	}
	}
	else{p <- p-1}
	}
}

# _______________________________________________________
# RUN PRESENCE/ABSENCE WITH THESE PARAMS
# _______________________________________________________
#start.time <- Sys.time()
## Extract ALL Parameters from your data
## Compute Poisson Parameters
## Compute probability matrices in parallel
#cluster <- makeCluster(3)
#clusterExport(cluster, c("Data_parameters", "Trap_Poisson_Params"))
#stopCluster(cluster)
#end.time <- Sys.time()
#time.taken <- end.time - start.time
#time.taken
###########################################################################################################################################
# PLOTTING
#hits <- geoData(Data_parameters$Hits_Only$Longitude, Data_parameters$Hits_Only$Latitude)
#misses <- geoDataSource(Data_parameters$Miss_Only$Longitude, Data_parameters$Miss_Only$Latitude)

# the map alone
#x11(display = "Original")
#geoPlotMap(data = hit_data, source = s, params = hit_params, breakPercent = seq(0, 10, 1), mapType = "roadmap", contourCols =c("darkred", "red", "orange", "yellow"),#
#          crimeCol = "darkgreen", crimeCex = 2, sourceCol = "blue", sourceCex = 2, surface = m$geoProfile)

#Source_Probabilities$Source_Both <- expandMatrix(Source_Probabilities$Source_Both, 500, 500)
#Source_Probabilities$Source_Hits <- expandMatrix(Source_Probabilities$Source_Hits, 500, 500)
#x11()
#geoPlotMap(data = hit_data, source = s, params = master_params, breakPercent = seq(0, 10, 1), mapType = "roadmap", contourCols =c("darkred", "red", "orange", "yellow"),
#						crimeCol = "darkgreen", crimeCex = 2, sourceCol = "red", sourceCex = 2, surface = rank(-Source_Probabilities$Source_Both))

# hits
#x11()
#geoPlotMap(data = hit_data, params = master_params, source = s, breakPercent = seq(0, 10, 1), mapType = "roadmap", contourCols =c("red", "orange", "yellow", "white"),
#           crimeCol = "darkgreen", crimeCex = 5, sourceCol = "red", sourceCex = 5, surface = rank(-Source_Probabilities$Source_Hits))

# misses
#x11("Misses")
#geoPlotMap(data = hits, source = misses, params = params, breakPercent = seq(50, 100, 10), mapType = "roadmap", contourCols =c("red", "orange", "yellow", "white"),
#           crimeCol = "darkgreen", crimeCex = 5, sourceCol = "red", sourceCex = 5, surface = rank(-Source_Probabilities$Source_Miss))

###########################################################################################################################################
