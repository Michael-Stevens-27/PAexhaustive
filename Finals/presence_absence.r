#_____________PRESENCE ABSENCE 1.0.0______________##

rm(list = ls())

##_____________ PACKAGES ______________##
library(gRbase) ## For the "combnPrim" function, efficiently listing combinations

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

Extract_Params <- function(Trap_Data, x_grid_cells = 10, y_grid_cells = 10, Time = 1, Guard_Rail = 0.05, Trap_Radius = 0.6, n_sources = 1)
			{
			colnames(Trap_Data) <- c("Longitude", "Latitude", "Hits")
			Trap_Data <- data.frame(Trap_Data)
			MnMx_Long <- c(min(Trap_Data$Longitude), max(Trap_Data$Longitude))
			MnMx_Lat <- c(min(Trap_Data$Latitude), max(Trap_Data$Latitude))
			Long_Max_Bound <- MnMx_Long[2] + Guard_Rail
			Long_Min_Bound <- MnMx_Long[1] - Guard_Rail
			Lat_Max_Bound <-  MnMx_Lat[2] + Guard_Rail
			Lat_Min_Bound <- MnMx_Lat[1] - Guard_Rail
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
			Grid <- expand.grid(1:(x_grid_cells), 1:(y_grid_cells))
			names(Grid) = c("x", "y")
			pairwise <- pairwise_distance(Trap_Data)
			Sd_x <- 0.2*mean(pairwise$distance, na.rm = TRUE)
			Sd_y <- 0.2*mean(pairwise$distance, na.rm = TRUE)

      params <- list(n_sources = n_sources, x_grid_cells = x_grid_cells, y_grid_cells = y_grid_cells, Sd_x = Sd_x, Sd_y = Sd_y,
				    Trap_Radius=Trap_Radius, Time=Time, Anchor_Points_Long=Anchor_Points_Long, Anchor_Points_Lat=Anchor_Points_Lat,
				    Grid = Grid, AP_allocation=AP_allocation, Hits_Only=Hits_Only, Miss_Only=Miss_Only)
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
					Hits_Prob <- mapply(dpois, Data_Params$Hits_Only$Hits, Po_Params$Po_Array_Hits[i,j,])
					Miss_Prob <- mapply(dpois, Data_Params$Miss_Only$Hits, Po_Params$Po_Array_Miss[i,j,])
					SS_Array_Hits[i,j,] <- Hits_Prob
					SS_Array_Miss[i,j,] <- Miss_Prob
					}
				}
				SS_Hits <- apply(SS_Array_Hits, c(1,2), prod)
				SS_Miss <- apply(SS_Array_Miss, c(1,2), prod)
				SS_Both <- SS_Hits * SS_Miss
				SS_Hits <- SS_Hits/sum(SS_Hits)
				SS_Miss <- SS_Miss/sum(SS_Miss)
				SS_Both <- SS_Both/sum(SS_Both)
				SS_Prob <- list(SS_Array_Hits = SS_Array_Hits, SS_Array_Miss = SS_Array_Miss, SS_Hits = SS_Hits, SS_Miss = SS_Miss, SS_Both = SS_Both)

				return(SS_Prob)

			}	 else{
					Sum_hit_po <- matrix(NA, ncol = length(Data_Params$Hits_Only$Longitude), nrow =length(Data_Params$AP_allocation[,1]))
  				Sum_miss_po <- matrix(NA, ncol = length(Data_Params$Miss_Only$Longitude), nrow =length(Data_Params$AP_allocation[,1]))

  				################################################################################################

					for(i in 1:length(Data_Params$AP_allocation[,1]))
					{
						for(j in 1:length(Data_Params$Hits_Only[,1]))
						{
							Sum_hit_po[i, j]  <- sum(Po_Params$Po_Array_Hits[,,j][Data_Params$AP_allocation[i, 1:Data_Params$n_sources]])

						}
						for(j in 1:length(Data_Params$Miss_Only[,1]))
						{
							Sum_miss_po[i, j]  <- sum(Po_Params$Po_Array_Miss[,,j][Data_Params$AP_allocation[i, 1:Data_Params$n_sources]])
						}
					}

					S_hit_prob <- c()
  			  S_miss_prob <- c()

					for(i in 1:length(Data_Params$AP_allocation[,1]))
					{
						S_hit_prob[i] <- prod(mapply(dpois, Data_Params$Hits_Only$Hits, Sum_hit_po[i,]))
						S_miss_prob[i] <- prod(mapply(dpois, Data_Params$Miss_Only$Hits, Sum_miss_po[i,]))
					}

					S_both_prob <- S_hit_prob*S_miss_prob
  				S_probs <- cbind(Data_Params$AP_allocation, S_hit_prob, S_miss_prob, S_both_prob)

				################################################################################################

				final_hits <- c()
				final_miss <- c()
				final_both <- c()

				for(i in 1:length(Data_Params$Grid[,1]))
				{
					a <- apply(TS_probs[,1:Data_Params$n_sources], FUN=function (x) any(x == i), MARGIN=1)
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
				Source_Prob <- list( Source_Hits = Source_Hits, Source_Miss = Source_Miss, Source_Both = Source_Both, S_probs = S_probs )
				return(Source_Prob)
				}
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
				 contour(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat, Probs$Source_Hits, col = "darkgreen", nlevels = 5)
				 contour(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat, Probs$Source_Miss,col = "red",add=TRUE, nlevels = 5)
				 contour(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat, Probs$Source_Both,col = "blue", add=TRUE, nlevels = 5)
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
				#points(DPM_SIM$source_lon, DPM_SIM$source_lat, pch = 16, col = "blue")

				#contour(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat, Two_Source$TwoD_Hits, col = "darkgreen", nlevels = 5)
				#	points(Data_Params$Hits_Only$Longitude, Data_Params$Hits_Only$Latitude , pch = 16, col = "green")
				#points(DPM_SIM$source_lon, DPM_SIM$source_lat, pch = 16, col = "blue")

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

##_____________USER INPUT ______________##

## The only input required is the Longitude/Latitude points of said traps,
## along with the number of hits associated with each (non negative integer)

setwd("/home/mstevens/Desktop/Presence Absence/Presence_Absence-master/Finals/data")

###########################################################################################################################################
My_trap_data <- read.table("Dummy_4_Source.txt", header = FALSE)
###########################################################################################################################################

############### ONE SOURCE ##############
## Extract ALL Parameters from your data
Data_params1S <- Extract_Params(My_trap_data, x_grid_cells = 100, y_grid_cells = 100, Trap_Radius = 5, Guard_Rail = 5, n_sources = 1)

## Compute Poisson Parameters
Trap_Po_Params <- Trap_Po_Parameters(Data_params1S)

## Compute probability matrices
#Single_Source_Prob <- Single_Source_Probability(Data_params1S, Trap_Po_Params)
Single_Source_Prob <- Multisource_probs(Data_params1S, Trap_Po_Params)

## Contour plot of the matrices
x11()
plot_sources(Data_params1S, Single_Source_Prob)

############## TWO SOURCE ###########
## Extract ALL Parameters from your data

Data_params2S <- Extract_Params(My_trap_data, x_grid_cells = 6, y_grid_cells = 6, Guard_Rail = 3, Trap_Radius = 1, n_sources = 4)

## Compute Poisson Parameters
Trap_Po_Params <- Trap_Po_Parameters(Data_params2S)

## Compute probability matrices
Two_Source_Prob <- Multisource_probs(Data_params2S, Trap_Po_Params)

## Contour plot of the matrices
x11()
plot_sources(Data_params2S, Two_Source_Prob)

###################
###################
