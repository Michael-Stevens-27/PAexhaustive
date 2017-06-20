rm(list = ls()) #remove all objects
library(flux)


# _______________________________________________________
# SIM DATA AND SET PARAMS
# _______________________________________________________
library(RgeoProfile)
# generate simulated data
sim <- rDPM(50, priorMean_longitude = -0.04217481, priorMean_latitude = 51.5235505, alpha = 1, sigma = 1, tau = 12)
sim
###########################################################################################################################################

# allocate hits and misses to crimes
hits <- sample(0:1,length(sim$longitude),replace=TRUE)

trap_data <-  cbind(sim$longitude,sim$latitude, hits)
colnames(trap_data) <- c("x", "y", "hits")
head(trap_data)

# _______________________________________________________
# RUN AS DPM TO EXTRACT PARAMS
# _______________________________________________________
# convert
d <- geoData(sim$longitude,sim$latitude)
s <- geoDataSource(sim$source_lon,sim$source_lat)

# params
params <- geoParams(data = d, sigma_mean = 1, sigma_squared_shape = 2, samples= 50000, chains = 20, burnin = 1000, priorMean_longitude = mean(d$longitude), priorMean_latitude = mean(d$latitude), guardRail = 0.05)
params$output$longitude_cells <- 50
params$output$latitude_cells <- 50

m <- geoMCMC(data=d,params= params)
#####################################################################################################################################

# _______________________________________________________
# FUNCTIONS
# _______________________________________________________
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

Extract_Params <- function(Trap_Data, x_grid_cells = 10, y_grid_cells = 10, Time = 1, Guard_Rail = 0.05, Trap_Radius = 0.6)
			{
			colnames(Trap_Data) <- c("Longitude", "Latitude", "Hits")
			Trap_Data <- data.frame(Trap_Data)
			MnMx_Long <- c(min(Trap_Data$Longitude), max(Trap_Data$Longitude))
			MnMx_Lat <- c(min(Trap_Data$Latitude), max(Trap_Data$Latitude))
			Long_Max_Bound <- MnMx_Long[2] + Guard_Rail
			Long_Min_Bound <- MnMx_Long[1] - Guard_Rail
			Lat_Max_Bound <-  MnMx_Lat[2] + Guard_Rail
			Lat_Min_Bound <- MnMx_Lat[1] - Guard_Rail
			Anchor_Points_Long <- seq(min(Long_Min_Bound,Long_Max_Bound), max(Long_Min_Bound,Long_Max_Bound), abs((Long_Max_Bound-Long_Min_Bound)/x_grid_cells))
			Anchor_Points_Lat <- seq(min(Lat_Min_Bound, Lat_Max_Bound), max(Lat_Min_Bound, Lat_Max_Bound), abs(Lat_Max_Bound-Lat_Min_Bound)/y_grid_cells)
			Anchor_Point_Grid <- expand.grid(Anchor_Points_Long, Anchor_Points_Lat)
			names(Anchor_Point_Grid) = c("Long", "Lat")
			AP_allocation <- cbind(rep(1:((x_grid_cells + 1)*(y_grid_cells + 1)),each=((x_grid_cells + 1)*(y_grid_cells + 1))),
						     rep(1:((x_grid_cells + 1)*(y_grid_cells + 1)),(x_grid_cells + 1)*(y_grid_cells + 1)))
			AP_allocation <- unique(t(apply(AP_allocation, 1, sort)))
			Hits_Only <- subset(Trap_Data, Hits != 0)
			Miss_Only <- subset(Trap_Data, Hits == 0)
			Grid <- expand.grid(1:(x_grid_cells + 1), 1:(y_grid_cells + 1))
			names(Grid) = c("x", "y")
			pairwise <- pairwise_distance(Trap_Data)
			Sd_x <- 0.5*mean(pairwise$distance, na.rm = TRUE)
			Sd_y <- 0.5*mean(pairwise$distance, na.rm = TRUE)
			#Trap_Radius <- (pairwise$TR)
		  #Time <- (Sd_x*Sd_y)/(Trap_Radius^2)
			params <- list(x_grid_cells = x_grid_cells, y_grid_cells = y_grid_cells, Sd_x = Sd_x, Sd_y = Sd_y,
				    Trap_Radius=Trap_Radius, Time=Time, Anchor_Points_Long=Anchor_Points_Long, Anchor_Points_Lat=Anchor_Points_Lat,
				    Anchor_Point_Grid=Anchor_Point_Grid, Grid = Grid, AP_allocation=AP_allocation, Hits_Only=Hits_Only, Miss_Only=Miss_Only)
			return(params)
			}

###########################################################################################################################################

Poisson_Parameter <- function(x, y, mu_x, mu_y, trap_radius, t, Sd_x, Sd_y)
                     	{
				co_efficient <- 4*(trap_radius^2)*t*(2*pi*Sd_x*Sd_y)^(-1)
				param <- co_efficient*exp( -0.5*( (x - mu_x)^2/(Sd_x*Sd_x) + (y - mu_y)^2/(Sd_y*Sd_y) ) )
				return(param)
                        }

########################################################################################################################################

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

Single_Source_Probability <- function(Data_Params, Trap_Params)
		{
		SS_Array_Hits <- array(NA, c(length(Trap_Params$Po_Array_Hits[,1,1]), length(Trap_Params$Po_Array_Hits[1,,1]), length(Trap_Params$Po_Array_Hits[1,1,]) ) )
		SS_Array_Miss <- array(NA, c(length(Trap_Params$Po_Array_Miss[,1,1]), length(Trap_Params$Po_Array_Miss[1,,1]), length(Trap_Params$Po_Array_Miss[1,1,]) ) )

			for(i in 1:length(Trap_Params$Po_Array_Hits[,1,1]))
				{
				for(j in 1:length(Trap_Params$Po_Array_Hits[1,,1]))
					{
					Hits_Prob <- mapply(dpois, Data_Params$Hits_Only$Hits, Trap_Params$Po_Array_Hits[i,j,])
					Miss_Prob <- mapply(dpois, Data_Params$Miss_Only$Hits, Trap_Params$Po_Array_Miss[i,j,])
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
	     }

###########################################################################################################################################

DS_Po_Params <- function(Data_Params, Po_Params)
		{
	Sum_hit_po <- matrix(NA, ncol = length(Data_Params$Hits_Only$Longitude), nrow =length(Data_Params$AP_allocation[,1]))
  Sum_miss_po <- matrix(NA, ncol = length(Data_Params$Miss_Only$Longitude), nrow =length(Data_Params$AP_allocation[,1]))

  ################################################################################################

	for(i in 1:length(Data_Params$AP_allocation[,1]))
	{

  Sum_hit_po[i, ] <- mapply( sum, Po_Params$Po_Array_Hits[   Data_Params$Grid[Data_Params$AP_allocation[i,1],1],
																	    Data_Params$Grid[Data_Params$AP_allocation[i,1],2] , ],
				          Po_Params$Po_Array_Hits[   Data_Params$Grid[Data_Params$AP_allocation[i,2],1],
																			Data_Params$Grid[Data_Params$AP_allocation[i,2],2], ])

	Sum_miss_po[i, ] <- mapply(sum, Po_Params$Po_Array_Miss[   Data_Params$Grid[Data_Params$AP_allocation[i,1],1],
																	    Data_Params$Grid[Data_Params$AP_allocation[i,1],2], ],
			            Po_Params$Po_Array_Miss[   Data_Params$Grid[Data_Params$AP_allocation[i,2],1],
																	    Data_Params$Grid[Data_Params$AP_allocation[i,2],2], ])
	}

	TS_hit_prob <- c()
  TS_miss_prob <- c()

	for(i in 1:length(Data_Params$AP_allocation[,1]))
	{
	TS_hit_prob[i] <- prod(mapply(dpois, Data_Params$Hits_Only$Hits, Sum_hit_po[i,]))
	TS_miss_prob[i] <- prod(mapply(dpois, Data_Params$Miss_Only$Hits, Sum_miss_po[i,]))
	}

	TS_both_prob <- TS_hit_prob*TS_miss_prob
  TS_probs <- cbind(Data_Params$AP_allocation, TS_hit_prob, TS_miss_prob, TS_both_prob)

	################################################################################################

	final_hits <- c()
	final_miss <- c()
	final_both <- c()
	for(i in 1:max(TS_probs[,1]))
	{
	a <- subset(TS_probs, TS_probs[,1] == i | TS_probs[,2] == i)
	final_hits[i] <- sum(a[,3])
	b <- subset(TS_probs, TS_probs[,1] == i | TS_probs[,2] == i)
	final_miss[i] <- sum(b[,4])
	c <- subset(TS_probs, TS_probs[,1] == i | TS_probs[,2] == i)
	final_both[i] <- sum(c[,5])
	}
	TwoD_Hits <- matrix(final_hits, ncol = length(Data_Params$Anchor_Points_Lat), nrow = length(Data_Params$Anchor_Points_Long))
	TwoD_Miss <- matrix(final_miss, ncol = length(Data_Params$Anchor_Points_Lat), nrow = length(Data_Params$Anchor_Points_Long))
	TwoD_Both <- matrix(final_both, ncol = length(Data_Params$Anchor_Points_Lat), nrow = length(Data_Params$Anchor_Points_Long))

	### NORMALISE ??????

	TwoD_Hits <-  TwoD_Hits/sum(TwoD_Hits)
	TwoD_Miss <-  TwoD_Miss/sum(TwoD_Miss)
	TwoD_Both <-  TwoD_Both/sum(TwoD_Both)
	TwoD_Prob <- list( TwoD_Hits = TwoD_Hits, TwoD_Miss = TwoD_Miss, TwoD_Both = TwoD_Both, TS_probs = TS_probs )
	return(TwoD_Prob)
}

###########################################################################################################################################

plot1source <- function(Data_Params, SS_Prob)
		   {
				 #x11()
				 contour(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat, SS_Prob$SS_Hits, col = "darkgreen", nlevels = 5)
				 contour(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat, SS_Prob$SS_Miss,col = "red",add=TRUE, nlevels = 5)
				 contour(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat, SS_Prob$SS_Both,col = "blue", add=TRUE, nlevels = 5)
				 points(Data_Params$Hits_Only$Longitude, Data_Params$Hits_Only$Latitude , pch = 16, col = "green")
				 points(Data_Params$Miss_Only$Longitude, Data_Params$Miss_Only$Latitude , pch = 16, col = "red")

				 #x11()
				 #par(mfrow=c(1,3))
				 #persp(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat,  SS_Prob$SS_Hits, col = "green")
				 #persp(Data_Params$Anchor_Points_Long,Data_Params$Anchor_Points_Lat,  SS_Prob$SS_Miss, col = "red")
				 #persp(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat,  SS_Prob$SS_Both, col = "blue")
			 }

plot2source <- function(Data_Params, Two_Source, DPM_SIM)
		   {
			#par(mfrow=c(2,2))

			contour(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat, Two_Source$TwoD_Hits, col = "darkgreen", nlevels = 5)
			contour(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat, Two_Source$TwoD_Miss,col = "red",add=TRUE, nlevels = 3)
			contour(Data_Params$Anchor_Points_Long, Data_Params$Anchor_Points_Lat, Two_Source$TwoD_Both,col = "blue",add=TRUE, nlevels = 5)
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

###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################

# _______________________________________________________
# RUN PRESENCE/ABSENCE WITH THESE PARAMS
# _______________________________________________________
trap_data <- as.data.frame(trap_data)

all_params <- Extract_Params(trap_data,  x_grid_cells = 49, y_grid_cells = 49)

po_params <- Trap_Po_Parameters(all_params)

one_D_probs <- Single_Source_Probability(all_params, po_params)

probs <- DS_Po_Params(all_params, po_params)

###########################################################################################################################################
dim(probs$TwoD_Both)

hits <- geoData(all_params$Hits_Only$Longitude, all_params$Hits_Only$Latitude)

misses <- geoDataSource(all_params$Miss_Only$Longitude,all_params$Miss_Only$Latitude)

s <- geoDataSource(sim$source_lon,sim$source_lat)

######

# the map alone

x11()
geoPlotMap(data = d, source = s, params = params, breakPercent = seq(0, 20, 2), mapType = "roadmap", contourCols =c("red", "orange", "yellow", "white"),
           crimeCol = "black", crimeCex = 2, sourceCol = "red", sourceCex = 2, surface = m$geoProfile)
######

# both
x11()
geoPlotMap(data = hits, source = misses, params = params, breakPercent = seq(0, 100, 10), mapType = "roadmap", contourCols =c("red", "orange", "yellow", "white"),
           crimeCol = "darkgreen", crimeCex = 5, sourceCol = "red", sourceCex = 5, surface = rank(probs$TwoD_Both))
# hits
x11()
geoPlotMap(data = hits, params = params, breakPercent = seq(0, 100, 10), mapType = "roadmap", contourCols =c("red", "orange", "yellow", "white"),
           crimeCol = "darkgreen", crimeCex = 5, sourceCol = "red", sourceCex = 5, surface = rank(probs$TwoD_Hits))


# misses
x11()
geoPlotMap(data = hits, source = misses, params = params, breakPercent = seq(0, 100, 10), mapType = "roadmap", contourCols =c("red", "orange", "yellow", "white"),
           crimeCol = "darkgreen", crimeCex = 5, sourceCol = "red", sourceCex = 5, surface = rank(-probs$TwoD_Miss))

#x11()
#perspGP(surface=t(probs$TwoD_miss),aggregate_size=5,surface_type="prob",phiGP=70,thetaGP=-10)


#x11()
#perspGP(surface=t(probs$TwoD_Both),aggregate_size=3,surface_type="prob")
#############################################################################################################################################
