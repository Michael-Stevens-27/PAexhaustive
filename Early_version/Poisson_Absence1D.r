rm(list = ls()) #remove all objects
library(flux)
###########sentinel_sites <-as.numeric(sentinel_sites)

#Specify the Trap locations, the number of 'hits' in this trap and the various anchor points

anchor_points <- seq(-5,10,0.05)

#hits and miss'

trap_location1 <- seq(-10, 0, 1)
trap_location2 <- seq(0.5, 2.5, 0.1)
trap_location3 <- seq(3, 10, 1)
trap_locations <- c(trap_location1, trap_location2, trap_location3)
length(trap_locations)
hit_or_miss <- c(1,1,1,1,1,1,1,1,1,1,2,4,5,7,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,1,1,0,0)
trap_data <- cbind(trap_locations, hit_or_miss)

#Specify the time period
t <- 10

#Define the trap radius and Standard deviation of source distribution
trap_radius <- 0.4
sd <- 3

#####################################
#####  PROBABILITIES - LOOPING  #####
#####################################

trap_po_params_combine <- c()
hit_probabilities <- c()
miss_probabilities <- c()

combine_hits <- c()
combine_miss <- c()
combine_both <- c()
combine <- c()

trap_data[,1]
for(ck in 1:length(anchor_points))
	{

	for(ci in 1:length(trap_data[,1]))
		{
		trap_po_params_combine[ci] <- (2*(trap_radius)*t*dnorm(trap_data[ci,1], mean = anchor_points[ck], sd))
		}

	for(cp in 1:length(trap_data[,1]))
	{
		if(trap_data[cp,2] == 0)
			{
			miss_probabilities[cp] <- exp(-trap_po_params_combine[cp])
			}
			else
			{
			hit_probabilities[cp] <- (((trap_po_params_combine[cp])^(trap_data[cp,2]))*exp(-trap_po_params_combine[cp]))/(gamma(trap_data[cp,2]+1))
			}
	}
	combine_hits[ck] <- prod(hit_probabilities, na.rm = TRUE)
	combine_miss[ck] <- prod(miss_probabilities, na.rm = TRUE)
	combine[ck] <- prod(hit_probabilities, miss_probabilities, na.rm = TRUE)
	}

#######################################
#####  PROBABILITIES - FUNCTIONS  #####
#######################################
#
#anchorMat <- matrix(	rep(anchor_points, length(trap_locations)),
#				ncol = length(anchor_points),
#				nrow = length(trap_locations))
#
#
#trap_po_params_combine <- matrix(NA, length(trap_locations), length(anchor_points))
#hit_probabilities <- matrix(NA, 1, length(anchor_points))
#miss_probabilities <- matrix(NA, 1, length(anchor_points))
#
#combine_hits <- matrix(NA, 1, length(anchor_points))
#combine_miss <- matrix(NA, 1, length(anchor_points))
#combine_both <- matrix(NA, 1, length(anchor_points))
#combine <- matrix(NA, 1, length(anchor_points))
#
#
#poisson_value <- function(lambda, k)
#           	     	{
#		     	return((lambda^k)*(exp(-lambda))/(gamma(k + 1)))
#		      }
#
#poisson_parameter <- function(trap_loc, trap_radius, t, mean, sd)
#                     	{
#				return(2*trap_radius*t*dnorm(trap_loc, mean = mean, sd = sd))
#                        }
#
#a <- poisson_parameter(trap_data[,1], trap_radius, t, anchor_points, sd)
#poisson_value(a, hit_or_miss)
#
#
#
##########################
######  NORMALISE   ######
##########################

norm_combine_hits <- combine_hits/auc(anchor_points, combine_hits)
norm_combine_miss <- combine_miss/auc(anchor_points, combine_miss)
norm_combine <- combine/auc(anchor_points, combine)

######################
######  PLOTS   ######
######################


#plot(anchor_points, norm_combine,type="l", col = "blue", xlab = "Possible Sources", ylab = "Probability")
#zeros <- rep(0, length(trap_data[,1]))
#text(trap_data[,1], zeros, labels = trap_data[,2])
#plot(anchor_points, norm_combine_hits,type="l", col = "green", xlab = "Possible Sources", ylab = "Probability")
#text(trap_data[,1], zeros, labels = trap_data[,2])
#plot(anchor_points, norm_combine_miss,type="l", col = "red", xlab = "Possible Sources", ylab = "Probability")
#text(trap_data[,1], zeros, labels = trap_data[,2])


x11()
plot(anchor_points, norm_combine,type="l", col = "blue", xlab = "Possible Sources", ylab = "Probability")
lines(anchor_points, norm_combine_hits, col="green")
lines(anchor_points, norm_combine_miss, col="red")

dev.copy(png,'myplot.png')
dev.off()
