rm(list = ls()) # Remove all objects


# LOAD ANY OTHER USEFUL FUNCTIONS
# _______________________________________________________
# allows rotation reflection etc of matrix
# note that "reflect_diagonal" only works with square matrices
matrix_manipulation <- function(my_matrix,my_operation)
	{
		if(my_operation=="reflect_y") {return(my_matrix[,ncol(my_matrix):1])}
		if(my_operation=="reflect_x") {return(my_matrix[nrow(my_matrix):1,])}
		if(my_operation=="rotate_180") {return(my_matrix[nrow(my_matrix):1,ncol(my_matrix):1])}
		if(my_operation=="reflect_diag") {return(t(my_matrix))}
		if(my_operation=="rotate_90") {return(t(my_matrix)[,nrow(my_matrix):1])}
		if(my_operation=="rotate_-90") {return(t(my_matrix)[ncol(my_matrix):1,])}
	}

################################################
## NUMBER OF TRAPS AND SOURCES, TRAP BORDERS ###
################################################
xminMax <- c(0, 6)
yminMax <- c(0, 3)
n_traps <- 10
k <- 2

#############################################
## LOCATIONS OF TRAPS - AND THEIR CONTENTS ##
#############################################
trap_locations <- matrix(runif(2*n_traps, xminMax[1], xminMax[2]), ncol = 2)
hits <- sample(0:1, n_traps, replace = TRUE)
trap_data <-  cbind(trap_locations, hits)
colnames(trap_data) <- c("x", "y", "hits")

################
## DUMMY DATA ##
################
#trap_locations <- matrix(c(1,5,1,5), ncol=2)
#hits <- c(1,1)
#trap_data <-  cbind(trap_locations, hits)
#colnames(trap_data) <- c("x", "y", "hits")

##########################################
## EXTRACT INDIVIDUAL HIT AND MISS DATA ##
##########################################
trap_data <- as.data.frame(trap_data)
just_hits <- subset(trap_data, hits != 0)
just_miss <- subset(trap_data, hits == 0)

#################################
#### ANCHOR POINT LOCATIONS ##### 
#################################
anchor_points_x <- seq(xminMax[1], xminMax[2], 0.5)
anchor_points_y <- seq(yminMax[1], yminMax[2], 0.5)
anchor_points <- expand.grid(anchor_points_x, anchor_points_y)

##############################################
## TRAP RADIUS, STANDARD DEVIATION AND TIME ##
##############################################
#trap_radius <- trap_radii[tr]
trap_radius <- 0.26 
Sd_x <- 1.1
Sd_y <- 1.1
t <- 100

############################################################################
##### PERMUTATIONS OF K INDISTINGUISHABLE SOURCES OUT OF N GRID SQUARE #####
############################################################################
assign_ap_position <- t(combn(length(anchor_points[,1]), k, simplify = TRUE))

###########################
## CREATE BLANK MATRICES ##
###########################

miss_probabilities <- c()
hit_probabilities <- c()
hit_miss_probabilities <- c()

total_hits <- c()
total_miss <- c()

sums_both <- c()
sums_hits <- c()
sums_miss <- c()

prods_both <- c()
prods_hits <- c()
prods_miss <- c()

final_both <- c()
final_hits <- c()
final_miss <- c()


#######################
#####  FUNCTIONS  #####
#######################
poisson_value <- function(lambda = 5 , k)
           	     	{
		     	return((lambda^k)*(exp(-lambda))/(factorial(k)))
		      }		

poisson_parameter <- function(x, y, mu_x, mu_y, trap_radius, t, Sd_x, Sd_y)
                     	{
				p <- 4*(trap_radius^2)*t*(2*pi*Sd_x*Sd_y)^(-1)
				return(p*exp( -0.5*( (x - mu_x)^2/(Sd_x*Sd_x) + (y - mu_y)^2/(Sd_y*Sd_y) ) ) )				
                        }	

#####################################################
#####  2D - PROBABILITIES - LOOPING - FUNCTIONS #####
#####################################################
for(ck in 1:length(anchor_points_x))
	{
	for(cb in 1:length(anchor_points_y))
		{
		trap_po_params_hits <- mapply(poisson_parameter, just_hits[,1], just_hits[,2], mu_x = anchor_points_x[ck], 
							   mu_y = anchor_points_y[cb], trap_radius = trap_radius, t = t, Sd_x = Sd_x, Sd_y = Sd_y)
		trap_po_params_miss <- mapply(poisson_parameter, just_miss[,1], just_miss[,2], mu_x = anchor_points_x[ck], 
							   mu_y = anchor_points_y[cb], trap_radius = trap_radius, t = t, Sd_x = Sd_x, Sd_y = Sd_y)

		miss_probabilities <- mapply(poisson_value, trap_po_params_miss, 0) 
		hit_probabilities <- mapply(poisson_value, trap_po_params_hits, 1)
		hit_miss_probabilities <- c(hit_miss_probabilities, hit_probabilities, miss_probabilities)
		total_hits <- c(total_hits, hit_probabilities)
		total_miss <- c(total_miss, miss_probabilities)
		}
	}
####################################################

hit_sources <- matrix(total_hits, ncol = length(anchor_points[,1]), nrow = length(just_hits[,1]))
#hit_sources <- matrix_manipulation(hit_sources, "rotate_-90")

miss_sources <- matrix(total_miss, ncol = length(anchor_points[,1]), nrow = length(just_miss[,1]))
#miss_sources <- matrix_manipulation(miss_sources, "rotate_-90")

both_sources <- matrix(hit_miss_probabilities, ncol = length(anchor_points[,1]), nrow = length(trap_data[,1]))
#both_sources <- matrix_manipulation(both_sources, "rotate_-90")

#####################################################

for(i in 1:length(assign_ap_position[,1]))
	{ 
	#i <- i - 1
	#lap <- length(assign_ap_position[,1])
	sums_hits <- mapply(sum, hit_sources[,assign_ap_position[i,1]], hit_sources[,assign_ap_position[i,2]])
	prods_hits[i] <- prod(sums_hits)	
	sums_miss <- mapply(sum, miss_sources[,assign_ap_position[i,1]], miss_sources[,assign_ap_position[i,2]])
	prods_miss[i] <- prod(sums_miss)
	sums_both <- mapply(sum, both_sources[,assign_ap_position[i,1]], both_sources[,assign_ap_position[i,2]])
	prods_both[i] <- prod(sums_both)		
	}

####################################################

combo_prob_hit <- cbind(assign_ap_position, prods_hits)
combo_prob_hit <- as.data.frame(combo_prob_hit)

combo_prob_miss <- cbind(assign_ap_position, prods_miss)
combo_prob_miss <- as.data.frame(combo_prob_miss)

combo_prob <- cbind(assign_ap_position, prods_both)
combo_prob <- as.data.frame(combo_prob)

for(i in 1:length(anchor_points[,1]))
	{
	a <- subset(combo_prob_hit, V1 == i | V2 == i)
	final_hits[i] <- sum(a[,k+1])
	b <- subset(combo_prob_miss, V1 == i | V2 == i)
	final_miss[i] <- sum(b[,k+1])
	c <- subset(combo_prob, V1 == i | V2 == i)
	final_both[i] <- sum(c[,k+1])
	}

#####################################################

ult_final_hits <- matrix(final_hits, nrow = length(anchor_points_x), ncol = length(anchor_points_y), byrow = TRUE)
#ult_final_hits <- matrix_manipulation(ult_final_hits, "rotate_-90")
ult_final_miss <- matrix(final_miss, nrow = length(anchor_points_x), ncol = length(anchor_points_y), byrow = TRUE)
#ult_final_miss <- matrix_manipulation(ult_final_miss, "rotate_-90")
ult_final_both <- matrix(final_both, nrow = length(anchor_points_x), ncol = length(anchor_points_y), byrow = TRUE)
#ult_final_both <- matrix_manipulation(ult_final_both, "rotate_-90")

#####################################################



##########################
######  NORMALISE   ######
##########################

two_source_hits <- ult_final_hits/sum(ult_final_hits)
two_source_miss <- ult_final_miss/sum(ult_final_miss)
two_source_both <- ult_final_both/sum(ult_final_both)

######################
######  PLOTS   ######
######################
## CONTOUR
#
windows()
#quartz()
contour(anchor_points_x, anchor_points_y, two_source_hits, col="darkgreen", nlevels = 5)
#contour(anchor_points_x, anchor_points_y, two_source_miss, col="red", add=TRUE,nlevels = 5)
#contour(anchor_points_x, anchor_points_y, two_source_both, col="blue", add= TRUE nlevels = 5)
points(just_hits[,1], just_hits[,2] , pch = 16, col = "green")
points(just_miss[,1], just_miss[,2] , pch = 16, col = "red")

## PERSPECTIVE
#windows()
#quartz()
#dev.new(width = 20, height = 15)
#par(mfrow=c(1,3))
#persp(two_source_hits,phi=30,theta=45,border="darkgreen")
#persp(two_source_miss,phi=30,theta=45,border="red")
#persp(two_source_both,phi=30,theta=45,border="blue")

###########################################################################################################################

