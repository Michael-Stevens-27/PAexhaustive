rm(list = ls())
n <- 10
points <- rnorm(n, 0, 1)

my_matrix <- matrix(points, ncol = 2)

pairwise_distance <- function(points)
{
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

pairwise <- pairwise_distance(my_matrix)

#pairwise_min <- apply(pairwise, 1, min, na.rm = TRUE)
#pairwise_min[pairwise_min =="Inf"] <- NA
#mean(pairwise_min, na.rm = TRUE)

###### diff between sd and tr?
##### sd 
##### half the mean nearest neighbour distance
##### tr 
##### minimum of matrix halved?
##### half the mean nearest neighbour 
##### minimum of each row? 
##### half the average (minimum of each row)






