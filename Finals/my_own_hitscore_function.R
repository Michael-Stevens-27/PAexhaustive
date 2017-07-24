
index <- order(-my_matrix)
my_matrix
Data_parameters <- Extract_Params(dhitmiss, x_grid_cells = 500, y_grid_cells = 500, Guard_Rail = 0.05, Trap_Radius = 0.015, n_sources = 1, n_cores = 1)
Trap_Poisson_Params <- Trap_Po_Parameters(Data_parameters)
Source_Probabilities <- Multisource_probs(Data_parameters, Trap_Poisson_Params)
source_location <- s

Michael_geoReportHitscores <- function(params, sources, probability_matrix)
 {
 source_matrix <- matrix(NA, ncol = params$y_grid_cells, nrow = params$x_grid_cells)
    for(i in 1:(length(Data_parameters$Anchor_Points_Long)-1))
    {
      a <- abs(source_location$source_longitude[1] - Data_parameters$Anchor_Points_Long[i])
      b <- abs(source_location$source_longitude[1] - Data_parameters$Anchor_Points_Long[i+1])
      if(a < abs(Data_parameters$Anchor_Points_Long[i+1] - Data_parameters$Anchor_Points_Long[i]) & b < abs(Data_parameters$Anchor_Points_Long[i+1] - Data_parameters$Anchor_Points_Long[i]))
      {
        s_grid_long <- Data_parameters$Anchor_Points_Long[i]
        s_grid_long_index <- i
      }
      else{
      }
    }
    for(i in 1:(length(Data_parameters$Anchor_Points_Lat)-1))
    {
      a <- abs(source_location$source_latitude[1] - Data_parameters$Anchor_Points_Lat[i])
      b <- abs(source_location$source_latitude[1] - Data_parameters$Anchor_Points_Lat[i+1])
      if(a < abs(Data_parameters$Anchor_Points_Lat[i+1] - Data_parameters$Anchor_Points_Lat[i]) & b < abs(Data_parameters$Anchor_Points_Lat[i+1] - Data_parameters$Anchor_Points_Lat[i]))
      {
        s_grid_lat <- Data_parameters$Anchor_Points_Long[i]
        s_grid_lat_index <- i
      }
      else{
      }
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
source_matrix <- Michael_geoReportHitscores(Data_parameters, source_location, Source_Probabilities$Source_Both)
