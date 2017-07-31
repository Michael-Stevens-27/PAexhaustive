Michael_geoReportHitscores <- function(params, sources, probability_matrix)
 {
 source_matrix <- matrix(NA, ncol = params$y_grid_cells, nrow = params$x_grid_cells)
    for(i in 1:(length(Data_parameters$Anchor_Points_Long)-1))
    {
      a <- abs(sources$source_longitude[1] - Data_parameters$Anchor_Points_Long[i])
      b <- abs(sources$source_longitude[1] - Data_parameters$Anchor_Points_Long[i+1])
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
      a <- abs(sources$source_latitude[1] - Data_parameters$Anchor_Points_Lat[i])
      b <- abs(sources$source_latitude[1] - Data_parameters$Anchor_Points_Lat[i+1])
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
