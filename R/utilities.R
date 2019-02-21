#Function to transform directional variables X to cos(X), sin(X). They are used to compute
#the Mardia-Shutton model efficiently.
get_linear_dataset <- function(data, is_directional)
{
  num_cols <- sum(is_directional + 1)
  linear_dataset <- matrix(0, nrow(data), num_cols)
  
  col_position <- 1
  for(i in 1:ncol(data))
  {
    if(is_directional[i] == 1)
    {
      linear_dataset[,col_position] = cos(data[,i])
      col_position <- col_position + 1
      linear_dataset[,col_position] = sin(data[,i])
    }else{
      linear_dataset[,col_position] = data[,i];
    }
    col_position <- col_position + 1
  }
  
  return(linear_dataset);
}


get_node_data <- function(data, linear_data, is_directional, parents_id, node_id)
{
  ncols <- 1;
  for(parent_id in parents_id)
  {
    if(is_directional[parent_id] == 1)
    {
      ncols <- ncols + 2;
    }else{
      ncols <- ncols + 1;
    }
  }
  
  node_data <- matrix(0, nrow(data), ncols)

  col_position <- 1;
  for(parent_id in parents_id)
  {
    if(is_directional[parent_id] == 1)
    {
      node_data[ ,col_position] <- linear_data[ ,parent_id + (parent_id-1)]
      col_position <- col_position + 1
      node_data[ ,col_position] <- linear_data[ ,parent_id + (parent_id-1) + 1]
    }else{
      node_data[ ,col_position] <- data[ ,parent_id]
    }
    col_position <- col_position + 1
  }
  
  node_data[ ,col_position] <- data[ ,node_id]
  
  return(node_data);
}
