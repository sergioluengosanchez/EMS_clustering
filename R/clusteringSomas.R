#' Function to set which variables of the somas dataset are linear or directional
#'
#' This function gets a dataset and assigns the label directional depending on the prefix theta and phi.
#'
#' @param data a data.frame with the morphological measures over the surface of the somas
#'
#' @return a list of two elements, the first one is the data with the columns ordered amd the second element is a binary vector denoting if a column is directional or not
#'
#' @examples
#' data<-read.csv("data/soma_dataset_short.csv")[,-1]
#' dataset <- set_directional_variables(data)
#'
#' @export
set_directional_variables<-function(data)
{
  idx_dir <- c(grep("theta\\_*",colnames(data)),grep("inst\\_Phi\\_*",colnames(data)))
  ordered_colnames <- c(colnames(data)[idx_dir],colnames(data)[setdiff(1:ncol(data),idx_dir)])
  data <- data[,ordered_colnames]
  is_dir <- matrix(0,ncol=ncol(data)+1,nrow=1)
  is_dir[1:length(idx_dir)] <- 1
  is_dir <- as.data.frame(is_dir)
  colnames(is_dir) <- c(ordered_colnames, "class")

  row.names(data) <- NULL
  return(list(data = data, is_dir = is_dir))
}

#' Function that performs clustering using the Extended Mardia-Sutton distribution
#'
#' This function gets a dataset and assigns the label directional depending on the prefix theta and phi.
#'
#' @param data a data.frame
#' @param is_dir a binary vector denoting if the variable at each column is directional or not
#' @param num_parents is an interger denoting the maximum number of parents for a variable
#' @param directional_parents is a boolean value denoting if the linear variables can have directional parents. Setting false is equivalent to assume independence among linear and directional variables
#' @param num_clusters is an intenger number denoting number of clusters that the algorithm should try for model selection. If the number of clusters is 0 the last column is assumed to be an initialisation of each somas in each cluster
#' @param cores denotes the number of cores used for parallel computation. Set 1 core to avoid parallel computation
#'
#' @return Structure as an adjacency matrix, weights and score of the best model found by the algorithm
#'
#' @examples
#' model <- EMS_clustering(dataset$data, dataset$is_dir, num_clusters = 3, cores = 10)
#'
#' @export
EMS_clustering <- function(data, is_dir, num_clusters = 10, num_parents = 5, directional_parents = T, cores = 1)
{
  #Generate a list where each element corresponds to a number of clusters
  SEM_cl_results<-list()
  count <- 1

  #If the user want to run the algorithm in parallel
  if(cores>1)
  {
    cores <- min(c(num_clusters-1,cores))
    cl<- makeCluster(cores, outfile="")
    registerDoParallel(cl)
  }

  #If the initiation should be at random
  if(num_clusters>0)
  {
    if(cores>1)#Parallel running
    {
      SEM_cl_results <-foreach(j=2:num_clusters, .packages=c("SomaSimulation"))%dopar%{
        c_hybrid_SEM(data_r = data.matrix(data), is_directional_r = data.matrix(is_dir), num_clusters = j, max_parents = num_parents, directional_parents = directional_parents)
      }
    }else{ #One process running
      for(i in 2:num_clusters)
      {
        SEM_result <- c_hybrid_SEM(data_r = data.matrix(data), is_directional_r = data.matrix(is_dir), num_clusters = i, max_parents = num_parents, directional_parents = directional_parents)
        SEM_cl_results[[count]] <- SEM_result
        count <- count + 1
      }
    }
  }else{#If the initiation is done according to a initial assignment provided by the user
    SEM_result <- c_hybrid_SEM(data_r = data.matrix(data), is_directional_r = data.matrix(is_dir), num_clusters = 0, max_parents = num_parents, directional_parents = directional_parents)
    SEM_cl_results[[count]] <- SEM_result
    count <- count + 1
  }


  if(cores>1)
  {
    stopCluster(cl)
  }

    #Get the BN from previous results.
    model <- getBestBN(SEM_cl_results, data, is_dir)
    return(model)
}
