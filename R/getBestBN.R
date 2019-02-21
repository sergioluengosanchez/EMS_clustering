#Given SEM clustering results get the parameters and structure of the best model found
getBestBN <- function(SEM_cl_results, data, is_dir)
{
  best_BIC <- -Inf
  best_model <- c()
  for(i in 1:length(SEM_cl_results))
  {
    if(!is.nan(SEM_cl_results[[i]]$BIC))
    {
      if(SEM_cl_results[[i]]$BIC > best_BIC)
      {
        best_model <- SEM_cl_results[[i]]
        best_BIC <- SEM_cl_results[[i]]$BIC
      }
    }
  }
  if(best_model$weights==-1)
  {
    return(list(priori = -1, structure = -1, gaussParams = -1, vMParams = -1, weights = -1, BIC = -1))
  }else{
    priori <- colSums(best_model$weights) / nrow(best_model$weights)
    structure <- getBNStructure(best_model, data)
    params <- getParams(best_model, structure, data, is_dir)

    return(list(priori = priori, structure = structure, gaussParams = params$gaussParams, vMParams = params$vMParams, weights = best_model$weights, BIC = best_BIC))
  }
}

getBNStructure <- function(best_model, data)
{
  var_names <- colnames(data)
  var_names <- var_names[var_names!=""]
  adj_matrix <- matrix(c(data.matrix(best_model$adj_m)),ncol=length(var_names),nrow=length(var_names),dimnames=list(var_names,var_names))
  structure <- empty.graph(var_names)
  amat(structure)<-adj_matrix
  return(structure)
}

getParams <- function(best_model, structure, data, is_dir)
{
  gaussParams<-list()
  #Get the params of the Gaussian variables for each k cluster
  for(k in 1:ncol(best_model$weights))
  {
    gaussParams[[k]] <- maximizationGauss(structure, data, best_model$weights[ ,k], is_dir)
  }

  #Get the params of the von Mises variables for each cluster
  vMParams<-NULL
  if(sum(is_dir)>0)
  {
    vMParams <- maximizationVM(as.data.frame(data[,which(is_dir==1)]), best_model$weights, colnames(data)[which(is_dir==1)])
  }

  return(list(gaussParams = gaussParams, vMParams = vMParams))
}
