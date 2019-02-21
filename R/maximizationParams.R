#' Compute E-step
#'
#' Compute responsabilities of each cluster on each instance
#'
#' @param priori is a vector that represents the a priori probability of an instance belonging to an instace without further information
#' @param gaussParams is a bn.fit object that saves the parameters of the multivariate gaussian Bayesian network
#' @param vMParams is a list saving objects of vonMises Class for each cluster with their parameters mu and kappa
#' @param linearData a dataframe of the variables that represents linear variables
#' @param angularData a dataframe of the variables that represents angular variables
#'
#' @return list with two entries, the first one is the matrix of responsabilities the second one is de log-likelihood of the model
expectation<-function(priori,gaussParams,vMParams,linearData,angularData)
{
  K<-length(gaussParams)
  logLikelihood<-matrix(0,nrow=nrow(linearData),ncol=K)

  for(i in 1:K)
  {
    for(j in 1:ncol(angularData))
    {
      logLikelihood[,i]<-logLikelihood[,i]+vMParams[[i]][[j]]$density(angularData[,j],log=T)
    }
    logLikelihood[,i]<-logLik(gaussParams[[i]], linearData,by.sample = T)+logLikelihood[,i]
  }

  #Compute the weights or responsabilities of each cluster in each instance with logSumExp to avoid underflow
  prioriLoglik <- sweep(logLikelihood, 2, log(priori), "+")
  denominator <- logSumExp(prioriLoglik)
  weights <- exp(sweep(prioriLoglik, 1, denominator, "-"))

  return (list(weights = weights, log_lik = sum(denominator)))
}


#' M-step Gaussian variables
#'
#' Compute M-step for Gaussian variables
#'
#' @param gaussParams is a bn.fit object that saves the parameters of the multivariate gaussian Bayesian network
#' @param linearData a dataframe of the variables that represents linear variables
#' @param responsabilities a matrix representing the probability of each instance to belong to each cluster. Computed from e-step
#'
#' @return an optimization of gaussParams
#'
#' @noRd
maximizationGauss<-function(structure, data, responsabilities, is_dir)
{
  linear_data <- get_linear_dataset(data, is_dir)
  gaussParams <- list()
  var_names <- names(structure$nodes)
  for(i in 1:length(var_names))
  {
    #If the variable is not directional
    if(is_dir[i] == 0)
    {
      parents <- which(var_names %in% structure$nodes[[i]]$parents)
      if (length(parents) == 0) #If it has not got parents
      {
        mean <- (matrix(responsabilities, nrow = 1) %*% matrix(data[,i], ncol = 1)) / sum(responsabilities)
        coefs <- c("(Intercept)" = mean)
        resid <- data[ ,i] - as.vector(mean)
        sd <- sqrt((matrix(responsabilities, nrow = 1) %*% matrix((resid)^2, ncol = 1)) / sum(responsabilities))
        gaussParams[[var_names[i]]] <- list(coef = as.numeric(coefs), sd = as.numeric(sd))
      }else{ #If there is at least one parent
        node_data <- get_node_data(data, linear_data, is_dir, parents, i)
        coefs <- estimateCoefficients(node_data, responsabilities)
        resid <- data.matrix(node_data)[,ncol(node_data)] - matrix(coefs, nrow = 1) %*% t(as.matrix(cbind(1, node_data[,1:(ncol(node_data)-1)])))
        sd <- sqrt((matrix(responsabilities, nrow = 1) %*% matrix((resid)^2, ncol = 1)) / sum(responsabilities))
        gaussParams[[var_names[i]]] <- list(coef = as.numeric(coefs), sd = as.numeric(sd))
      }#ELSE
    }
  }#FOR


  return(gaussParams)
}

#' M-step von Mises
#'
#' Compute M-step for von Mises variables
#'
#' @param vMParams is a list saving objects of vonMises Class for each cluster with their parameters mu and kappa
#' @param angularData a dataframe of the variables that represents angular variables
#' @param responsabilities a matrix representing the probability of each instance to belong to each cluster. Computed from e-step
#'
#' @return an optimization of gaussParams
#'
#' @noRd
maximizationVM <- function(angularData, responsabilities, var_names)
{
  num_weighted_instances <- colSums(responsabilities)
  K <- ncol(responsabilities)
  params <- list()
  clusterParams <- list()
  varList <- list()

  for(i in 1:K)
  {
    for(j in 1:ncol(angularData))
    {
      params$mu <- atan2((responsabilities[,i] %*% sin(angularData[,j])), (responsabilities[,i] %*% cos(angularData[,j])))
      params$kappa <- A1inv((responsabilities[,i] %*% cos(angularData[,j] - as.vector(params$mu))) / num_weighted_instances[i])
      varList[[var_names[j]]] <- vonMises$new(params$mu, params$kappa)
    }
    clusterParams[[i]] <- varList
  }
  return (clusterParams)
}


#Estimate coefficients of the parents for the maximization step of Gaussian bayesian network
#A Linear Gaussian distribution is represented as N(x|B_0+B_1*x_p1+B_2*x_p2,sigma)
#B_0+B_1*x_p1+B_2*x_p2 is the mean of the gaussian, B_ is a coefficient
#and x_p are the column vector of each parent of the actual node
#To compute the mle, the partial derivative of each coefficient is computed giving as result a linear system
estimateCoefficients <- function(data, responsabilities)
{
  num_coef <- ncol(data)
  #Initialize coefficient and independent term matrix
  coefficient_matrix <- matrix(NA, nrow = num_coef, ncol = num_coef)
  independent_term_matrix <- matrix(rep(NA, num_coef), ncol = 1)

  #The column of the data that corresponds to the actual node is represented as a row to multiply faster
  node_column_values <- matrix(data[ ,num_coef], nrow = 1)
  #The same as before but with the parents nodes. We added a column of 1s to multiply the intercept
  parent_columns_values <- as.matrix(cbind(1, data[ ,1:(num_coef-1)]))

  #Responsabilities (weights) multiplying each parent
  comb_respons_parents <- matrix(sweep(data.matrix(data[ ,1:(num_coef-1)]), 1, responsabilities, "*"), ncol = (num_coef - 1))#Se computa el vector de responsabilidades por los datos del padre que se esta estudiando.

  #Compute the ecuation of the partial derivative of the intercept
  coefficient_matrix[1, ] <- matrix(responsabilities, nrow = 1) %*% parent_columns_values
  independent_term_matrix[1] <- node_column_values %*% matrix(responsabilities, ncol = 1)


  #Compute the ecuation of the partial derivative of each one of the parents
  for(parent in 2:num_coef)
  {
    coefficient_matrix[parent, ] <- matrix(comb_respons_parents[ ,parent-1], nrow = 1) %*% parent_columns_values
    independent_term_matrix[parent] <- node_column_values %*% comb_respons_parents[ ,parent-1]
  }

  #Solve for the coefficients
  coefs<-solve(coefficient_matrix, independent_term_matrix, tol = 1e-40)
  return(coefs)
}


