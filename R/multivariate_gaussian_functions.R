#Compute Joint mean
# compute_joint_mean<-function(structure, vM_p, gauss_p, vM_q, gauss_q, is_dir)
compute_joint_mean<-function(structure, vM_p, gauss_p, is_dir)
{
  col_names <- node.ordering(structure)
  dir_names <- names(vM_p)
  dir_names <- col_names[col_names %in% dir_names]
  dir_names_trig <- c(paste0(dir_names,"-cos"),paste0(dir_names,"-sin"))
  linear_names <- names(gauss_p)
  linear_names <- col_names[col_names %in% linear_names]
  col_names_trig <- c(dir_names_trig,linear_names)
  node_name_pos <- gsub("-cos|-sin","",col_names_trig)

  # col_names <- col_names[col_names %in% names(gauss_p)]
  n <- length(col_names_trig)

  #Definition of joint mean vector
  joint_mean_p <- matrix(data=0,nrow = length(col_names_trig), ncol = 1)
  names(joint_mean_p) <- col_names_trig
  # joint_mean_q <- joint_mean_p

  #Insert directional means
  for(j in 1:length(dir_names))
  {
    idx <- grep(dir_names[j],dir_names_trig)
    joint_mean_p[idx[1]] <- cos(vM_p[[dir_names[j]]]$get_mean())
    joint_mean_p[idx[2]] <- sin(vM_p[[dir_names[j]]]$get_mean())
    # joint_mean_q[idx[1]] <- cos(vM_q[[dir_names[j]]]$get_mean())
    # joint_mean_q[idx[2]] <- sin(vM_q[[dir_names[j]]]$get_mean())
  }

  #For each variable compute its mean
  for(i in 1:length(linear_names)){
    #Get the parents and the coefficients
    variable <- linear_names[i]
    parents <- structure$nodes[[variable]]$parents
    coeffs_p <- gauss_p[[variable]]$coef
    # coeffs_q <- gauss_q[[variable]]$coef

    #If the node does not have parents
    if(length(parents)==0)
    {
      joint_mean_p[linear_names[i]]<-as.numeric(coeffs_p)
      # joint_mean_q[linear_names[i]]<-as.numeric(coeffs_q)
    }else{
      #Get the mean of the parents an multiply by the coefficients. Coeffs must be ordered according to the topological order
      idx_parents <- c()
      for(j in 1:length(parents))
      {
        idx_parents <- c(idx_parents,grep(paste0("^",parents[j],"$"),node_name_pos))
      }
      joint_mean_p[linear_names[i]] <- coeffs_p %*% c(1,joint_mean_p[idx_parents])
      # joint_mean_q[linear_names[i]] <- coeffs_q %*% c(1,joint_mean_q[idx_parents])
    }
  }
  # return(list(mean_p=joint_mean_p, mean_q=joint_mean_q))
  return(joint_mean_p)
}

#Compute Covariance matrix. variance is Koller implementation
compute_covariance <- function(structure, vM_p, gauss_p, is_dir)
{
  col_names <- node.ordering(structure)
  dir_names <- names(vM_p)
  dir_names <- col_names[col_names %in% dir_names]
  dir_names_trig <- c(paste0(dir_names,"-cos"),paste0(dir_names,"-sin"))
  linear_names <- names(gauss_p)
  linear_names <- col_names[col_names %in% linear_names]
  col_names_trig <- c(dir_names_trig,linear_names)
  node_name_pos <- gsub("-cos|-sin","",col_names_trig)

  # col_names <- col_names[col_names %in% names(gauss_p)]
  n <- length(col_names_trig)

  covariance_p <- matrix(0,ncol=n,nrow=n)
  colnames(covariance_p) <- col_names_trig
  rownames(covariance_p) <- col_names_trig

  for(i in 1:length(dir_names))
  {
    idx <- grep(dir_names[i],col_names_trig)
    for(j in idx)
    {
      covariance_p[j,j] <- 1/vM_p[[dir_names[i]]]$get_kappa()
    }
  }

  for(i in 1:length(linear_names))
  {
    variable <- linear_names[i]
    parents <- structure$nodes[[variable]]$parents
    coeffs_p <- gauss_p[[variable]]$coef
    variance_p <- gauss_p[[variable]]$sd^2

    #If the node does not have parents
    if(length(parents)==0)
    {
      covariance_p[linear_names[i],linear_names[i]] <- variance_p
    }else{
      #Get sigma
      idx_parents <- c()
      for(j in 1:length(parents))
      {
        idx_parents <- c(idx_parents,grep(paste0("^",parents[j],"$"),node_name_pos))
      }

      parent_coeffs_p <- coeffs_p[-1]
      covariance_p[linear_names[i],linear_names[i]]<- variance_p + t(parent_coeffs_p) %*% covariance_p[idx_parents,idx_parents] %*% parent_coeffs_p

      idx_node <- which(col_names_trig==variable)
      #Get covariates
      covariance_p[linear_names[i],1:(idx_node-1)] <- as.matrix(covariance_p[1:(idx_node-1),idx_parents])%*%as.matrix(parent_coeffs_p)
      covariance_p[1:(idx_node-1),linear_names[i]] <- covariance_p[linear_names[i],1:(idx_node-1)]
    }
  }
  return (covariance_p)
}

#Compute conditional parameters
compute_conditional_params<-function(mu, sigma, idx_evidence)
{
  #Variable definition
  sigma_yx <- as.matrix(sigma[-idx_evidence,idx_evidence])
  sigma_xx <- as.matrix(solve(sigma[idx_evidence,idx_evidence]))
  sigma_yy <- as.matrix(sigma[-idx_evidence,-idx_evidence])

  #Compute conditional mu and conditional covariance matrix
  beta_0 <- mu[-idx_evidence] + (sigma_yx %*% sigma_xx %*% mu[idx_evidence])
  beta <- sigma_yx %*% sigma_xx
  sigma_y <- sigma_yy - (sigma_yx %*% sigma_xx %*% t(sigma_yx))

  col_names <- colnames(sigma_yy)

  rownames(sigma_y)<-col_names
  colnames(sigma_y)<-col_names

  return(list(beta_0=beta_0, beta = beta, sigma_y = sigma_y))
}
