#' Compute KL divergence between two von Mises distributions
#'
#' @examples
#'vM_p_var <- model$vMParams[[1]]$theta_2
#'vM_q_var <- model$vMParams[[2]]$theta_2
vM_KL <- function(vM_p_var,vM_q_var)
{
  k_p <- vM_p_var$get_kappa()
  k_q <- vM_q_var$get_kappa()
  mu_q <- vM_q_var$get_mean()-vM_p_var$get_mean()
  return(log(besselI(k_q,0))-log(besselI(k_p,0))+A1(k_p)*(k_p-k_q*cos(mu_q)))
}

#' Compute KL divergence of all directional variables
#'
#' @examples
#' vM_p <- model
#' vM_q <- model
#' cluster_p <- 1
#' cluster_q <- 2
#' compute_directional_KL(vM_p,vM_q)
compute_directional_KL<-function(model_p, model_q, cluster_p, cluster_q)
{
  total_KL <- 0
  vM_p <- model_p$vMParams[[cluster_p]]
  vM_q <- model_q$vMParams[[cluster_q]]
  for(i in 1:length(vM_p))
  {
    total_KL <- total_KL + vM_KL(vM_p[[i]],vM_q[[i]])
  }
  return(total_KL)
}

#' Compute KL divergence of linear variables
#'
#' @examples
#' cluster_p <- 2
#' cluster_q <- 2
#' is_dir <- dataset$is_dir
#' compute_linear_KL(model, clustering_model, cluster_p, cluster_q, is_dir)
compute_linear_KL<-function(model_p, model_q, cluster_p, cluster_q, is_dir)
{
  structure_p <- model_p$structure
  structure_q <- model_q$structure
  gauss_p <- model_p$gaussParams[[cluster_p]]
  gauss_q <- model_q$gaussParams[[cluster_q]]

  vM_p <- model_p$vMParams[[cluster_p]]
  vM_q <- model_q$vMParams[[cluster_q]]

  v_kappa <- matrix(0,nrow=1,ncol= length(vM_p))
  for(i in 1:length(vM_p))
  {
    v_kappa[i] <- vM_p[[i]]$get_kappa()
  }
  names(v_kappa) <- names(vM_p)
  v_ratioI <- c(besselI(v_kappa,2)/besselI(v_kappa,0))
  v_A1 <- c(A1(v_kappa))

  #Compute the parameters of the multivariate conditional gaussian conditioned to the directional variables
  params_p <- obtain_cond_params(structure_p, vM_p, gauss_p, is_dir)
  params_q <- obtain_cond_params(structure_q, vM_q, gauss_q, is_dir)
  s_q <- params_q$sigma_y
  s_p <- params_p$sigma_y[rownames(s_q),colnames(s_q)]

  inv_s_q <- solve(s_q)#Inverse of the covariance matrix of Q

  gauss_t <- list()
  gauss_t$beta_0 <- params_q$beta_0 - params_p$beta_0[rownames(params_q$beta_0),]
  beta_1_idx <- grep("cos",colnames(params_q$beta))
  beta_2_idx <- grep("sin",colnames(params_q$beta))
  gauss_t$beta_1 <- params_q$beta[,beta_1_idx] - params_p$beta[rownames(params_q$beta),beta_1_idx]
  gauss_t$beta_2 <- params_q$beta[,beta_2_idx] - params_p$beta[rownames(params_q$beta),beta_2_idx]
  #The constant
  C <- sum(diag(inv_s_q%*%s_p))-length(gauss_p)+log(abs(det(s_q)))-log(abs(det(s_p)))
  for(i in 1:length(gauss_p))
  {
    for(j in 1:length(gauss_p))
    {
      c1 <- gauss_t$beta_0[i]*gauss_t$beta_0[j]
      c2 <- 2*gauss_t$beta_0[i]*(gauss_t$beta_1[j,]%*%v_A1)
      c3 <- sum(gauss_t$beta_1[i,]*gauss_t$beta_1[j,]/2)+(gauss_t$beta_1[i,]*gauss_t$beta_1[j,]/2)%*%v_ratioI
      temp_c4 <- (gauss_t$beta_1[i,]*v_A1)%*%t(gauss_t$beta_1[j,]*v_A1) #Generate all the combinations of the multiplications
      diag(temp_c4) <- 0 #Remove elements in the diagonal because they are c3
      c4 <- sum(temp_c4)
      c5 <- sum(gauss_t$beta_2[i,]*gauss_t$beta_2[j,]/2)-(gauss_t$beta_2[i,]*gauss_t$beta_2[j,]/2)%*%v_ratioI
      C <- C + inv_s_q[i,j]*(c1 + c2 + c3 + c4 + c5)
    }
  }
  return(C/2)
}

discrete_KL_divergence <- function(prior_p, prior_q)
{
  return(sum(prior_p*(log(prior_p)-log(prior_q))))
}


obtain_cond_params <- function(structure, vM_p, gauss_p, is_dir)
{
  mean <- compute_joint_mean(structure, vM_p, gauss_p, is_dir)
  covariance <- compute_covariance(structure, vM_p, gauss_p, is_dir)

  col_names <- node.ordering(structure)
  dir_names <- names(vM_p)
  dir_names <- col_names[col_names %in% dir_names]
  dir_names_trig <- c(paste0(dir_names,"-cos"),paste0(dir_names,"-sin"))
  idx <- c()
  for(j in 1:length(dir_names))
  {
    idx <- c(idx,grep(paste0(dir_names[j],"-"),dir_names_trig))
  }
  params_p <- compute_conditional_params(mean, covariance, idx)
  return(params_p)
}


entropy_vM <- function(model, cluster)
{
  entropy <- 0
  for(i in 1:length(model$vMParams[[cluster]]))
  {
    kappa <- model$vMParams[[cluster]][[i]]$get_kappa()
    entropy <- entropy + log(2*pi*besselI(kappa,0))-kappa*besselI(kappa,1)/besselI(kappa,0)
  }
  return(entropy)
}

#' @examples
#' entropy <- entropy_gaussian(models$original_model,1)
#' cluster_p <- 1
#' cluster_q <- 3
#' KL <- compute_linear_KL(models$original_model, models$SEM_model, cluster_p, cluster_q, is_dir)
#' KL / entropy
entropy_gaussian <-function(model, cluster)
{
  structure <- model$structure
  gauss_p <- model$gaussParams[[cluster]]
  gauss_q <- model$gaussParams[[1]]

  vM_p <- model$vMParams[[cluster]]
  vM_q <- model$vMParams[[1]]

  v_kappa <- matrix(0,nrow=1,ncol= length(vM_p))
  for(i in 1:length(vM_p))
  {
    v_kappa[i] <- vM_p[[i]]$get_kappa()
  }
  names(v_kappa) <- names(vM_p)
  v_ratioI <- c(besselI(v_kappa,2)/besselI(v_kappa,0))
  v_A1 <- c(A1(v_kappa))

  #Compute the parameters of the multivariate conditional gaussian conditioned to the directional variables
  params_p <- obtain_cond_params(structure, vM_p, gauss_p, is_dir)
  s_p <- params_p$sigma_y
  k <- length(model$gaussParams[[1]])
  entropy <- k/2+k/2*log(2*pi)+0.5*log(det(s_p))
  return(entropy)
}

compute_entropy <- function(model, cluster)
{
  entropy <- entropy_vM(model, cluster)
  entropy <- entropy + entropy_gaussian(model, cluster)
  return(entropy)
}

clustering_correspondence <- function(model_p, model_q)
{
  weights <- matrix(0,ncol=ncol(model_q$weights),nrow=ncol(model_q$weights))

  for(p in 1:ncol(model_q$weights))
  {
    for(q in 1:ncol(model_q$weights)){
      cluster_p <- p
      cluster_q <- q
      weights[p,q] <- compute_directional_KL(model_p, model_q, cluster_p, cluster_q)
      weights[p,q] <- weights[p,q] + compute_linear_KL(model_p, model_q, cluster_p, cluster_q, is_dir)
    }
  }
  correspondence <- solve_LSAP(weights)
  return(list(weights = weights, correspondence = correspondence, similarity = diag(weights[1:length(model_p$gaussParams),c(correspondence)])))
}

#' @examples
#' s<-clustering_correspondence(models$original_model,models$SEM_model)
#' KL1 <- model_KL_divergence(models$original_model,models$SEM_model, s$correspondence)
#' s<-clustering_correspondence(models$original_model,models$simple_SEM_model)
#' KL2 <- model_KL_divergence(models$original_model,models$simple_SEM_model, s$correspondence)
#' entropy <- entropy_model(models$original_model)
model_KL_divergence <- function(model_p, model_q, correspondence){
  prior_p <- rep(1,length(model_p$gaussParams))/length(model_p$gaussParams)
  prior_q <- colSums(model_q$weights)/sum(model_q$weights)
  total_KL <- discrete_KL_divergence(prior_p,prior_q)
  for(p in 1:length(correspondence))
  {
    cluster_p <- p
    cluster_q <- c(correspondence)[p]
    total_KL <- total_KL + prior_p[p] * compute_directional_KL(model_p, model_q, cluster_p, cluster_q)
    total_KL <- total_KL + prior_p[p] * compute_linear_KL(model_p, model_q, cluster_p, cluster_q, is_dir)
  }
  return(total_KL)
}

entropy_model <- function(model)
{

  if(is.null(model$weights))
  {
    prior <- rep(1,length(model$gaussParams))/length(model$gaussParams)
  }else{
    prior <- colSums(model$weights)/sum(model$weights)
  }

  entropy <- - sum(prior*log(prior))
  for(cluster in 1:length(prior))
  {
    entropy <- entropy + prior[cluster] * entropy_vM(model, cluster)
    entropy <- entropy + prior[cluster] * entropy_gaussian(model, cluster)
  }

  return(entropy)
}


# KL <- clusters_KL_divergence(model, cluster_p = 1, cluster_q = 2, is_dir = dataset$is_dir)
clusters_KL_divergence <- function(model, cluster_p, cluster_q, is_dir){
  total_KL <- 0
  total_KL <- total_KL +  compute_directional_KL(model, model, cluster_p, cluster_q)
  total_KL <- total_KL +  compute_linear_KL(model, model, cluster_p, cluster_q, is_dir)
  return(total_KL)
}
