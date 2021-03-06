sampleSomas<-function(model, is_dir, cluster_id, num_simulations)
{
  node_order <- node.ordering(model$structure)

  vMParams <- model$vMParams[[cluster_id]]
  gaussParams <- model$gaussParams[[cluster_id]]

  var_names <- names(model$structure$nodes)
  linear_names <- names(gaussParams)

  simulations <- matrix(0, nrow=num_simulations, ncol=length(var_names))
  colnames(simulations) <- var_names
  for(var in node_order)#Simulate each variable according to the topological order
  {
    if(var %in% linear_names)#Check if the variable is linear or directional
    {
      parents <- model$structure$nodes[[var]]$parents
      if(length(parents) == 0)#If the variable does not have parents
      {
        if(grepl("cos\\_phi\\_*", var))
        {
          simulations[,var] <- as.numeric(rtruncnorm(num_simulations, gaussParams[[var]]$coef, gaussParams[[var]]$sd,a=-1, b=-0.7))
        }else{
          if(grepl("h\\_*", var) | grepl("B\\_r\\_*", var) | grepl("B\\_R\\_*", var))
          {
            simulations[,var] <- as.numeric(rtruncnorm(num_simulations, gaussParams[[var]]$coef, gaussParams[[var]]$sd,a=0, b=Inf))
          }else{
          simulations[,var] <- as.numeric(rnorm(num_simulations, gaussParams[[var]]$coef, gaussParams[[var]]$sd))
          }
        }

      }else{
        linear_dataset <- get_linear_dataset(simulations, is_dir)
        node_data <- get_node_data(simulations, linear_dataset, is_dir, which(var_names %in% parents), which(var == var_names))
        parent_data <- cbind(1, node_data[,-ncol(node_data)])
        if(grepl("cos\\_phi\\_*", var))
        {
          simulations[,var] <- as.numeric(rtruncnorm(num_simulations, parent_data %*% gaussParams[[var]]$coef, gaussParams[[var]]$sd,a=-1, b=-0.8))
        }else{
          if(grepl("^h\\_*", var) | grepl("B\\_r\\_*", var) | grepl("B\\_R\\_*", var))
          {
            simulations[,var] <- as.numeric(rtruncnorm(num_simulations, parent_data %*% gaussParams[[var]]$coef, gaussParams[[var]]$sd,a=0, b=Inf))
          }else{
            simulations[,var] <- as.numeric(rnorm(num_simulations, parent_data %*% gaussParams[[var]]$coef, gaussParams[[var]]$sd))
          }
        }
      }
    }else{
      simulations[,var] <- as.numeric(rvonmises(num_simulations, vMParams[[var]]$get_mean(), vMParams[[var]]$get_kappa()))
    }
  }

  return(simulations)
}


#' 3D virtual spine
#'
#' Generate the 3D representation of the virtual spines
#'
#' @param simulated_somas a matrix where each row is a virtual spine
#' @param idx is number of the row number of the spine to representate
#'
#' @return rgl object 3D mesh3d
#'
#' @examples
#' simulated_somas <- sampleSomas(bestModel, dataset$is_dir, cluster_id = 1, num_simulations = 20)
#' mesh <- simulation3Dmesh(simulated_somas,1,iterations=0)
#' shade3d(mesh,col="red")
#' mesh <- simulation3Dmesh(simulated_somas,2,iterations=0)
#' shade3d(mesh,col="red")
#' mesh <- simulation3Dmesh(simulated_somas,7,iterations=0)
#' shade3d(mesh,col="red")
#'
#'
#'
#' @export
simulation3Dmesh<-function(simulated_somas,idx,iterations = 2)
{
  ellipsesSpine <- simulation3Dellipses(simulated_somas,idx)
  spine <- simulation3Dfaces(ellipsesSpine$ellipses,ellipsesSpine$skeleton)

  mesh <- list()
  mesh$vb <- rbind(t(spine$vertices),1)
  mesh$it <- t(spine$faces)
  class(mesh) <- "mesh3d"
  # shade3d(mesh, col="red")
  subdivision<-vcgSubdivide(mesh,type="Loop",looptype="continuity",silent=TRUE,iterations = iterations)
  return(subdivision)
}

#' Generate the 3D representation of the spine represented uniquely by the ellipses
#'
#' Generate the 3D representation of the spine represented uniquely by the ellipses
#'
#' @param simulated_somas a matrix where each entry represents the values of a simulated spines
#' @param idx an integer that determines which row of simulated_somas must me selected to simulate that spine
#'
#' @return a list of three elements, first one is a matrix with all the vertices of all the ellipses,
#' second is a list saving in each element a matrix with the vertices of a ellipse, third are the centroids of the ellipses.
#'
#' @noRd
simulation3Dellipses<-function(simulated_somas,idx)
{
  curve_points<-data.frame(x = 0, y = 0, z = 0)
  cartesianEllipses<-list()
  cartesianEllipses[[1]]<-c(0,0,0)

  #Column index of the variables
  height_position <- grep("^h\\_*", colnames(simulated_somas))
  theta_position <- grep("theta\\_*", colnames(simulated_somas))
  cos_phi_position <- grep("cos\\_phi\\_*", colnames(simulated_somas))
  B_r_position <- grep("B\\_r\\_*", colnames(simulated_somas))
  B_R_position <- grep("B\\_R\\_*",colnames(simulated_somas))
  PCA_azimuth_position <- grep("inst\\_Phi\\_*", colnames(simulated_somas))
  PCA_theta_position <- grep("inst\\_Theta\\_", colnames(simulated_somas))

  vector_length <- as.numeric(simulated_somas[idx, height_position])
  phi <- as.numeric(simulated_somas[idx, theta_position])
  cosAngle <- as.numeric(simulated_somas[idx, cos_phi_position])
  B_r <- as.numeric(simulated_somas[idx, B_r_position])
  B_R <- as.numeric(simulated_somas[idx, B_R_position])
  perp_vectors_sph <- matrix(as.numeric(simulated_somas[idx, c(PCA_azimuth_position, PCA_theta_position)]), ncol = 2)

  skeleton <- matrix(0, nrow = length(height_position), ncol = 3)
  skeleton[1,]=c(0,0,1)*vector_length[1];

  for(i in 1:length(cosAngle))
  {
    if(i==1)
    {
      x<-c(0,0,-1);
    }else{
      x<-skeleton[i-1,]-skeleton[i,];
      x<-x/c(normv(x));
    }

    y<-c(0,0,0);

    u <- vcrossp(x,c(1,0,0));
    u <- as.numeric(u) / c(normv(u))
    theta <- acos(x %*% c(1,0,0));
    rotationMatrix<-rotateMatrix(u,theta);

    xyPosition<-c(cosAngle[i],sqrt(1-cosAngle[i]^2),0);
    rotationMatrix2<-rotateMatrix(c(1,0,0),phi[i]);
    original<-t(rotationMatrix2%*%xyPosition);
    lastArray<-as.numeric(t(t(rotationMatrix)%*%t(original))) * vector_length[i+1];
    skeleton[i+1,]<-lastArray+skeleton[i,];
  }

  major_axis <- B_R
  minor_axis <- B_r

  perp_PCA_matrix <- spherical2Cartesian(perp_vectors_sph[ ,1], perp_vectors_sph[ ,2], 1)
  for (i in 1:nrow(perp_PCA_matrix))
  {
    u <- vcrossp(perp_PCA_matrix[i,],c(0,0,1));
    u <- as.numeric(u) / c(normv(u))
    rotation_angle <- acos(perp_PCA_matrix[i,] %*% c(0,0,1));
    rotationMatrix<-rotateMatrix(u,rotation_angle);

    centers<-as.numeric(rotationMatrix%*%skeleton[i,]);

    fit <- list()
    fit$angle <- 1
    fit$maj <- max(major_axis[i],minor_axis[i])
    fit$min <- min(major_axis[i],minor_axis[i])
    fit$center <- centers[1:2]

    points_ellipse <- get.ellipse(fit)
    x <- points_ellipse[,1]
    y <- points_ellipse[,2]
    z <- centers[3]

    u = vcrossp(c(cos(1), sin(1), 0),c(1,0,0));
    u = as.numeric(u) / c(normv(u))
    rotationZ<-rotateMatrix(u,1);

    rotated_ellipse<-matrix(c(x, y, matrix(0,nrow=length(x),ncol=1)),nrow=length(x),ncol=3);
    rotated_ellipse<-t(rotationZ%*%t(rotated_ellipse));

    rotated_ellipse[,3]<-centers[3];
    elipse<-t(t(rotationZ)%*%t(rotated_ellipse));

    recovered_ellipse <- t(t(rotationMatrix) %*% t(elipse));
    cartesianEllipses[[i+1]]<-recovered_ellipse
    curve_points <- rbind(curve_points, data.frame(x = recovered_ellipse[,1], y = recovered_ellipse[,2], z = recovered_ellipse[,3]))
  }

  curve_points <- rbind(curve_points,skeleton[nrow(skeleton),])
  return(list(ellipses=cartesianEllipses, vertices=curve_points,skeleton=skeleton))
}

#' Generate the surface between ellipses
#'
#' Generate the surface between each pair of consecutive ellipses
#'
#' @param ellipses is a list where each element is a matrix of dimension nx3 representing the X,Y,Z coordinates of the ellipse
#' @param skeleton is a matrix with the centroids of each ellipse
#'
#' @return a list of two elements where each element is a matrix, first element is the set of vertices defining the spine and the second the faces joining those vertices
#'
#' @noRd
simulation3Dfaces<-function(ellipses,skeleton)
{
  num_points<-nrow(ellipses[[2]])
  points<-rbind(ellipses[[1]],ellipses[[2]])
  faces<-cbind(1,seq(2,(nrow(points)-1),by=1),seq(3,nrow(points),by=1))
  faces<-rbind(faces,c(1,nrow(points),2))

  for(i in 2:(length(ellipses)-1))
  {
    first_point_1<-nrow(points)-num_points+1
    first_point_2<-nrow(points)+1

    idx_first_ellipse<-seq(first_point_1,first_point_1+num_points-1,by=1)
    idx_snd_ellipse<-seq(first_point_2,first_point_2+num_points-1,by=1)

    first_triang<-matrix(c(idx_first_ellipse,idx_snd_ellipse,idx_snd_ellipse+1),ncol=3,byrow=F)
    first_triang[nrow(first_triang),3]<-idx_snd_ellipse[1]

    snd_triang<-matrix(c(idx_snd_ellipse,idx_first_ellipse,idx_first_ellipse-1),ncol=3,byrow=F)
    snd_triang[1,3]<-idx_first_ellipse[length(idx_first_ellipse)]

    faces2<-rbind(first_triang,snd_triang)

    faces<-rbind(faces,first_triang,snd_triang)

    points<-rbind(points,ellipses[[i+1]])
  }

  first_point_1 <- nrow(points)-num_points+1
  first_point_2 <- nrow(points)+1

  temp_faces<-cbind(first_point_2,seq(first_point_1,nrow(points),by=1)+1,seq(first_point_1,nrow(points),by=1))
  temp_faces[nrow(temp_faces),2]<-first_point_1[1]

  faces<-rbind(faces,temp_faces)

  points<-rbind(points,skeleton[nrow(skeleton),])

  return(list(vertices = points, faces = faces))
}
