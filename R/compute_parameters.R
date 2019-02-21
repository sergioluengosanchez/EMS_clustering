check_orientation <- function(soma_model, curve)
{
  curve_centroids <- zeros(1,3)
  z_axis <- c(0,0,1)
  #Compute parameters of each curve.
  PCA <- prcomp(curve[[2]])
  perp_PCA_matrix <- t(PCA$rotation[,3])

  if(isTRUE(all.equal(perp_PCA_matrix, z_axis)))
  {
    rotated_curve <- curve[[2]][,1:2]
  }else{
    u <- vcrossp(perp_PCA_matrix, z_axis)
    u <- u / c(normv(u))
    theta <- acos(perp_PCA_matrix %*% z_axis)
    rotationMatrix <- rotateMatrix(u, theta);

    rotated_curve <- t(rotationMatrix %*% t(curve[[2]]))

    #Check if the elipse is over the XY plane. It is needed to avoid
    #rotation in the wrong direction
    if(mean(rotated_curve[ ,3]) < 0)
    {
      perp_PCA_matrix < --perp_PCA_matrix;

      u <- vcrossp(perp_PCA_matrix, z_axis);
      u <- u / c(normv(u))
      theta <- acos(perp_PCA_matrix %*% z_axis)
      rotationMatrix <- rotateMatrix(u, theta)

      rotated_curve <- t(rotationMatrix %*% t(curve[[2]]))
    }
    ellipse <- fit.ellipse(rotated_curve[,1],rotated_curve[,2])
    curve_centroids <- t(rotationMatrix) %*% c(ellipse$center, mean(rotated_curve[,3]))
  }

  #Rotate the soma to align the vector between de insertion point and
  #the centroid of the first curve with the Z axis. Also, rotate the
  #curves and the centroids.
  #u=rotation axis, theta=rotation angle
  unit_vector <- c(curve_centroids) / c(normv(curve_centroids))
  u <- vcrossp(unit_vector, z_axis);
  u <- u / c(normv(u))
  theta <- acos(unit_vector %*% z_axis)

  rotationMatrix <- rotateMatrix(u, theta)

  farest_point_rotation <- (rotationMatrix %*% curve[[length(curve)]])

  if(farest_point_rotation[3] < 0){
    insertion_point <- 2 * t(curve_centroids)
  }else{
    insertion_point <- c(0, 0, 0)
  }

  return(insertion_point)
}



compute_curve_parameters <- function(soma_model, insertion_point_idx, curve)
{
  z_axis <- c(0, 0, 1)
  #Soma is translated to place the insertion point as the origin
  insertion_point <- soma_model$vertices[insertion_point_idx,]
  soma_model$vertices <- sweep(soma_model$vertices, 2, insertion_point, "-")
  for (i in 1:size(curve,2))
  {
    if(length(curve[[i]]) == 3)
    {
      curve[[i]] <- curve[[i]] - insertion_point
    }else{
      curve[[i]] <- sweep(curve[[i]], 2, insertion_point, "-")
    }

  }

  insertion_point <- check_orientation(soma_model,curve)
  if(sum(insertion_point) != 0)
  {
    soma_model$vertices <- sweep(soma_model$vertices, 2, insertion_point, "-")
    for (i in 1:length(curve))
    {
      if(length(curve[[i]]) == 3)
      {
        curve[[i]] <- curve[[i]] - insertion_point
      }else{
        curve[[i]] <- sweep(curve[[i]], 2, insertion_point, "-")
      }
    }

    insertion_point <- c(0, 0, 0)
  }

  #Compute the parameters of the curve, i.e., centroids, semiaxes, angle
  #of the curve, perpendicular vector to the curve and coefficients for
  #the z axis.
  perp_vectors_sph <- matrix(0, length(curve) - 2, 2)   #Perpendicular vector to curve in spherical coords.
  ellipse_angles <- matrix(0, length(curve) - 2, 1)     #Angle of the elipse in the XY plane.
  major_axis <- matrix(0, length(curve) - 2, 1)     #Relation between semiaxes.
  curve_centroids <- matrix(0, length(curve), 3)       #Centroids of the curves.
  perp_PCA_matrix <- matrix(0, length(curve) - 2, 3)   #Perpendicular vector to curve in cartesian coords.
  ellipse_angles_vectors <- matrix(0, length(curve) - 2, 3)  #Vector representation of the angle of the ellipse.
  minor_axis <- matrix(0, length(curve) -2, 1)            #One of the minor_axis of the ellipse.


  #Compute parameters of each curve.
  for (i in 2:(length(curve)-1))
  {
    PCA <- prcomp(curve[[i]])$rotation
    perp_PCA_matrix[i-1, ] <- PCA[,3]
    if(isTRUE(all.equal(perp_PCA_matrix[i-1,], z_axis)))
    {
      rotated_curve <- curve[[i]]
      ellipse <- fit.ellipse(rotated_curve[,1],rotated_curve[,2])
      ellipse_angles_vectors[i-1,] <- c(cos(ellipse$angle), sin(ellipse$angle), 0)
      curve_centroids[i,] <- c(ellipse$center, mean(rotated_curve[,3]))
    }else{
      u <- vcrossp(perp_PCA_matrix[i-1,], z_axis)
      u <- u / c(normv(u))
      theta <- acos(perp_PCA_matrix[i-1,] %*% z_axis)
      rotationMatrix <- rotateMatrix(u, theta)

      rotated_curve <- rotated_curve <- t(rotationMatrix %*% t(curve[[i]]))
      #Aseguramos que la elipse se encuentra por encima del plano XY.
      #Es necesario para rotar siempre en la misma direccion, si no rotariamos al reves en algunos casos.
      if(mean(rotated_curve[,3]) < 0)
      {
        perp_PCA_matrix[i-1, ]  <- -perp_PCA_matrix[i-1, ]

        u <- vcrossp(perp_PCA_matrix[i-1,], z_axis)
        u <- u / c(normv(u))
        theta <- acos(perp_PCA_matrix[i-1,] %*% z_axis)
        rotationMatrix <- rotateMatrix(u,theta);

        rotated_curve <- rotated_curve <- t(rotationMatrix %*% t(curve[[i]]))
      }

      ellipse <- fit.ellipse(rotated_curve[,1],rotated_curve[,2])

      ellipse_angles_vectors[i-1,] <- t(rotationMatrix) %*% c(cos(ellipse$angle), sin(ellipse$angle), 0)
      curve_centroids[i,] <- t(rotationMatrix) %*% c(ellipse$center, mean(rotated_curve[,3]))
    }
    major_axis[i-1] <- ellipse$major
    minor_axis[i-1] <- ellipse$minor
    ellipse_angles[i-1] <- ellipse$angle
  }

  curve_centroids[nrow(curve_centroids),] <- curve[[length(curve)]]

  # Rotate the spine to align the vector between de inserction point and
  # the centroid of the first curve with the Z axis. Also, rotate the
  # curves and the centroids.
  # u=rotation axis, theta=rotation angle
  unit_vector <- curve_centroids[2,] / c(normv(curve_centroids[2,]))
  u <- vcrossp(unit_vector, z_axis)
  u <- u / c(normv(u))
  theta <- acos(unit_vector %*% z_axis)

  rotationMatrix <- rotateMatrix(u, theta)
  soma_model$vertices <- t(rotationMatrix %*% t(soma_model$vertices))

  for (i in 1:length(curve))
  {
    if(length(curve[[i]]) == 3)
    {
      curve[[i]] <- t(rotationMatrix %*% c(curve[[i]]))
    }else{
      curve[[i]] <- t(rotationMatrix %*% t(curve[[i]]))
    }
  }

  curve_centroids <- t(rotationMatrix %*% t(curve_centroids))

  perp_PCA_matrix <- t(rotationMatrix %*% t(perp_PCA_matrix));
  ellipse_angles_vectors <- t(rotationMatrix %*% t(ellipse_angles_vectors))

  #Rotate the spine to place the farest point of the spine parallel to the
  #X axis. It is rotated around the Z axis.
  unit_vector <- c(curve[[length(curve)]][1], curve[[length(curve)]][2], 0) / c(normv(c(curve[[length(curve)]][1], curve[[length(curve)]][2], 0)))

  u <- vcrossp(unit_vector, c(1,0,0));
  u <- u / c(normv(u))
  theta <- acos(unit_vector %*% c(1,0,0))

  rotationZ <- rotateMatrix(u, theta)

  soma_model$vertices <- t(rotationZ %*% t(soma_model$vertices))

  for (i in 1:length(curve))
  {
    if(length(curve[[i]]) == 3)
    {
      curve[[i]] <- t(rotationZ %*% c(curve[[i]]))
    }else{
      curve[[i]] <- t(rotationZ %*% t(curve[[i]]))
    }
  }

  curve_centroids <- t(rotationZ %*% t(curve_centroids))
  perp_PCA_matrix <- t(rotationZ %*% t(perp_PCA_matrix))
  ellipse_angles_vectors <- t(rotationZ %*% t(ellipse_angles_vectors))

  for (i in 1:nrow(perp_PCA_matrix))
  {
    u <- vcrossp(perp_PCA_matrix[i,], z_axis);
    u <- u / c(normv(u))
    theta <- acos(perp_PCA_matrix[i,] %*% z_axis)
    rotationMatrix <- rotateMatrix(u, theta)
    rotated_angle_matrix <- t(rotationMatrix %*% ellipse_angles_vectors[i,])
    ellipse_angles[i] <- atan2(rotated_angle_matrix[2], rotated_angle_matrix[1])
  }


  spherical_coords <- as.data.frame(cartesian2Spherical(perp_PCA_matrix[,1], perp_PCA_matrix[,2], perp_PCA_matrix[,3]))
  perp_vectors_sph[,1] <- spherical_coords$azimuth
  perp_vectors_sph[,2] <- spherical_coords$elevation

  #Compute the parameters needed to build the skeleton of the spine.
  #|h|, phi, theta, alpha
  vector_length <- matrix(0, nrow(curve_centroids) - 1, 1)
  phi <- matrix(0, nrow(curve_centroids) - 2, 1)
  cosAngle <- matrix(0, nrow(curve_centroids) - 2, 1)

  for (i in 1:(nrow(curve_centroids)-1))
  {
    section_vector <- curve_centroids[i+1,] - curve_centroids[i,]
    vector_length[i] <- c(normv(section_vector))
  }


  for (i in 2:(nrow(curve_centroids)-1))
  {
    translation <- sweep(curve_centroids[(i-1):(i+1),], 2, curve_centroids[i,])
    translation[1,] <- translation[1,] / c(vector_length[i-1])
    translation[3,] <- translation[3,] / c(vector_length[i])

    u <- vcrossp(translation[1,], c(1,0,0))
    u <- u / c(normv(u))
    theta <- acos(translation[1,] %*% c(1,0,0))
    rotationMatrix <- rotateMatrix(u, theta)

    transRotated <- t(rotationMatrix %*% t(translation))
    cosAngle[i-1] <- transRotated[1,] %*% transRotated[3,]
    xyPosition <- c(cosAngle[i-1], sqrt(1-cosAngle[i-1]^2), 0)

    r <- xyPosition[2:3] / c(normv(xyPosition[2:3]))
    s <- transRotated[3,2:3] / c(normv(transRotated[3,2:3]))
    phi[i-1]  <- atan2(s[2], r %*% s)
  }

  ellipse_area <- pi * minor_axis * major_axis #Compute the area of the ellipses
  #Compute how many times is ellipse j bigger or smaller than ellipse i.
  #Results are saved in the lower triangular matrix
  proportional_areas <- ellipse_area %*% t((c(1)/ellipse_area))
  proportional_areas[upper.tri(proportional_areas)] <- 0
  proportional_areas <- proportional_areas * (matrix(1, length(ellipse_area), length(ellipse_area)) - diag(length(ellipse_area)))
  proportional_areas_vector <- proportional_areas[which(proportional_areas != 0)]

  instance <- c(vector_length,phi,cosAngle,minor_axis,major_axis,perp_vectors_sph[,1],perp_vectors_sph[,2])

  return(instance)
}

