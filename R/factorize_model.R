#Soma is customize mesh that has the propery distance which is the geodesic distance from the insertion vertex to the remaining vertices
#num_regions are the number of regions the mesh should be divided
factorize_model<-function(soma, num_regions)
{
  minval <- min(soma$distances)
  maxval <- max(soma$distances)
  range <- maxval - minval
  limits <- seq(from = minval, to = maxval, by = range/num_regions)

  # Divide
  result<-list()
  result$newModel <- soma
  result$newFuncVals <- soma$distances
  result$newModel$remvert <- NULL
  result$newModel$normals <- NULL
  for (f in 2:num_regions)
  {
    limitValue <- limits[f] # Get limit
    result <- limitModel(result$newModel, result$newFuncVals, limitValue)
  }
  factorized_soma <- result$newModel
  factorized_soma$distances <- result$newFuncVals
  return(list(soma=factorized_soma, limits=limits))
}



limitModel <- function(model, funcVals, limitValue)
{
  faces <- t(model$it)
  vertices <- model$vertices
  nV = nrow(vertices)

  cond <- funcVals >= limitValue;
  factors <- matrix(0, ncol = 1, nrow = nV)
  factors[cond] <- 1

  #Check all triangles
  newFuncVals = funcVals
  newvertices = vertices
  newfaces = faces

  #Get faces thar are cut in different regions
  cut_faces <- (factors[faces[,1]] != factors[faces[,2]]) | (factors[faces[,2]] != factors[faces[,3]]) | (factors[faces[,1]] != factors[faces[,3]])
  erasefaces <- which(cut_faces)

  vsoloIndex <- matrix(0, ncol = 1, nrow = length(erasefaces))
  vpairIndex <- matrix(0, ncol = 2, nrow = length(erasefaces))

  for(i in 1:length(erasefaces))
  {
    if (factors[faces[erasefaces[i],2]] == factors[faces[erasefaces[i],3]])
    {
      vsoloIndex[i] <- faces[erasefaces[i],1]
      vpairIndex[i,] = faces[erasefaces[i],2:3]
    }else if(factors[faces[erasefaces[i],1]] == factors[faces[erasefaces[i],3]])
    {
      vsoloIndex[i] <- faces[erasefaces[i],2]
      vpairIndex[i,] = faces[erasefaces[i],c(1,3)]
    }else{
      vsoloIndex[i] <- faces[erasefaces[i],3]
      vpairIndex[i,] = faces[erasefaces[i],1:2]
    }
  }

  vsolo <- newvertices[vsoloIndex,]
  vpair <- newvertices[t(vpairIndex),]


  v1 <- vpair[seq(from = 1, to = nrow(vpair), by=2),] - vsolo
  length_new_point <- abs(limitValue - newFuncVals[vsoloIndex]) / sqrt(rowSums(v1^2))
  vnew1 = (length_new_point * v1) + vsolo

  vnew1Index = nrow(newvertices) + 1:nrow(vnew1);
  newFuncVals<- c(newFuncVals, rep(limitValue, nrow(vnew1)))
  newvertices <- rbind(newvertices, vnew1)
  factors <- matrix(c(factors, rep(0, nrow(vnew1))), ncol=1)

  v2 <- vpair[seq(from = 2, to = nrow(vpair), by=2),] - vsolo
  length_new_point <- abs(limitValue - newFuncVals[vsoloIndex]) / sqrt(rowSums(v2^2))
  vnew2 = (length_new_point * v2) + vsolo

  vnew2Index = nrow(newvertices) + 1:nrow(vnew2);
  newFuncVals<- c(newFuncVals, rep(limitValue, nrow(vnew2)))
  newvertices <- rbind(newvertices, vnew2)
  factors <- matrix(c(factors, rep(0, nrow(vnew2))), ncol=1)

  newfaces <- rbind(newfaces, cbind(vsoloIndex,vnew1Index,vnew2Index))
  newfaces <- rbind(newfaces, cbind(vnew1Index,vpairIndex[,1], vpairIndex[,2]))
  newfaces <- rbind(newfaces, cbind(vnew2Index,vnew1Index, vpairIndex[,2]))

  cond <- rep(T, nrow(newfaces))
  cond[erasefaces] <- F
  newfaces <- newfaces[cond,]

  result <- list()
  result$newModel <- model
  result$newModel$vb <- t(cbind(newvertices,1))
  result$newModel$it <- t(newfaces)
  result$newModel$vertices <- newvertices
  class(result$newModel) <- "mesh3d"
  result$newFuncVals <- newFuncVals

  return(result)
}
