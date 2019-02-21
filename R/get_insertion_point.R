#' Get insertion point of somas
#'
#' Apply Gaussian mixture models to segment soma from dendrites
#'
#' @param directory path to the folder where PLY files of the somas are saved
#' @param output_shape_diameter path to the folder where somas with shape diameter will be saved
#' @param broken_mesh path to the folder where somas will be saved after that the GMM threshold removes the dendrites
#' @param output_poisson_reconstruction path to the folder where somas will be saved after the mesh closing process
#' @param final_result path to the folder where final segmented somas will be saved after the isolated pieces are removed
#' @return None
#'
#' @useDynLib SomaSimulation
#'
#' @examples
#' path <- "/home/universidad/RStudio Projects/SomaSimulation/data"
#' soma_parameters <- get_insertion_point(path_complete_neuron = file.path(path,"ao_results"), path_soma = file.path(path, "sdf_results"), path_off_export = file.path(tempdir(),"soma_off_files"), num_curves = 6)
#'
#' @export
get_insertion_point <- function(path_complete_neuron, path_soma, path_off_export, num_curves) {
  PLY_ao_files <- list.files(path_complete_neuron)[which(file_ext(list.files(path_complete_neuron)) == "ply")]
  PLY_sdf_files <- list.files(path_soma)[which(file_ext(list.files(path_soma)) == "ply")]

  PLY_files <- PLY_ao_files[PLY_ao_files %in% PLY_sdf_files]

  soma_parameters <- list()
  list_soma_names <- list()
  count <- 1
  for (file in PLY_files)
  {
    soma_name <- file_path_sans_ext(file)
    list_soma_names[[count]] <-soma_name

    print(paste0("Processing soma:", file))
    print(paste0("Soma number ",count, " of ", length(PLY_files)))
    ply_ao <- vcgPlyRead(file.path(path_complete_neuron,file))
    ply_sdf <- vcgPlyRead(file.path(path_soma,file))

    soma <- vcgIsolated(ply_sdf)
    soma<-vcgClean(soma,sel=0:7)

    ply_ao$vertices <- data.table(t(ply_ao$vb[-4,]))
    soma$vertices <- t(soma$vb[-4,])
    bounding_box_sdf <- apply(soma$vertices,2,function(x){return(c(min(x),max(x)))})

    apical_dendrite <- ply_ao$vertices[bounding_box_sdf[1,1]<V1 & bounding_box_sdf[2,1]>V1 & bounding_box_sdf[1,3]<V3 & bounding_box_sdf[2,3]>V3 & bounding_box_sdf[1,2]>V2]
    closest_points <- vcgClostKD(as.matrix(apical_dendrite), soma)
    closest_points$vertices <- t(closest_points$vb[-4,])
    mean_closest_point <- colMeans(closest_points$vertices)

    insertion_point <- which.min(rowSums((sweep(soma$vertices,2,mean_closest_point))^2))
    vcgOffWrite(soma, file.path(path_off_export,"temp.off"))
    soma$distances <- geodesic_distance(file.path(path_off_export,"temp.off"), insertion_point-1)

    factorized_model <- factorize_model(soma, num_curves)
    #Corregir direccion de las normales
    soma<-vcgClean(factorized_model$soma, sel=num_curves)
    soma$vertices <- factorized_model$soma$vertices
    soma$distances <- factorized_model$soma$distances

    curve <- list()
    for(i in 1:length(factorized_model$limits))
    {
      if(i==length(factorized_model$limits))
      {
        curve[[i]] = soma$vertices[which(soma$distances == max(soma$distances)), ]
      }else{
        curve[[i]] = soma$vertices[which(soma$distances == factorized_model$limits[i]), ]
      }
    }

    soma_parameters[[count]] <- compute_curve_parameters(soma, insertion_point, curve)
    count <- count +1
  }

  soma_parameters <- matrix(unlist(soma_parameters), ncol = length(soma_parameters[[1]]), nrow = length(soma_parameters), byrow = T)
  header <- c(paste0("h_",1:num_curves), paste0("theta_",2:num_curves), paste0("cos_phi_",2:num_curves), paste0("B_r_",1:(num_curves-1)),paste0("B_R_",1:(num_curves-1)),
    paste0("inst_Phi_",1:(num_curves-1)),paste0("inst_Theta_",1:(num_curves-1)))

  soma_parameters<-as.data.frame(soma_parameters)
  colnames(soma_parameters) <- header
  row.names(soma_parameters) <- unlist(list_soma_names)
  return(soma_parameters)
}
