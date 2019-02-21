example<-function()
{
  path <- "data"

  #Compute the parameters of the somas. It can take some time, you can avoid this computation loading directly the data from the csv
  tdir<-tempdir()
  data <- get_insertion_point(path_complete_neuron = file.path(path,"ao_results"), path_soma = file.path(path, "sdf_results"), path_off_export = tdir, num_curves = 6)

  #Avoid preprocessing loading the data from the csv file
  data <- read.csv(file.path(path,"soma_dataset.csv"))

  #Generate a list of two elements, the first one is the data ordered and the second one is a binary vector denoting if the variables are directional
  dataset <- set_directional_variables(data = data)

  #Perform the clustering of according to the Extended Mardia-Sutton mixture model
  model <- EMS_clustering(data = dataset$data, is_dir = dataset$is_dir, num_clusters = 4, num_parents = 2, cores = 4)

  #To load the precomputed model used for the article load the following model
  model <- readRDS("inst/model_somas.rds")

  #To compute the KL divergence between the clusters 1 and 2 run
  KL <- clusters_KL_divergence(model = model, cluster_p = 1, cluster_q = 2, is_dir = dataset$is_dir)

  #Simulate 20 somas from the third cluster of the model
  simulated_somas <- sampleSomas(model = model, is_dir = dataset$is_dir, cluster_id = 3, num_simulations = 20)

  #Generate the three-dimensional representation of the first simulated soma and do not smooth. Changing number of iterations the surface of the soma smooths
  mesh <- simulation3Dmesh(simulated_somas = simulated_somas, idx = 1)

  #Generate the three-dimensional representation of the soma
  shade3d(mesh,col="red")
}
