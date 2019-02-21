## Introduction
Neural somas perform most of the metabolic activities in the neuron and support the chemical process that generates the basic elements of the synapses, and consequently, of the brain activity. The morphology of the somas is one of the fundamental features for classifying neurons and their functionality. In this study, we characterize the morphology of 39 three-dimensional reconstructed human pyramidal somas in terms of their multiresolutional Reeb graph representation, from which we extract a set of directional and linear variables to perform model-based clustering. To deal with this dataset, we introduce the Extended Mardia-Sutton mixture model whose mixture components are distributed according to a newly proposed multivariate probability density function that is able to capture the directional-linear correlations. We exploit the capabilities of Bayesian networks in combination with the Structural Expectation-Maximization algorithm to learn the finite mixture model that clusters the neural somas by their morphology and the conditional independence constraints between variables. We also derive the Kullback-Leibler divergence of the Extended Mardia-Sutton distribution to be used as a measure of similarity between soma clusters. The proposed finite mixture model discovered three subtypes of human pyramidal somas. We performed Weltch t-tests and Watson-Williams tests, as well as rule-based identification of clusters to characterize each group by its most prominent features. Furthermore, the resulting model allowed us to simulate 3D virtual representations of somas from each cluster, which can be a useful tool for neuroscientists to reason and suggest new hypotheses.


## Prerequirements and installing guide
This software has been developed as a R > 3.0.0 package. To speed-up the Structural Expectation-Maximisation algorithm, the package uses some of the C++ 14 standard features, so it needs a compiler with C++ 14 support. Consequently, it is needed an R environment and internet connectivity to download additional package dependencies. R software can be downloaded from <https://cran.r-project.org/>. 

The C++ compiler we have used is
|Compiler|Version|Tested|
|--------|-------|------|
|     g++|     >5| 5.4.1|

Also Rstudio IDE is needed to compile the source files.
|IDE     |Version|Tested |
|--------|-------|-------|
| Rstudio|   >0.9| 1.0.136|

We provide the steps for a clean installation in Ubuntu 16.04. This software has not been tried under Windows.

The package also uses the following dependencies: 
 
To install **CGAL** follow the next steps
  1. Run sudo apt-get install libcgal-dev
Also can be downloaded from the official [website](https://www.cgal.org/index.html)

To install **Eigen3** follow the next steps
  
  1. Run sudo apt install libeigen3-dev

  OR
  
  1. Download Eigen3 from its [website](http://eigen.tuxfamily.org)
  2. In the directory where you want put the Eigen3 installation uncompressing the downloaded file
  3. Make a directory named *build* inside the installation folder
  4. Run `cmake ../` from the shell in the recently created folder *build*
  5. Run `sudo make install` from the shell
  6. Edit the Makevars of the R project to include the path of the headers files

The compiler search the library in the path /usr/include/eigen3. If it is installed in other path update the file src/Makevars with the new path.

Some R packages are needed to perform some specific tasks releated with 3D processing, data management, or modeling. They must be installed thorugh the command `install.packages("name_of_the_package")` to use this package. The R dependencies of the package are:
|Package    |Version|
|-----------|-------|
|foreach    |  1.4.4|
|gplots     |3.0.1.1|
|mclust     |  5.4.2|
|Morpho     |    2.6|
|Rcpp       |  1.0.0| 
|rgl        |0.99.16|
|Rvcg       |   0.18|
|parallel   |  3.5.2|
|tools      |  3.5.2|
|matrixStats| 0.54.0|
|pracma     |  2.2.2|
|geometry   |  0.4.0|
|misc3d     |  0.8.4|
|bnlearn    |    4.4|
|caret      | 6.0-81|
|truncnorm  |  1.0-8|
|doParallel | 1.0.14|
|data.table | 1.12.0|
|R6         |  2.4.0|
|circular   | 0.4-93|
|msm        |  1.6.6|
|clue       | 0.3-56|

Updated versions of the R dependencies packages should be supported.

To compile the project, open **SomaSimulation.Rproj**. When Rstudio is open press Cntrl+Shift+B to compile the project.

## example.R
File **R/example.R** provides a demo that shows how to use the code to: compute the insertion points and the morphometric features of the somas, cluster the data according to the Extended Mardia-Sutton finite mixture model, compute the KL divergence between clusters, and simulate three-dimensional somas. 

## Three-dimensional examples of the cluster
The three-dimensional representations of the clusters assigned to each cluster can be found in the file **PDF_3DClusters.zip**
