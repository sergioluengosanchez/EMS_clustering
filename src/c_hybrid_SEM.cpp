/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   main.cpp
 * Author: sergio
 *
 * Created on July 22, 2017, 3:44 PM
 */


#include <Rcpp.h>
#include <boost/config.hpp>
#include <iostream>                         // for std::cout
#include <utility>                          // for std::pair
#include <algorithm>                        // for std::for_each
#include <boost/utility.hpp>                // for boost::tie
#include <boost/graph/graphviz.hpp>


#include <chrono>
#include <iostream>
#include <Eigen/Dense>

#include "include/utilities.hpp"
#include "include/validation/check_DAG.h"
#include "include/readCSV/matrix2graph.h"
#include "include/learn_structure/learn_structure.hpp"
#include "include/distributions/Gaussian.hpp"
#include "include/learn_structure/structure_initialization.hpp"
#include "include/EM.hpp"

using namespace boost;
using namespace Rcpp;

NumericMatrix MatrixXd2NumericMatrix(Eigen::MatrixXd matrix)
{
  int m = matrix.rows(), n = matrix.cols();
  NumericMatrix ans(m,n);
  for(int i = 0; i < m; ++i)
  {
    for(int j = 0; j < n; ++j)
    {
      ans(i,j) = matrix(i,j);
    }
  }
  return(ans);
}

// [[Rcpp::export]]
List c_hybrid_SEM(Rcpp::NumericMatrix data_r, Rcpp::NumericVector is_directional_r, int num_clusters, int max_parents, bool directional_parents)
{
    //Cast Rcpp object to Eigen
    Eigen::MatrixXd data_csv = Eigen::MatrixXd::Map(data_r.begin(), data_r.nrow(), data_r.ncol());
    Eigen::VectorXd is_directional = Eigen::VectorXd::Map(is_directional_r.begin(), is_directional_r.size());

    int num_data_columns = data_csv.cols();
    if(num_clusters > 0)
    {
      ++num_data_columns;
    }
    Eigen::MatrixXd data(data_csv.rows(), num_data_columns);

    //If the number of clusters is bigger than 0 the initialization is performed in the next lines. If not,
    //it is assumed that the initialization is included as part of the dataset in the last column.
    if(num_clusters > 0)
    {
      //Complete the data randomly with a number of categories
      data.topLeftCorner(data_csv.rows(),data_csv.cols()) = data_csv;
      Eigen::VectorXi category_vector = get_random_categories(num_clusters,data.rows());
      data.col(data.cols()-1) = category_vector.cast<double>();
    }else{
      data = data_csv;
    }

    //Function to transform directional variables X to cos(X), sin(X). They are used to compute
    //the Mardia-Shutton model efficiently.
    Eigen::MatrixXd linear_data = get_linear_dataset(data, is_directional);
    int num_categories = data.col(data.cols()-1).maxCoeff() + 1;

    //Function to transform the labels to a multinomial matrix.
    Eigen::MatrixXd weights = class2multinomial(data,num_categories);

    //Initialization of nodes and structure
    std::vector<Distribution*> nodes = node_initialization(data, linear_data, weights, num_categories, is_directional);
    Graph structure = structure_initialization(data, nodes);

    double current_BIC, previous_BIC;
    current_BIC = compute_BIC_observed_data(structure, data, linear_data, weights, is_directional);
    previous_BIC = current_BIC;

    auto start = std::chrono::high_resolution_clock::now();
    double observed_loglik, temp_observed_loglik;
    Eigen::MatrixXd temp_weights;

    weights = expectation_step(structure, data, linear_data, weights, is_directional, observed_loglik);
    temp_observed_loglik = observed_loglik;
    std::cout <<observed_loglik << std::endl;
    bool improved = true, improved_structure =true;
    int count = 0, count_structure = 0, count_EM = 0;

    while(improved_structure && !isnan(observed_loglik) && !isnan(current_BIC))
    {
    	learn_structure(structure, data, linear_data, weights, is_directional, max_parents, directional_parents);
    	current_BIC = compute_BIC_observed_data(structure, data, linear_data, weights, is_directional);
    	if(std::abs((current_BIC-previous_BIC)/previous_BIC)<0.001)
  		{
  			improved_structure = false;
  		}
    	count_EM = 0;
  		while(improved && count_EM < 1000 && !isnan(observed_loglik))
  		{
  			maximization_step(structure, data, linear_data, weights, is_directional);
  			temp_weights = expectation_step(structure, data, linear_data, weights, is_directional, observed_loglik);
  			if(std::abs((observed_loglik-temp_observed_loglik)/temp_observed_loglik)<0.001)
  			{
  				improved=false;
  			}else{
  				weights=temp_weights;
  			}
  			temp_observed_loglik = observed_loglik;
  			++count; ++count_EM;
  		}
      std::cout << "Current BIC is: " << current_BIC << std::endl;
  		previous_BIC = current_BIC;
  		++count_structure;
    }

    if(!isnan(observed_loglik) && !isnan(current_BIC))
    {
      auto stop = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed = stop - start;
      std::cout << "Elapsed time is: " << elapsed.count() << std::endl;
      std::cout << "Num times structure is learnt " << count_structure << " Num EM iterations:" << count<< std::endl;

      double finalBIC = compute_BN_BIC(structure, data, linear_data, weights, is_directional);
      std::cout << "Final BIC is: "<< finalBIC << std::endl;
      Eigen::MatrixXd adjacency_matrix = graph2matrix(structure);

      // Cast Eigen::MatrixXd to R type object
      NumericMatrix ans = MatrixXd2NumericMatrix(adjacency_matrix);
      NumericMatrix ans2 = MatrixXd2NumericMatrix(weights);

      //Return the structure and the responsabilities of each cluster.
      return (List::create(_["adj_m"]=ans, _["weights"] = ans2, _["BIC"] = finalBIC));
    }else{
      std::cout << "Error. Not enough instances in at least one cluster to estimate the parameters." << std::endl;
      return (List::create(_["adj_m"]=-1, _["weights"] = -1, _["BIC"] = -std::numeric_limits<double>::max()));
    }

}
