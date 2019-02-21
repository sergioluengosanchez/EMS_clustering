/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   utilities.hpp
 * Author: sergio
 *
 * Created on September 18, 2017, 2:16 PM
 */

#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "validation/check_DAG.h"

void print_edges(Graph structure)
{
    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = edges(structure); ei != ei_end; ++ei) std::cout<<"("<<source(*ei,structure)<<", "<<target(*ei,structure)<<"), ";
    std::cout << std::endl;
}

double A1inv(double x)
{
  if(0 <= x & x < 0.53) {
    return(2 * x + pow(x,3) + (5 * pow(x,5))/6);
  }else if(x < 0.85)
  {
    return(-0.4 + 1.39 * x + 0.43/(1 - x));
  } else
  {
    return(1/(pow(x,3) - 4 * pow(x,2) + 3 * x));
  }
};

double Exp(double x) // the functor we want to apply
{
    return std::exp(x);
}

Eigen::VectorXd logSumExp(Eigen::MatrixXd &x)
{
	Eigen::VectorXd max_by_row = x.rowwise().maxCoeff();
	Eigen::MatrixXd y = x.colwise()-max_by_row;
	return(max_by_row.array() + y.unaryExpr(&Exp).rowwise().sum().array().log());
}

Eigen::MatrixXd class2multinomial(Eigen::MatrixXd &data, int num_categories)
{
	Eigen::MatrixXd weights = Eigen::MatrixXd::Zero(data.rows(),num_categories);

	for(int i = 0; i < data.rows(); ++i)
	{
		weights(i,data(i,data.cols()-1))=1;
	}

	return(weights);
}

//Function to transform directional variables X to cos(X), sin(X). They are used to compute
//the Mardia-Shutton model efficiently.
Eigen::MatrixXd get_linear_dataset(Eigen::MatrixXd &data, Eigen::VectorXd &is_directional)
{
	int ncols = (is_directional.array() + 1).sum();
	Eigen::MatrixXd linear_dataset(data.rows(),ncols);

	int col_position = 0;
	for(int i = 0; i < data.cols(); ++i)
	{
		if(is_directional(i) == 1)
		{
			linear_dataset.col(col_position) = data.col(i).array().cos();
			++col_position;
	 		linear_dataset.col(col_position) = data.col(i).array().sin();
	 	}else{
	 		 linear_dataset.col(col_position) = data.col(i);
		}
		++col_position;
	}

	return(linear_dataset);
}

Eigen::MatrixXd get_node_data(Distribution* node, Eigen::MatrixXd &data, Eigen::MatrixXd &linear_data, Eigen::VectorXd &is_directional)
{
	std::vector<int> parents = node->get_idx_parents();
	int ncols = 1; //The number of columns are the parents + the node

	//Get the size of the dataset
	for(auto parent = parents.begin(); parent != parents.end(); ++parent)
	{
		if(is_directional(*parent) == 1)
		{
			ncols += 2;
		}else{
			++ncols;
		}
	}

	//Generate the dataset
	Eigen::MatrixXd node_data(data.rows(),ncols);
	int col_position = 0;
	for(auto parent = parents.begin(); parent != parents.end(); ++parent)
	{
		if(is_directional(*parent) == 1)
		{
			node_data.col(col_position) = linear_data.col((*parent) * 2);
			++col_position;
			node_data.col(col_position) = linear_data.col(((*parent) * 2) + 1);

		}else{
			node_data.col(col_position) = data.col(*parent);
		}
		++col_position;
	}

	node_data.col(col_position) = data.col(node->get_id());

	return(node_data);
}

Eigen::VectorXi get_random_categories(int num_categories, int length_vector)
{
	  Eigen::VectorXd categories = (Eigen::VectorXd::Random(length_vector)+Eigen::VectorXd::Ones(length_vector))*num_categories/2;
	  Eigen::VectorXi vi = categories.cast<int>();
	  return (vi);
}


void writeToCSVfile(std::string name, Eigen::MatrixXd matrix)
{
	Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

    std::ofstream file(name.c_str());
    file << matrix.format(CSVFormat);
 }

#endif /* UTILITIES_HPP */

