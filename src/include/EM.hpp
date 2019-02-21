/*
 * EM.h
 *
 *  Created on: May 30, 2018
 *      Author: universidad
 */

#ifndef INCLUDE_EM_HPP_
#define INCLUDE_EM_HPP_


#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/config.hpp>
#include "distributions/Gaussian.hpp"
#include "scores.hpp"
#include "utilities.hpp"
#include <iomanip>
#include <math.h>
#include <unistd.h>

using namespace boost;


//Eigen::MatrixXd expectation_step(Graph &structure, Eigen::MatrixXd &data, Eigen::MatrixXd &weights, double &observed_loglik)
//{
//	Eigen::MatrixXd log_lik = Eigen::MatrixXd::Zero(weights.rows(),weights.cols());
//
//	int id = 0;
//	for (auto node = structure.m_vertices.begin(); node != structure.m_vertices.end(); ++node) {
//		id = node->m_property->get_id();
//		log_lik += node->m_property->density(data, weights).array().log().matrix();
//	}
//
//	Eigen::VectorXd priori = (weights.colwise().sum() / data.rows()).array().log();
//	log_lik.rowwise() += priori.transpose();
//
//	Eigen::VectorXd denominator = logSumExp(log_lik);
//	Eigen::MatrixXd new_weights = (log_lik.colwise() - denominator).array().exp();
//	observed_loglik = denominator.sum();
//	return(new_weights);
//}
//
////Compute the MLE for each variable
//void maximization_step(Graph &structure, Eigen::MatrixXd &data, Eigen::MatrixXd &weights)
//{
//	for (auto node = structure.m_vertices.begin(); node != structure.m_vertices.end(); ++node) {
//		node->m_property->mle(data, weights);
//	}
//}

Eigen::MatrixXd expectation_step(Graph &structure, Eigen::MatrixXd &data, Eigen::MatrixXd &linear_data, Eigen::MatrixXd &weights, Eigen::VectorXd &is_directional, double &observed_loglik)
{
	Eigen::MatrixXd log_lik = Eigen::MatrixXd::Zero(weights.rows(),weights.cols());
	Eigen::MatrixXd temp_loglik;
	Eigen::MatrixXd node_data;
	Distribution* node;
	for (auto id = 0; id < num_vertices(structure); ++id)
	{
		node = structure[id];
		node_data = get_node_data(node, data, linear_data, is_directional);
		temp_loglik = node->log_density(node_data, weights);
		log_lik += temp_loglik;
	}
	Eigen::VectorXd priori = (weights.colwise().sum() / data.rows()).array().log();
	log_lik.rowwise() += priori.transpose();

	Eigen::VectorXd denominator = logSumExp(log_lik);
	Eigen::MatrixXd new_weights = (log_lik.colwise() - denominator).array().exp();
	observed_loglik = denominator.sum();
	return(new_weights);
}

//Compute the MLE for each variable
void maximization_step(Graph &structure, Eigen::MatrixXd &data, Eigen::MatrixXd &linear_data, Eigen::MatrixXd &weights, Eigen::VectorXd &is_directional)
{
	Eigen::MatrixXd node_data;
	Distribution* node;
	for (auto id = 0; id < num_vertices(structure); ++id)
	{
		node = structure[id];
		node_data = get_node_data(node, data, linear_data, is_directional);
		node->mle(node_data, weights);
	}
}


#endif /* INCLUDE_EM_HPP_ */
