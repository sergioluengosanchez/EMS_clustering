///*
// * To change this license header, choose License Headers in Project Properties.
// * To change this template file, choose Tools | Templates
// * and open the template in the editor.
// */
//
///*
// * File:   scores.hpp
// * Author: sergio
// *
// * Created on September 18, 2017, 2:34 PM
// */
//
#ifndef SCORES_HPP
#define SCORES_HPP

#include "distributions/Gaussian.hpp"
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/config.hpp>
#include "utilities.hpp"

using namespace boost;

Eigen::VectorXd compute_node_loglik(Graph &structure, Eigen::MatrixXd &data, Eigen::MatrixXd &weights){
    Eigen::VectorXd loglik (data.cols()-1);

	int id = 0;
	for (auto node = structure.m_vertices.begin(); node != structure.m_vertices.end(); ++node) {
		id = node->m_property->get_id();
		loglik[id] = node->m_property->log_likelihood(data, weights).sum();
	}
    return(loglik);
}

double compute_loglik_observed_data(Graph &structure, Eigen::MatrixXd &data, Eigen::MatrixXd &weights){
	Eigen::MatrixXd log_lik = Eigen::MatrixXd::Zero(weights.rows(),weights.cols());

	int id = 0;
	for (auto node = structure.m_vertices.begin(); node != structure.m_vertices.end(); ++node) {
		id = node->m_property->get_id();
		log_lik += node->m_property->density(data, weights).array().log().matrix();
	}

	Eigen::VectorXd priori = (weights.colwise().sum() / data.rows()).array().log();
	log_lik.rowwise() += priori.transpose();

	Eigen::VectorXd denominator = logSumExp(log_lik);
	return(denominator.sum());
}

//double compute_BIC_observed_data(Graph &structure, Eigen::MatrixXd &data, Eigen::MatrixXd &weights){
//	Eigen::MatrixXd log_lik = Eigen::MatrixXd::Zero(weights.rows(),weights.cols());
//	int num_params_BN = 0;
//	double BIC = 0, penalization = 0;
//
//	int id = 0;
//	for (auto node = structure.m_vertices.begin(); node != structure.m_vertices.end(); ++node) {
//		id = node->m_property->get_id();
//		log_lik += node->m_property->density(data, weights).array().log().matrix();
//		num_params_BN += node->m_property->get_num_params();
//	}
//
//	Eigen::VectorXd priori = (weights.colwise().sum() / data.rows()).array().log();
//	log_lik.rowwise() += priori.transpose();
//
//	Eigen::VectorXd denominator = logSumExp(log_lik);
//	penalization = num_params_BN * log(data.rows()) * 0.5;
//	BIC = denominator.sum() - penalization;
//	return(BIC);
//}

//double compute_BIC_observed_data(Graph &structure, Eigen::MatrixXd &data, Eigen::MatrixXd &linear_data, Eigen::MatrixXd &weights, Eigen::VectorXd &is_directional){
//	Eigen::MatrixXd log_lik = Eigen::MatrixXd::Zero(weights.rows(),weights.cols());
//	int num_params_BN = 0;
//	double BIC = 0, penalization = 0;
//
//	Eigen::MatrixXd node_data;
//	Distribution* node;
//	for (auto id = 0; id < num_vertices(structure); ++id)
//	{
//		node = structure[id];
//		node_data = get_node_data(node, data, linear_data, is_directional);
//		log_lik += node->density(node_data, weights).array().log().matrix();
//		num_params_BN += node->get_num_params();
//	}
//
//	Eigen::VectorXd priori = (weights.colwise().sum() / data.rows()).array().log();
//	log_lik.rowwise() += priori.transpose();
//
//	Eigen::VectorXd denominator = logSumExp(log_lik);
//	penalization = num_params_BN * log(data.rows()) * 0.5;
//	BIC = denominator.sum() - penalization;
//	return(BIC);
//}

double compute_BIC_observed_data(Graph &structure, Eigen::MatrixXd &data, Eigen::MatrixXd &linear_data, Eigen::MatrixXd &weights, Eigen::VectorXd &is_directional){
	Eigen::MatrixXd log_lik = Eigen::MatrixXd::Zero(weights.rows(),weights.cols());
	int num_params_BN = 0;
	double BIC = 0, penalization = 0;

	Eigen::MatrixXd node_data;
	Distribution* node;
	for (auto id = 0; id < num_vertices(structure); ++id)
	{
		node = structure[id];
		node_data = get_node_data(node, data, linear_data, is_directional);
		log_lik += node->log_density(node_data, weights);
		num_params_BN += node->get_num_params();
	}

	Eigen::VectorXd priori = (weights.colwise().sum() / data.rows()).array().log();
	log_lik.rowwise() += priori.transpose();

	Eigen::VectorXd denominator = logSumExp(log_lik);
	penalization = num_params_BN * log(data.rows()) * 0.5;
	BIC = denominator.sum() - penalization;
	return(BIC);
}


//Eigen::VectorXd compute_node_BIC(Graph &structure, Eigen::MatrixXd &data, Eigen::MatrixXd &weights){
//    Eigen::VectorXd BIC (data.cols()-1);
//	int id = 0;
//	for (auto node = structure.m_vertices.begin(); node != structure.m_vertices.end(); ++node) {
//		id = node->m_property->get_id();
//		BIC[id] = node->m_property->BIC(data, weights);
//	}
//    return(BIC);
//}

Eigen::VectorXd compute_node_BIC(Graph &structure, Eigen::MatrixXd &data, Eigen::MatrixXd &linear_data, Eigen::MatrixXd &weights, Eigen::VectorXd &is_directional){
    Eigen::VectorXd BIC (data.cols()-1);
    Eigen::MatrixXd node_data;
	Distribution* node;
	for (auto id = 0; id < num_vertices(structure); ++id)
	{
		node = structure[id];
		node_data = get_node_data(node, data, linear_data, is_directional);
		BIC[id] = node->BIC(node_data, weights);
	}

    return(BIC);
}

double compute_priori_loglik(Eigen::MatrixXd &data, Eigen::MatrixXd &weights)
{
	Eigen::VectorXd priori = weights.colwise().sum() / data.rows();
	Eigen::VectorXd weighted_priori = weights * priori;
	return(weighted_priori.array().log().sum());
}

double compute_priori_BIC(Eigen::MatrixXd &data, Eigen::MatrixXd &weights)
{
	return(compute_priori_loglik(data, weights) - (log(data.rows()) / 2 * (weights.cols()-1)));
}

double compute_BN_BIC(Graph &structure, Eigen::MatrixXd &data, Eigen::MatrixXd &linear_data, Eigen::MatrixXd &weights, Eigen::VectorXd &is_directional)
{
	return (compute_priori_BIC(data, weights) + compute_node_BIC(structure, data, linear_data, weights, is_directional).sum());
}


#endif /* SCORES_HPP */

