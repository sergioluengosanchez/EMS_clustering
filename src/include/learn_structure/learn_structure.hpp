///*
// * To change this license header, choose License Headers in Project Properties.
// * To change this template file, choose Tools | Templates
// * and open the template in the editor.
// */
//
///*
// * File:   learn_structure.hpp
// * Author: sergio
// *
// * Created on August 21, 2017, 2:05 PM
// */
//
#ifndef LEARN_STRUCTURE_HPP
#define LEARN_STRUCTURE_HPP

#include <boost/graph/adjacency_matrix.hpp>
#include <tuple>

#include "../scores.hpp"
#include "update_cache.hpp"
#include "hc_one_step.hpp"

void learn_structure(Graph& structure, Eigen::MatrixXd &data, Eigen::MatrixXd &linear_data, Eigen::MatrixXd &weights, Eigen::VectorXd &is_directional, int max_parents, bool directional_parents)
{
	Eigen::MatrixXd cache = Eigen::MatrixXd::Constant(num_vertices(structure), num_vertices(structure),-std::numeric_limits<double>::max());
	Eigen::MatrixXd node_data;
	std::vector<int> update; // vector with 100 ints.
	//Only non linear variables can be updated.

	for(int i = 0; i < num_vertices(structure); ++i)
	{
		if(is_directional(i) == 0)
		{
			update.push_back(i);
		}
	}
	//std::iota (std::begin(update), std::end(update), 0); //Fill with 0, 1, 2, ...

	std::vector<int> best_operation, parents;

	Distribution* node;
	int i, j, op_id;
	bool improve=true;

	//Get initial loglik of each node
    Eigen::VectorXd node_score = compute_node_BIC(structure, data, linear_data, weights, is_directional);
//    std::cout << "----Last node_score---: " << node_score << std::endl;
    while(improve)
    {
    	update_cache(structure, data, linear_data, weights, is_directional, cache, update, node_score, max_parents, directional_parents);
    	update.clear();
		best_operation = hc_one_step(structure, cache, node_score);
		i = best_operation.at(0);
		j = best_operation.at(1);
		op_id = best_operation.at(2);

		//std::cout<<"Operation :" << op_id << "Node j: "<< j << "Node i:" << i << std::endl;

		if(op_id == 0)
		{
			boost::add_edge(j, i, structure);
			update.push_back(i);
			node = structure[i];
			parents = node->get_idx_parents();
			parents.push_back(j);
			node->set_idx_parents(parents);
			node_data = get_node_data(node, data, linear_data, is_directional);
			node->mle(node_data, weights);

			node_score(i) = node_score(i) + cache(i,j);
		}

		if(op_id == 1)
		{
			boost::remove_edge(j, i, structure);
			update.push_back(i);

			node = structure[i];
			parents = node->get_idx_parents();
			parents.erase(std::remove(parents.begin(), parents.end(), j), parents.end());
			node->set_idx_parents(parents);
			node_data = get_node_data(node, data, linear_data, is_directional);
			node->mle(node_data, weights);

			node_score(i) = node_score(i) + cache(i,j);
		}

		if(op_id == 2)
		{
			boost::remove_edge(j, i, structure);
			boost::add_edge(i, j, structure);

			node = structure[i];
			parents = node->get_idx_parents();
			parents.erase(std::remove(parents.begin(), parents.end(), j), parents.end());
			node->set_idx_parents(parents);
			node_data = get_node_data(node, data, linear_data, is_directional);
			node->mle(node_data, weights);

			node = structure[j];
			parents = node->get_idx_parents();
			parents.push_back(i);
			node->set_idx_parents(parents);
			node_data = get_node_data(node, data, linear_data, is_directional);
			node->mle(node_data, weights);

			node_score(i) = node_score(i) + (cache(i,j) + cache(j,i));

			update.push_back(i);
			update.push_back(j);
		}

		if(op_id == -1)
		{
			improve=false;
		}

//		std::cout << "BIC is: " << compute_BN_BIC(structure, data, linear_data, weights, is_directional) << std::endl;
    }
}
#endif /* LEARN_STRUCTURE_HPP */

