/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   structure_initialization.hpp
 * Author: sergio
 *
 * Created on September 18, 2017, 1:07 PM
 */

#ifndef STRUCTURE_INITIALIZATION_HPP
#define STRUCTURE_INITIALIZATION_HPP

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include "../distributions/Distribution.hpp"
#include "../distributions/Gaussian.hpp"
#include "../distributions/VonMises.hpp"

using namespace boost;

std::vector<Distribution*> node_initialization(Eigen::MatrixXd &data, Eigen::MatrixXd &linear_data, Eigen::MatrixXd &weights, int num_categories, Eigen::VectorXd is_directional)
{
    std::vector<Distribution*> nodes;
    Eigen::MatrixXd node_data;
    for(int id = 0; id < (data.cols()-1); ++id)
    {
    	if(is_directional(id) == 1)
    	{
    		nodes.push_back(new VonMises{id, num_categories});
    	}else{
    		nodes.push_back(new Gaussian{id, num_categories});
    	}

    	node_data = get_node_data(nodes.at(id), data, linear_data, is_directional);
        nodes.at(id)->mle(node_data, weights);
    }

    return(nodes);
}

Graph structure_initialization(Eigen::MatrixXd &data, std::vector<Distribution*> &nodes)
{
    Graph initial_structure;
    Vertex v;
    for(int id = 0; id < (data.cols() - 1); ++id){
        v = boost::add_vertex(nodes.at(id), initial_structure);
    }

    return (initial_structure);
}




#endif /* STRUCTURE_INITIALIZATION_HPP */

