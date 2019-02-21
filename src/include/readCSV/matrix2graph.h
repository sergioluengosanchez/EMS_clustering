/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   matrix2graph.h
 * Author: sergio
 *
 * Created on August 18, 2017, 10:48 AM
 */

#ifndef MATRIX2GRAPH_H
#define MATRIX2GRAPH_H
#include <Eigen/Dense>
#include <iostream>
#include "../validation/check_DAG.h"

//typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::directedS > Graph;
//  typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;



Graph matrix2graph(Mat adjacency_matrix)
{
    Graph g(adjacency_matrix.cols());
    
    for(int i=0;i<adjacency_matrix.cols();++i)
    {
        for(int j=0;j<adjacency_matrix.cols();++j)
        {
            if(adjacency_matrix(i,j)==1)
            {
                    add_edge((Vertex) i,(Vertex) j,g);   
            }
        }
    }
    
    if(check_DAG(g))
    {
        return(g);
    }else{
        std::cout<<"The adjacency matrix is not a DAG, a disconnected graph is returned instead"<<std::endl;
       Graph empty(adjacency_matrix.cols());
       return(empty); 
    }

};

Eigen::MatrixXd graph2matrix(Graph &structure)
{
    Eigen::MatrixXd adjacency_matrix = Eigen::MatrixXd::Zero(num_vertices(structure), num_vertices(structure));

    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = edges(structure); ei != ei_end; ++ei){
    	adjacency_matrix(source(*ei,structure),target(*ei,structure)) = 1;
    }

    return (adjacency_matrix);
};


#endif /* MATRIX2GRAPH_H */

