/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   cycles.h
 * Author: sergio
 *
 * Created on August 18, 2017, 12:40 PM
 */

#ifndef CHECK_DAG_H
#define CHECK_DAG_H

#include <boost/functional/hash.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/graph_utility.hpp> 
#include "../distributions/Distribution.hpp"

typedef boost::adjacency_list<boost::listS, boost::vecS, boost::bidirectionalS, Distribution*> Graph;
typedef Eigen::MatrixXd Mat;
typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;

  struct cycle_detector : public boost::dfs_visitor<>
  {
    cycle_detector(bool& has_cycle) 
      : _has_cycle(has_cycle) { }

    template <class Edge, class Graph>
    void back_edge(Edge, Graph&) {
      _has_cycle = true;
    }
    bool get_cycle(){
        return _has_cycle;
    }
  protected:
    bool& _has_cycle;
  };

  bool check_DAG(const Graph& g)
  {
    bool cycle = false;

    /*Check if there is a cycle in the graph*/
    cycle_detector vis(cycle);
    boost::depth_first_search(g, visitor(vis));
    return (!cycle); 
  }; 
  
#endif /* CHECK_DAG_H */

