/*
 * hc_one_step.hpp
 *
 *  Created on: May 24, 2018
 *      Author: sergio
 */

#ifndef HC_ONE_STEP_HPP
#define HC_ONE_STEP_HPP

#include "../distributions/Gaussian.hpp"
#include "../scores.hpp"

//This function return the node i, node j,and operation id
std::vector<int> hc_one_step(Graph &structure, Eigen::MatrixXd &cache, Eigen::VectorXd &node_score)
{
	std::vector<int> best_operation(3,-1);
	std::vector<int> parents;
	Distribution *node;
	boost::graph_traits<Graph>::vertex_iterator i, j, end_i, end_j;
	double best_score = 0;

	//For each vertex in the graph
    for(boost::tie(i,end_i) = vertices(structure); i != end_i; ++i) {
    	node = structure[*i];
    	parents = node->get_idx_parents();
    	//For each vertex in the graph
    	for(boost::tie(j,end_j) = vertices(structure); j != end_j; ++j) {
    		//If i and j are not the same vertex
    		if(*i != *j)
    		{
    			//If j is not a parent of i check if adding as parent yields a cycle in the graph.
    			//If not, check if the value in the cache is better than the best_score found until this point.
            	if(!(std::find(parents.begin(), parents.end(), *j) != parents.end()))
            	{
            		boost::add_edge(*j, *i, structure);
            		if(check_DAG(structure))//Check if the addition generates cycles in the structure
            		{
            			if(cache(*i,*j) > best_score)
            			{
            				best_operation[0] = *i;
            				best_operation[1] = *j;
            				best_operation[2] = 0;//0 is the ID of add an arc from best_operation[1] to best_operation[0]
            				best_score = cache(*i,*j);
            			}//End IF
            		}
            		boost::remove_edge(*j, *i, structure);
            	}else{//If j is parent of i then the arc can be removed or reversed
            		//Remove arc
            		if(cache(*i,*j) > best_score)
            		{
						best_operation[0] = *i;
						best_operation[1] = *j;
						best_operation[2] = 1;//0 is the ID of add an arc from best_operation[1] to best_operation[0]
						best_score = cache(*i,*j);
					}//End IF

            		//Reverse arc
            		boost::add_edge(*i, *j, structure);
            		if(check_DAG(structure))
            		{
            			if((cache(*i,*j)+cache(*j,*i)) > best_score)
						{
							best_operation[0] = *i;
							best_operation[1] = *j;
							best_operation[2] = 2;//0 is the ID of add an arc from best_operation[1] to best_operation[0]
							best_score = cache(*i,*j)+cache(*j,*i);
						}//End IF
            		}
            		boost::remove_edge(*i, *j, structure);
            	}//End IF-ELSE
    		}//End IF
    	}//End FOR
    }//End FOR
	return(best_operation);
}


#endif /* INCLUDE_LEARN_STRUCTURE_HC_ONE_STEP_HPP_ */
