#ifndef LOCAL_SEARCH_HPP
#define LOCAL_SEARCH_HPP


#include "../distributions/Gaussian.hpp"
#include "../scores.hpp"

void update_cache(Graph &structure, Eigen::MatrixXd &data, Eigen::MatrixXd &linear_data, Eigen::MatrixXd &weights, Eigen::VectorXd &is_directional, Eigen::MatrixXd &cache, std::vector<int> update, Eigen::VectorXd &node_score, int max_parents, bool directional_parents)
{
	Distribution* node;
	std::vector<int> parents, parents_temp;
	std::vector<Eigen::VectorXd> temp_params;
	Eigen::MatrixXd node_data;

	//For each node to be updated
	for(auto i = update.begin(); i != update.end(); ++i)
	{
		node = structure[*i];
		parents = node->get_idx_parents();
		parents_temp = parents;

		//Try to set the remaining nodes as parent, or remove them if their are already a parent.
		//Compute the score of the operation
		for(int j = 0; j < cache.rows(); ++j)
		{
			if(*i != j && !(is_directional[j] == 1 && directional_parents == false))
			{
				temp_params = node->get_params();
				if((!(std::find(parents.begin(), parents.end(), j) != parents.end())) && (parents.size() < max_parents))//If j is not a parent of i. Compute the score of adding as parent
				{
          //Set j as parent of i
					parents_temp.push_back(j);
					node->set_idx_parents(parents_temp);
					node_data = get_node_data(node, data, linear_data, is_directional);
					node->mle(node_data, weights);

					//Save the BIC score of node i after adding j as parent
					cache(*i, j) = node->BIC(node_data, weights) - node_score[*i];

					//Reset the state of node i
					parents_temp = parents;
					node->set_idx_parents(parents);
				}else{
					parents_temp.erase(std::remove(parents_temp.begin(), parents_temp.end(), j), parents_temp.end());
					node->set_idx_parents(parents_temp);
					node_data = get_node_data(node, data, linear_data, is_directional);
					node->mle(node_data, weights);

					//Save the BIC score of node i after removing j as parent
					cache(*i, j) = node->BIC(node_data, weights) - node_score[*i];

					//Reset the state of node i
					parents_temp = parents;
					node->set_idx_parents(parents);
				}//End IF-ELSE
				node->set_params(temp_params);
			}//End IF
		}//End FOR
	}//End FOR
}
#endif /* LOCAL_SEARCH_HPP */
