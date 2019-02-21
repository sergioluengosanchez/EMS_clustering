/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   VonMises.hpp
 * Author: universidad
 *
 * Created on July 25, 2017, 12:59 PM
 */

#ifndef VONMISES_HPP
#define VONMISES_HPP

#include "Distribution.hpp"


class VonMises : public Distribution{
	static constexpr double denom_const = sqrt(2 * M_PI);
	public:
		VonMises (int id, int num_categories);
		VonMises (int id, std::vector<Eigen::VectorXd> m, int num_categories);
		VonMises (int id, std::vector<Eigen::VectorXd> m, Eigen::VectorXd sd, int num_categories);
		VonMises (int id, std::vector<Eigen::VectorXd> m, Eigen::VectorXd sd, std::vector<int> parents, int num_categories);

		~VonMises(){}

		//Distribution function
		Eigen::MatrixXd density(Eigen::MatrixXd &data, Eigen::MatrixXd &weights);
		Eigen::MatrixXd log_density(Eigen::MatrixXd &data, Eigen::MatrixXd &weights);
		Eigen::MatrixXd log_likelihood(Eigen::MatrixXd &data, Eigen::MatrixXd &weights);
		double BIC(Eigen::MatrixXd &data, Eigen::MatrixXd &weights);
		void mle(Eigen::MatrixXd &data, Eigen::MatrixXd &weights);

		//Get functions
		int get_num_params();
		std::vector<Eigen::VectorXd> get_params();
		int get_id();
		std::vector<int> get_idx_parents();

		//Set functions
		void set_params(std::vector<Eigen::VectorXd> &params);
		void set_idx_parents(std::vector<int> &idx_parents);

	private:
		double pvm_mu0(double theta, double acc);
		double A1inv(double x);
		std::vector<Eigen::VectorXd> mu;
		Eigen::VectorXd kappa;
		std::vector<int> idx_parents;
		int id;
}; //Class VonMises


#endif /* VONMISES_HPP */

