/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Gaussian.hpp
 * Author: universidad
 *
 * Created on July 25, 2017, 12:58 PM
 */

#ifndef GAUSSIAN_HPP
#define GAUSSIAN_HPP
 
#include "Distribution.hpp"

class Gaussian : public Distribution{
    static constexpr double denom_const = sqrt(2 * M_PI);
    public:
        Gaussian (int id, int num_categories);
        Gaussian (int id, std::vector<Eigen::VectorXd> m, int num_categories);
        Gaussian (int id, std::vector<Eigen::VectorXd> m, Eigen::VectorXd sd, int num_categories);
        Gaussian (int id, std::vector<Eigen::VectorXd> m, Eigen::VectorXd sd, std::vector<int> parents, int num_categories);
    
        ~Gaussian(){}
  
        Eigen::MatrixXd density(Eigen::MatrixXd &data, Eigen::MatrixXd &weights);
        Eigen::MatrixXd log_density(Eigen::MatrixXd &data, Eigen::MatrixXd &weights);
        Eigen::MatrixXd log_likelihood(Eigen::MatrixXd &data, Eigen::MatrixXd &weights);
        double BIC(Eigen::MatrixXd &data, Eigen::MatrixXd &weights);
        void mle(Eigen::MatrixXd &data, Eigen::MatrixXd &weights);
        int get_num_params();
        std::vector<Eigen::VectorXd> get_params();
        void set_params(std::vector<Eigen::VectorXd> &params);

        //Get functions
        std::vector<Eigen::VectorXd> get_mu();
        Eigen::VectorXd get_sd();
        int get_id();
        std::vector<int> get_idx_parents();
//
//        //Set functions
//        void set_mu(Eigen::VectorXd mu);
//        void set_sd(double sd);
        void set_idx_parents(std::vector<int> &idx_parents);
    
    private:
        std::vector<Eigen::VectorXd> mu; //\Beta_0, \Beta_1, ... as parent it has
        Eigen::VectorXd sd; //\sigma
        std::vector<int> idx_parents; //An integer denoting the column in the data of the parents
        int id;
}; //Class vonMises


#endif /* GAUSSIAN_HPP */

