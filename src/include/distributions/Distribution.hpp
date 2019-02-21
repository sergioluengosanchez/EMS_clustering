/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Distribution.hpp
 * Author: universidad
 *
 * Created on July 25, 2017, 12:58 PM
 */

#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP


#include <vector>
#include <cmath>
#include <numeric>
#include <random>
#undef Success
#include <Eigen/Dense>

class Distribution{
public:
  virtual Eigen::MatrixXd density(Eigen::MatrixXd &data, Eigen::MatrixXd &weights)=0;
  virtual Eigen::MatrixXd log_density(Eigen::MatrixXd &data, Eigen::MatrixXd &weights)=0;
  virtual Eigen::MatrixXd log_likelihood(Eigen::MatrixXd &data, Eigen::MatrixXd &weights)=0;
  virtual double BIC(Eigen::MatrixXd &data, Eigen::MatrixXd &weights)=0;
  virtual void mle(Eigen::MatrixXd &data, Eigen::MatrixXd &weights)=0;
  virtual int get_num_params()=0;
  virtual std::vector<Eigen::VectorXd> get_params()=0;
  virtual void set_params(std::vector<Eigen::VectorXd> &params)=0;
  virtual std::vector<int> get_idx_parents()=0;
  virtual void set_idx_parents(std::vector<int> &idx_parents)=0;
  virtual int get_id()=0;
  virtual ~Distribution(){}
};

#endif /* DISTRIBUTION_HPP */

