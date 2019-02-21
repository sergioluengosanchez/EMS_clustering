#include <iostream>
#include <limits>
#include <boost/math/special_functions/bessel.hpp>
#include "include/distributions/VonMises.hpp"

using namespace std;

//Constructor
VonMises::VonMises(int id, int num_categories) : id(id) {
    Eigen::VectorXd temp(1);
    temp<<1;
    for(auto i=0; i<num_categories; ++i)
    {
    	this->mu.push_back(temp);
    }
    this->kappa=Eigen::VectorXd::Ones(num_categories);
}


VonMises::VonMises(int id, std::vector<Eigen::VectorXd> m, int num_categories) : id(id), mu(m) {
	this->kappa=Eigen::VectorXd::Constant(num_categories,1);
}

VonMises::VonMises(int id, std::vector<Eigen::VectorXd> m, Eigen::VectorXd kappa, int num_categories) : id(id), mu(m), kappa(kappa) {}

VonMises::VonMises(int id, std::vector<Eigen::VectorXd> m, Eigen::VectorXd kappa, std::vector<int> idx_parents, int num_categories) : id(id), mu(m), kappa(kappa), idx_parents(idx_parents) {}

//Compute density of the von Mises Distribution
Eigen::MatrixXd VonMises::density(Eigen::MatrixXd &data, Eigen::MatrixXd &weights){
	double denom_const = 2*M_PI;
	int num_categories = weights.cols();
	Eigen::MatrixXd density(data.rows(),num_categories);

	//For each category compute the density
	for(int categories = 0; categories < num_categories; ++categories)
	{
		density.col(categories) = (data.col(data.cols()-1).array() - this->mu.at(categories)(0)).array().cos() * this->kappa(categories);
		density.col(categories) = density.col(categories).array().exp() / (denom_const * boost::math::cyl_bessel_i(0,this->kappa(categories)));
	}
	return(density);
};

Eigen::MatrixXd VonMises::log_density(Eigen::MatrixXd &data, Eigen::MatrixXd &weights){
	double denom_const = log(2*M_PI);
	int num_categories = weights.cols();
	Eigen::MatrixXd loglik(data.rows(),num_categories);

	//For each category compute the density
	for(int categories = 0; categories < num_categories; ++categories)
	{
		loglik.col(categories) = (data.col(data.cols()-1).array() - this->mu.at(categories)(0)).array().cos() * this->kappa(categories);
	  if(this->kappa(categories) < 700)
	  {
	    loglik.col(categories) = loglik.col(categories).array() - denom_const - log(boost::math::cyl_bessel_i(0,this->kappa(categories)));
	  }else{
	    loglik.col(categories) = loglik.col(categories).array() - std::numeric_limits<double>::quiet_NaN();
	  }
	}
	return(loglik);
};

Eigen::MatrixXd VonMises::log_likelihood(Eigen::MatrixXd &data, Eigen::MatrixXd &weights){
	double denom_const = log(2*M_PI);
	int num_categories = weights.cols();
	Eigen::MatrixXd loglik(data.rows(),num_categories);

	//For each category compute the density
	for(int categories = 0; categories < num_categories; ++categories)
	{
		loglik.col(categories) = (data.col(data.cols()-1).array() - this->mu.at(categories)(0)).array().cos() * this->kappa(categories);
		loglik.col(categories) = loglik.col(categories).array() - denom_const - log(boost::math::cyl_bessel_i(0,this->kappa(categories)));
		loglik.col(categories) = weights.col(categories).cwiseProduct(loglik.col(categories));
	}
	return(loglik);
};

double VonMises::BIC(Eigen::MatrixXd &data, Eigen::MatrixXd &weights){
    double log_lik = VonMises::log_likelihood(data, weights).sum();
    double penalization = this->get_num_params() * log(data.rows()) * 0.5;
    return (log_lik - penalization);
};

//Aux function to compute mle with Eigen. It is an approximation to compute kappa
double VonMises::A1inv(double x)
{
  if(0 <= x && x < 0.53) {
    return(2 * x + pow(x,3) + (5 * pow(x,5))/6);
  }else if(x < 0.85){
    return(-0.4 + 1.39 * x + 0.43/(1 - x));
  } else{
    return(1/(pow(x,3) - 4 * pow(x,2) + 3 * x));
  }
};

void VonMises::mle(Eigen::MatrixXd &data, Eigen::MatrixXd &weights) {
	int num_categories = weights.cols();
	Eigen::VectorXd priori = weights.colwise().sum();
	Eigen::VectorXd residuals;
	for(int category = 0; category < num_categories; ++category)
	{
		this->mu[category](0) = std::atan2(weights.col(category).dot(data.col(data.cols()-1).array().sin().matrix()),weights.col(category).dot(data.col(data.cols()-1).array().cos().matrix()));
		residuals = data.col(data.cols()-1).array() - this->mu[category](0);
		this->kappa[category] = this->A1inv(weights.col(category).dot(residuals.array().cos().matrix()) / priori(category));
	}
};

int VonMises::get_num_params(){
	int num_params = this->mu.size();
	return(num_params);
};


std::vector<Eigen::VectorXd> VonMises::get_params(){
	std::vector<Eigen::VectorXd> params;
	for (int i = 0; i < this->mu.size(); ++i)
	{
		params.push_back(this->mu.at(i));
	}
	params.push_back(this->kappa);
	return(params);
};


void VonMises::set_params(std::vector<Eigen::VectorXd> &params){
	Eigen::VectorXd temp_mu (params.size()-1);
	for (int i = 0; i < (params.size()-1); ++i)
	{
		this->mu[i] = params[i];
	}
	this->kappa = params[params.size()-1];
};

std::vector<int> VonMises::get_idx_parents(){
    return (this->idx_parents);
};

void VonMises::set_idx_parents(std::vector<int> &idx_parents){
    this->idx_parents=idx_parents;
};

int VonMises::get_id(){
    return (this->id);
};

//Adapted from CircStats
//float VonMises::pvm_mu0 (float theta, float acc) {
//  bool flag = true;
//  int p = 1;
//  float sum = 0;
//  float term=0;
//  while (flag) {
//    term = (boost::math::cyl_bessel_i(p, this->sd) * sin(p * theta))/p;
//    sum = sum + term;
//    p = p + 1;
//    if (std::abs(term) < acc)
//      flag = false;
//  }
//  return(theta/(2 * M_PI) + sum/(M_PI * boost::math::cyl_bessel_i(0,this->sd)));
//};



