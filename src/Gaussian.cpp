#include <iostream>
#include <chrono>
#include "include/distributions/Gaussian.hpp"



Gaussian::Gaussian(int id, int num_categories) : id(id) {
    Eigen::VectorXd temp(1);
    temp<<1;
    for(auto i=0; i<num_categories; ++i)
    {
    	this->mu.push_back(temp);
    }
    this->sd=Eigen::VectorXd::Ones(num_categories);
}

Gaussian::Gaussian(int id, std::vector<Eigen::VectorXd> m, int num_categories) : id(id), mu(m) {
	this->sd=Eigen::VectorXd::Constant(num_categories,1);
}

Gaussian::Gaussian(int id, std::vector<Eigen::VectorXd> m, Eigen::VectorXd sd, int num_categories) : id(id), mu(m), sd(sd) {}

Gaussian::Gaussian(int id, std::vector<Eigen::VectorXd> m, Eigen::VectorXd sd, std::vector<int> idx_parents, int num_categories) : id(id), mu(m), sd(sd), idx_parents(idx_parents) {}

 /*Compute density an instance according to linear Gaussian. */
Eigen::MatrixXd Gaussian::density(Eigen::MatrixXd &data, Eigen::MatrixXd &weights){
    double denom_const = std::sqrt(2*M_PI);
    int num_coef = this->idx_parents.size() + 1;
    int num_categories = weights.cols();
	Eigen::MatrixXd parent_data(data.rows(), num_coef), density(data.rows(),num_categories);
	Eigen::VectorXd v_mean(data.rows());

	//Generate the dataset of the parents
	parent_data.col(0) = Eigen::VectorXd::Ones(data.rows());
	for(int i = 0; i < this->idx_parents.size(); ++i)
	{
		parent_data.col(i+1) = data.col(this->idx_parents.at(i));
	}

	//Compute the density of the variable
	for(int categories = 0; categories < num_categories; ++categories)
	{
		v_mean = parent_data * this->mu.at(categories);
		density.col(categories) = -(data.col(this->id) - v_mean).array().pow(2) / (2*std::pow(this->sd(categories),2));
		density.col(categories) = density.col(categories).array().exp() / (denom_const * this->sd(categories));
	}

	return(density);
};

/*Compute density an instance according to linear Gaussian. */
Eigen::MatrixXd Gaussian::log_density(Eigen::MatrixXd &data, Eigen::MatrixXd &weights){
	double denom_const = log(2*M_PI)/2;
	int num_coef = data.cols();
	int num_categories = weights.cols();
	Eigen::MatrixXd parent_data(data.rows(), num_coef), loglik(data.rows(), num_categories);
	Eigen::VectorXd v_mean(data.rows());

	//Generate the dataset of the parents
	parent_data.col(0) = Eigen::VectorXd::Ones(data.rows());
	for(int i = 0; i < (data.cols()-1); ++i)
	{
		parent_data.col(i+1) = data.col(i);
	}

	//Compute the density of the variable
	for(int categories = 0; categories < num_categories; ++categories)
	{
		v_mean = parent_data * this->mu.at(categories);
		loglik.col(categories) = -(data.col(data.cols()-1) - v_mean).array().pow(2) / (2*std::pow(this->sd(categories),2));
		loglik.col(categories) = loglik.col(categories).array() - denom_const - log(this->sd(categories));
	}

	return(loglik);
};

/*Compute loglik of an instance according to linear Gaussian. */
Eigen::MatrixXd Gaussian::log_likelihood(Eigen::MatrixXd &data, Eigen::MatrixXd &weights){
	double denom_const = log(2*M_PI)/2;
	int num_coef = data.cols();
	int num_categories = weights.cols();
	Eigen::MatrixXd parent_data(data.rows(), num_coef), loglik(data.rows(), num_categories);
	Eigen::VectorXd v_mean(data.rows());

	//Generate the dataset of the parents
	parent_data.col(0) = Eigen::VectorXd::Ones(data.rows());
	for(int i = 0; i < (data.cols()-1); ++i)
	{
		parent_data.col(i+1) = data.col(i);
	}

	//Compute the density of the variable
	for(int categories = 0; categories < num_categories; ++categories)
	{
		v_mean = parent_data * this->mu.at(categories);
		loglik.col(categories) = -(data.col(data.cols()-1) - v_mean).array().pow(2) / (2*std::pow(this->sd(categories),2));
		loglik.col(categories) = loglik.col(categories).array() - denom_const - log(this->sd(categories));
		loglik.col(categories) = weights.col(categories).cwiseProduct(loglik.col(categories));
	}

	return(loglik);
};

double Gaussian::BIC(Eigen::MatrixXd &data, Eigen::MatrixXd &weights){
    double log_lik = Gaussian::log_likelihood(data, weights).sum();
    double penalization = this->get_num_params() * log(data.rows()) * 0.5;
    return (log_lik - penalization);
};

//Compute mle assuming that the last column corresponds to the class variable
void Gaussian::mle(Eigen::MatrixXd &data, Eigen::MatrixXd &weights){
  int num_coef = data.cols(); //The number of coefficient is equal to the number of variables in data, which is the data of the parents and the variable itself
	int num_categories = weights.cols();
  int category = -1;
	Eigen::VectorXd num_instances = weights.colwise().sum();

	//Definition of the variables
	std::vector<Eigen::MatrixXd> v_coef_matrix;
	std::vector<Eigen::VectorXd> v_independent_terms;
	Eigen::MatrixXd parent_data(data.rows(), num_coef);

	//Initialization of the matrices
	for(int i = 0; i < num_categories; ++i)
	{
		v_coef_matrix.push_back(Eigen::MatrixXd::Zero(num_coef,num_coef));
		v_independent_terms.push_back(Eigen::VectorXd::Zero(num_coef));
	}

	//Generate the dataset of the parents. First column is the intercept
	parent_data.col(0) = Eigen::VectorXd::Ones(data.rows());
	for(int i = 0; i < (data.cols()-1); ++i)
	{
		parent_data.col(i+1) = data.col(i);
	}

	//For each instance compute the value of the coefficients and aggregate.
	//Carefull!! The values can be extremely big.
	for(int i = 0; i < data.rows(); ++i)
	{
		for(int category = 0; category < num_categories; ++category)
		{
			v_coef_matrix[category] += (weights(i,category) * parent_data.row(i).transpose() * parent_data.row(i));
			v_independent_terms[category] += (weights(i,category) * data(i,data.cols()-1) * parent_data.row(i));
		}

	}

	//Compute the parameters for each category
	for(int category = 0; category < num_categories; ++category)
	{
		v_coef_matrix[category] = v_coef_matrix[category] / num_instances(category);
		v_independent_terms[category] = v_independent_terms[category] / num_instances(category);

		//Solve the system of equations
		this->mu[category] = v_coef_matrix[category].colPivHouseholderQr().solve(v_independent_terms[category]);

		//Compute the residuals
		this->sd(category) = 0;
		for(int i = 0; i < data.rows(); ++i)
		{
			this->sd(category) += (weights(i,category) * std::pow(data(i,data.cols()-1) - this->mu[category].dot(parent_data.row(i)),2));
		}

		this->sd(category) = sqrt(this->sd(category)/(num_instances(category) - this->idx_parents.size() - 1));
	}
};

//Get the number of parameters in the node. It is computed as the number of (categories * num_coeffs) + num_coeffs.
int Gaussian::get_num_params(){
	int num_params = this->mu.size() * this->mu.at(0).size();
	return(num_params);
};

std::vector<Eigen::VectorXd> Gaussian::get_params(){
	std::vector<Eigen::VectorXd> params;
	for (int i = 0; i < this->mu.size(); ++i)
	{
		params.push_back(this->mu.at(i));
	}
	params.push_back(this->sd);
	return(params);
};

void Gaussian::set_params(std::vector<Eigen::VectorXd> &params){
	Eigen::VectorXd temp_mu (params.size()-1);
	for (int i = 0; i < (params.size()-1); ++i)
	{
		this->mu[i] = params[i];
	}
	this->sd = params[params.size()-1];
};

std::vector<Eigen::VectorXd> Gaussian::get_mu(){
    return (this->mu);
};

Eigen::VectorXd Gaussian::get_sd(){
    return (this->sd);
};

int Gaussian::get_id(){
    return (this->id);
};

std::vector<int> Gaussian::get_idx_parents(){
    return (this->idx_parents);
};

void Gaussian::set_idx_parents(std::vector<int> &idx_parents){
    this->idx_parents=idx_parents;
};
