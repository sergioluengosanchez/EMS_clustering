///*
// * To change this license header, choose License Headers in Project Properties.
// * To change this template file, choose Tools | Templates
// * and open the template in the editor.
// */
//
///* 
// * File:   GaussianMardia.hpp
// * Author: universidad
// *
// * Created on August 14, 2017, 1:30 PM
// */
//
//#ifndef GAUSSIANMARDIA_HPP
//#define GAUSSIANMARDIA_HPP
//  #include "Distribution.hpp"
//
//
//class GaussianMardia : public Distribution{
//  public:
//    GaussianMardia ();
//    GaussianMardia (float m);
//    GaussianMardia (float m, float c);
//    
//    ~GaussianMardia(){}
//    
//    float density (std::vector<float> x);
//    float log_likelihood (std::vector<std::vector<float>> data);  
//    void mle (std::vector<std::vector<float>> data);
//    float get_mu ();
//    float get_sd ();
//    
//  private:
//    float mu, mu0, sd, ro1, ro2, sigma, kappa;
//}; //Class GaussianMardia
//
//#endif /* GAUSSIANMARDIA_HPP */
//
