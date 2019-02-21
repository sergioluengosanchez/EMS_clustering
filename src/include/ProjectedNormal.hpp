///*
// * To change this license header, choose License Headers in Project Properties.
// * To change this template file, choose Tools | Templates
// * and open the template in the editor.
// */
//
///* 
// * File:   ProjectedNormal.hpp
// * Author: sergio
// *
// * Created on September 4, 2017, 11:54 AM
// */
//
//#ifndef PROJECTEDNORMAL_HPP
//#define PROJECTEDNORMAL_HPP
//
//#include <boost/math/distributions/normal.hpp>
//#include "Distribution.hpp"
//
////It is a bivariate normal but the covariance matrix is constraint to be the identity matrix
//class ProjectedNormal : public Distribution{
//    public:
//        ProjectedNormal ();
//        ProjectedNormal (float m);
//        ProjectedNormal (float m, float c);
//    
//        ~ProjectedNormal(){}
//  
//        float density(std::vector<float> x);
//        float cumulative(float x);
//        float log_likelihood(std::vector<std::vector<float>> data);  
//        void mle(std::vector<std::vector<float>> data);
//        float get_mu();
//        float get_sd();
//  
//    private:
//    boost::math::normal_distribution<float> normal;
//        
//}; //Class vonMises
//
//
//#endif /* PROJECTEDNORMAL_HPP */
//
