/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ReadMatrix.hpp
 * Author: sergio
 *
 * Created on August 17, 2017, 1:46 PM
 */

#ifndef READMATRIX_HPP
#define READMATRIX_HPP

#include <string>
#include <Eigen/Dense>
#include "CSVIterator.hpp"

typedef Eigen::MatrixXd Mat;

Mat read_csv_matrix(std::string file_path)
{
    std::ifstream file(file_path);
    std::vector<double> param;
    
    int num_elem=0,num_row=0;
    
    std::string::size_type sz; 
    

    CSVIterator loop(file);
    num_elem=loop->size();
    loop++;
    
    for(loop; loop != CSVIterator(); ++loop)
    {
        num_elem=(*loop).size();
        for(auto i=0; i < num_elem; ++i)
        {
            auto prueba= std::stod ((*loop)[i],&sz);
             param.push_back(std::stod ((*loop)[i],&sz));
        }
        
        ++num_row;
    }
    
    Mat matrix(num_elem,num_row);

    for(auto i = 0; i < param.size(); i++) matrix(i) = param.at(i);
    matrix.transposeInPlace();
    
    return matrix;
}


#endif /* READMATRIX_HPP */

