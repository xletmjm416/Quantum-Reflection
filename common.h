#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>
#include "Eigen/Dense"

#define PI_CONST 3.14159265358979323846
typedef std::complex<double> dcplx;
const dcplx dcplx_i(0.0,1.0);

typedef Eigen::MatrixXcd MatrixC; 	// complex matrix
typedef Eigen::VectorXcd VectorC; 	// complex vector
typedef Eigen::MatrixXd MatrixR; 	// real matrix
typedef Eigen::VectorXd VectorR; 	// real vector

#endif