#ifndef BASE
#define BASE

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

namespace utils {
	VectorR linspace(double step_x, int start = 0, int end = 1);
	
	MatrixR three_pt_diff(int size);
	
	VectorC linear_solver(MatrixC coefficients, MatrixC rhs);
}

class QMSystem {
	/* Whole QM System with the following parameters:
	space
	hamiltonian
	step x
	step t
	number of steps x
	number of steps t
	system size
	h bar
	mass
	*/
	public:
		QMSystem(VectorC, VectorR, double, double, double h_bar = 1, double mass = 1, int x_size = 1, int t_size = 1);
		MatrixC hamiltonian(VectorR);
		VectorC cranknicolson();
		VectorC get_state() { return wavefunction; }
	private:
		// physical
		VectorC wavefunction;
		VectorR potential;
		MatrixC hamilton;
		double h_bar, mass;			// physical parameters
		
		// computational
		double step_x, step_t;		// discretisation parameters
		int N_space, N_time;	// number of discrete steps
		int x_size;					// total system size: 1 is default
		int t_size;					// time unit: 1 is default
};

dcplx sin_state(dcplx x);
dcplx gauss_state(dcplx x);
#endif
