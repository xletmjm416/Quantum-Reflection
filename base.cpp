#include <iostream>
#include <cmath>
#include <complex>
#include "Eigen/Dense"

typedef std::complex<double> dcplx;
const dcplx dcplx_i(0.0,1.0);

dcplx sin(dcplx x) {
	return std::sin(x);
}

Eigen::VectorXd linspace(double step_x) {
	/*	Initialises the linearly divided space.
	*	step_x		- length of a step
	*/
	Eigen::VectorXd ans;
	int N_steps = (int)(1/step_x);
	ans = Eigen::VectorXd::LinSpaced(N_steps, 0, 1);
	return ans;
}

Eigen::VectorXcd map(dcplx (*func)(dcplx x), Eigen::VectorXcd vec) {
	/*	Maps a vector to a vector to which coefficient-wise function "func" was applied.
	*	func		- function f(x) which defines the wavefunction at t=0;
	*	Eigen::VectorXcd vec		- the vector to be mapped
	*/
	Eigen::VectorXcd ans;
	ans = vec.unaryExpr(std::ptr_fun(func));
	return ans;
}

Eigen::MatrixXd second_derv(int size) {
	/*	Produces a three-point difference algorithm matrix (tridiagonal with 1, -2, 1 values.
	*	size		- size of the matrix (size x size)
	*/
	Eigen::MatrixXd ans(size, size);
	ans = Eigen::MatrixXd::Zero(size, size);
	for(int row = 0; row<size; row++) {
		if(row != 0) ans(row,row-1) = 1;
		if(row != size-1) ans(row,row+1) = 1;
		ans(row,row) = -2;
	}
	return ans;
}

int main() {
	double step_t = 0.1;
	double step_x = 0.1;
	const double h_bar = 1;
	const double mass = 1;
	
	Eigen::VectorXcd space = linspace(step_x);
	Eigen::VectorXcd psi = map(sin,space);
	std::cout << psi << std::endl; //initial function
	
	int N_space = psi.rows(); //number of space points
	int N_time = (int)(1/step_t); //number of time points
	
	Eigen::MatrixXcd potential = Eigen::MatrixXcd::Zero(N_space, N_space); //to be implemented
	Eigen::MatrixXcd hamilton = -(h_bar)*(h_bar)/(2*mass) * second_derv(N_space) + potential; //hamiltonian
	
	Eigen::MatrixXcd factor_minus = Eigen::MatrixXcd::Identity(N_space,N_space) - (0.5/h_bar * dcplx_i * step_t) * hamilton; //RHS factor
	Eigen::MatrixXcd factor_plus = Eigen::MatrixXcd::Identity(N_space,N_space) + (0.5/h_bar * dcplx_i * step_t) * hamilton;	//coefficients matrix
	
	Eigen::VectorXcd RHS = factor_minus * psi;
	
	std::cout << factor_plus << std::endl << std::endl;
	std::cout << RHS << std::endl << std::endl;
	
	Eigen::VectorXcd next_psi = Eigen::VectorXcd(N_space);
	next_psi = factor_plus.fullPivLu().solve(RHS);
	std::cout << next_psi << std::endl << std::endl;
	/*
	Eigen::MatrixXcd multip = TPD*psi;
	std::cout << multip << std::endl;
	*/
	char c;
	std::cin >> c;
	return 0;
}