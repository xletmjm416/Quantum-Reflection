#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>
#include "Eigen/Dense"

typedef std::complex<double> dcplx;
const dcplx dcplx_i(0.0,1.0);

dcplx sin_state(dcplx x) {
	//if(std::abs(x)>0.9) return dcplx(1e100,0);
	if(std::abs(x)>0.4 && std::abs(x)<0.6) return std::sin(3.14159/0.2*(x-0.4));
	else return 0;
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
	/*	Produces a three-point difference algorithm matrix (tridiagonal with 1, -2, 1 values).
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

Eigen::VectorXcd linear_solver(Eigen::MatrixXcd coefficients, Eigen::MatrixXcd correction, Eigen::VectorXcd rhs) {
	/*	Solves for vector x such that A*x = B*c, where
	coefficients		- matrix A
	correction			- matrix B
	rhs 				- vector c
	*/
	Eigen::VectorXcd ans(rhs.rows());
	ans = coefficients.fullPivLu().solve(correction*rhs); // use full pivoted LU decomposition
	return ans;
}

Eigen::VectorXcd step(Eigen::VectorXcd initial, double step_x, double step_t, Eigen::VectorXcd pot) {
	//make sure rows of potental match
	const double h_bar = 1;
	const double mass = 1;
	
	int N_space = initial.rows(); //number of space points
	//int N_time = (int)(1/step_t); //number of time points

	Eigen::MatrixXcd identity = Eigen::MatrixXcd::Identity(N_space,N_space);
	//Eigen::MatrixXcd potential = pot.asDiagonal();
	Eigen::MatrixXcd hamilton = -(h_bar)*(h_bar)/(2*mass*(step_x*step_x)) * second_derv(N_space); //hamiltonian - derivative part

	Eigen::MatrixXcd factor_minus = identity - (0.5/h_bar * dcplx_i * step_t) * hamilton; //RHS factor
	Eigen::MatrixXcd factor_plus = identity + (0.5/h_bar * dcplx_i * step_t) * hamilton;	//coefficients matrix

	Eigen::VectorXcd ans(N_space);
	ans = linear_solver(factor_plus, factor_minus, initial); //solve implicit scheme
	return ans;
}

int main() {
	double step_t = 0.08;
	double step_x = 0.01;
	int N_space = (int)(1/step_x); //number of space points
	int N_time = (int)(1/step_t); //number of time points
	
	Eigen::VectorXcd space = linspace(step_x);
	Eigen::VectorXcd psi = space.unaryExpr(&sin_state);
	
	Eigen::VectorXcd next_psi = Eigen::VectorXcd(N_space);
	Eigen::VectorXd prob_distr = Eigen::VectorXd(N_space);
	Eigen::VectorXcd pot = Eigen::VectorXcd::Zero(N_space);
	
	std::ofstream output;
	output.open("out.csv");
	
	prob_distr = psi.cwiseAbs2();
	output << prob_distr.transpose() << std::endl;
	for(int t=0; t<N_time*10;t++) {
		next_psi = step(psi, step_x, step_t, pot);
		prob_distr = next_psi.cwiseAbs2();
		//std::cout << prob_distr.transpose() << std::endl;
		output << prob_distr.transpose() << std::endl;
		psi = next_psi;
	}
	std::cout << "finished";
	output.close();
	char c;
	std::cin >> c;
	return 0;
}