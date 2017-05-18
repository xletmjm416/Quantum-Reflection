#include "QMSystem.h"

namespace utils {
	VectorR linspace(double step_x, int start, int end) {
		/*	Initialises the linearly divided space.
		*	step_x		- length of a step
		*/
		VectorR ans;
		int N_steps = (int)((end-start)/step_x);
		ans = VectorR::LinSpaced(N_steps, start, end);
		return ans;
	}
	
	MatrixR three_pt_diff(int size) {
		/*	Produces a three-point difference algorithm matrix (tridiagonal with 1, -2, 1 values).
		*	size		- size of the matrix (size x size)
		*/
		MatrixR ans(size, size);
		ans = MatrixR::Zero(size, size);
		for(int row = 0; row<size; row++) {
			if(row != 0) ans(row,row-1) = 1;
			if(row != size-1) ans(row,row+1) = 1;
			ans(row,row) = -2;
		}
		return ans;
	}
	
	VectorC linear_solver(MatrixC coefficients, MatrixC rhs) {
		/*	Solves for vector x such that A*x = B*c, where
		coefficients		- matrix A
		rhs 				- vector c
		*/
		VectorC ans(rhs.rows());
		ans = coefficients.fullPivLu().solve(rhs); // use full pivoted LU decomposition
		return ans;
	}

}

QMSystem::QMSystem(VectorC initial, VectorR potential, double step_x, double step_t,
					double h_bar, double mass, int x_size, 
					int t_size) : 
					h_bar(h_bar), mass(mass), step_x(step_x), step_t(step_t), x_size(x_size), t_size(t_size) {
			this->N_space = (int)x_size/step_x;
			this->N_time = (int)t_size/step_t;
			this->wavefunction = initial;
			this->potential = potential;
			this->hamilton = hamiltonian(potential);
}

MatrixC QMSystem::hamiltonian(VectorR pot) {
	/*	Updates the hamiltonian of the system	*/
	double constant = -(h_bar)*(h_bar)/(2*mass*(step_x*step_x));
	MatrixC ans = constant * utils::three_pt_diff(N_space); 	// second derivative part
	ans += pot.asDiagonal(); 					// add potential
	return ans;
}

VectorC QMSystem::cranknicolson() {
	MatrixC identity = MatrixC::Identity(N_space,N_space);

	MatrixC factor_minus = identity - (0.5/h_bar * dcplx_i * step_t) * hamilton; //RHS factor
	MatrixC factor_plus = identity + (0.5/h_bar * dcplx_i * step_t) * hamilton;	//coefficients matrix

	VectorC ans(N_space);
	ans = utils::linear_solver(factor_plus, factor_minus*wavefunction); //solve implicit scheme
	this->wavefunction = ans;
	return ans;
}

dcplx sin_state(dcplx x) {
	//if(std::abs(x)>0.9) return dcplx(1e100,0);
	if(std::abs(x)>0.4 && std::abs(x)<0.6) return std::sin(PI_CONST/0.2*(x-0.4));
	else return 0;
}

dcplx gauss_state(dcplx x) {
	double a = 60;
	double mean = 0.5;
	double k = 10;
	dcplx exponent = dcplx_i*k*x - 1*a*std::pow(x - mean,2);
	dcplx constant = std::pow(2*a/PI_CONST,0.25);
	dcplx gauss = constant * std::exp(exponent);
	return gauss;
}