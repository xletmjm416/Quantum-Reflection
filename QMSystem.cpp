#include "QMSystem.h"

namespace utils {
	VectorR linspace(double step_x, double max) {
		/*	Initialises the linearly divided space.
		*	step_x		- length of a step
		*/
		VectorR ans;
		int N_steps = (int)(max/step_x);
		ans = VectorR::LinSpaced(N_steps, 0, max);
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
	
	dcplx gauss_state(dcplx x, double k,  double mean, double stddev) {
		dcplx exponent = dcplx_i * k * x - 0.5 * std::pow( (x - mean) / stddev , 2 );
		dcplx constant = 1/(stddev*std::pow(2*PI_CONST,0.5));
		dcplx gauss = constant * std::exp(exponent);
		return gauss;
	}

}

QMSystem::QMSystem(VectorC initial, VectorR potential, double step_x, double step_t,
					double h_bar, double mass, int t_size) : 
					h_bar(h_bar), mass(mass), step_x(step_x), step_t(step_t), t_size(t_size) {
			this->N_space = initial.rows();
			this->x_size = N_space*step_x;
			this->N_time = (int)t_size/step_t;
			this->wavefunction = initial;
			if(potential.rows() != initial.rows())throw;
			this->potential = potential;
			this->hamilton = hamiltonian(potential);
}

MatrixC QMSystem::hamiltonian(VectorR pot) {
	/*	Updates the hamiltonian of the system	*/
	double constant = -(h_bar)*(h_bar)/(2*mass*(step_x*step_x));
	MatrixC ans = constant * utils::three_pt_diff(N_space); 	// second derivative part
	std::cout << ans.rows() << " " << ans.cols() << " " << pot.rows() << std::endl;
	ans += pot.asDiagonal(); 					// add potential
	return ans;
}

VectorC QMSystem::cranknicolson() {
	MatrixC identity = MatrixC::Identity(N_space,N_space);

	MatrixC factor_minus = identity - (0.5/h_bar * dcplx_i * step_t) * hamilton; //RHS factor
	MatrixC factor_plus = identity + (0.5/h_bar * dcplx_i * step_t) * hamilton;	//coefficients matrix

	VectorC ans(N_space), rhs(N_space);
	rhs = factor_minus*wavefunction;
	ans = utils::linear_solver(factor_plus, rhs); //solve implicit scheme

	this->wavefunction = ans;
	return ans;
}

dcplx sin_state(dcplx x) {
	//if(std::abs(x)>0.9) return dcplx(1e100,0);
	if(std::abs(x)>0.4 && std::abs(x)<0.6) return std::sin(PI_CONST/0.2*(x-0.4));
	else return 0;
}
/*
dcplx gauss_state(dcplx x) {
	double mean = 5;
	double stddev = 1;
	double k = 10000;
	dcplx exponent = dcplx_i * k * x - 0.5 * std::pow( (x - mean) / stddev , 2 );
	dcplx constant = 1/(stddev*std::pow(2*PI_CONST,0.5));
	dcplx gauss = constant * std::exp(exponent);
	return gauss;
}
*/