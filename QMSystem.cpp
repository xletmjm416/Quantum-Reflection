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
	
	dcplx gauss_state(dcplx x, double k,  double mean, double spread) {
		double a = 1/(4*spread*spread); //exp(-kx - a x^2)
		dcplx exponent = dcplx_i * k * (x-mean) - a*(x-mean)*(x-mean);
		dcplx constant(std::pow(2*a/PI_CONST,0.25),0); //constant of normalisation = 4th_root(2a/pi)
		dcplx gauss = constant * std::exp(exponent);
		return gauss;
	}
		
} //namespace utils

QMSystem::QMSystem(VectorC initial, VectorR potential, double step_x, double step_t,
					double h_bar, double mass) : 
					h_bar(h_bar), mass(mass), step_x(step_x), step_t(step_t) {
			this->N_space = initial.rows();
			this->N_time = (int)TIME_SIZE/step_t;
			this->wavefunction = initial;
			if(potential.rows() != initial.rows()) throw;
			this->potential = potential;
			this->hamilton = hamiltonian(potential);
}

MatrixC QMSystem::hamiltonian(VectorR pot) {
	/*	Updates the hamiltonian of the system	*/
	double constant = -(h_bar)*(h_bar)/(2*mass*(step_x*step_x));
	MatrixC ans = constant * utils::three_pt_diff(N_space); 	// second derivative part
	std::cout << ans.rows() << " " << ans.cols() << " " << pot.rows() << std::endl;
	ans += pot.asDiagonal(); 					// add potential
	OUTPUT("Hamiltonian: ", std::endl << ans.block(0,0,3,3));
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