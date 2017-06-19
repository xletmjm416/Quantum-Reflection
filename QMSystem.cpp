#include "QMSystem.h"

namespace utils {
	VectorR linspace(double step_x, double max) {
		/*	Initialises the linearly divided space.
		*	step_x		- length of a step
		*/
		VectorR ans;
		int N_steps = (int)(max/step_x);
		ans = VectorR::LinSpaced(N_steps, 0, max-step_x);
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
	
	VectorC thomas_solver(MatrixC coefficients, VectorC rhs) {
		// from wikipedia https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
		int i=0, N=rhs.rows();
		if(coefficients.rows() != N) throw;
		dcplx temp;
		VectorC ans(N);
		
		//first step
		i = 0;
		coefficients(i,i+1) = coefficients(i,i+1)/coefficients(i,i); //c* = c/b
		rhs(i) = rhs(i)/coefficients(i,i); //d* = d/b
		
		// forward propagation
		for(i=1;i<(N-1);i++) {
			coefficients(i,i+1) = coefficients(i,i+1)/(coefficients(i,i) - coefficients(i,i-1) * coefficients(i-1,i)); //c* = c/(b - a . c_(i-1)* )
			
			rhs(i) = (rhs(i) - coefficients(i,i-1)*rhs(i-1)) / (coefficients(i,i) - coefficients(i,i-1)*coefficients(i-1,i)); // d* = [ d - a * d_(i-1)* ] / [ b - a * c_(i-1)* ]
		}
		i = N-1;
		rhs(i) = (rhs(i)-coefficients(i,i-1)*rhs(i-1))/(coefficients(i,i)-coefficients(i,i-1)*coefficients(i-1,i));
		
		//back substitution
		ans(i) = rhs(i); //i == N-1 at this point
		for(i=(N-2);i>=0;i--) {
			ans(i) = rhs(i) - coefficients(i,i+1)*ans(i+1); //x = d* - c* . x_i+1
		}
		
		return ans;
	}
	
	dcplx gauss_function(dcplx x, double k,  double mean, double spread) {
		dcplx exponent_1(dcplx_i * k * x ); //ikx
		dcplx exponent_2 = -1.0*(x-mean)*(x-mean)/(4*spread*spread); //-ax^2
		dcplx constant=std::pow(1/(2*PI_CONST*spread*spread),0.25); //constant of normalisation = 4th_root(1/(2pi sigma^2))
		dcplx gauss = constant * std::exp(exponent_1+exponent_2);
		return gauss;
	}
	
	dcplx linear_function(dcplx x, double stretch,  double translation, double offset) {
		dcplx ans = stretch*(x-offset)+translation;
		return ans;
	}
	
	dcplx heaviside_function(dcplx x, double centre, double height, short orient) {
		if(orient == 0) {
			if(x.real() <= centre) return dcplx(height,0);
			if(x.real() > centre) return dcplx(0,0);
		}
		else {
			if(x.real() >= centre) return dcplx(height,0);
			if(x.real() < centre) return dcplx(0,0);
		}
		return dcplx(0,0);
	}
	
	dcplx step_function(dcplx x, double pos, double width) {
		double height = 1/width;
		if(std::abs(x) > pos+width/2 || std::abs(x) < pos-width/2) return 0;
		else return height;
	}
		
} //namespace utils

QMSystem::QMSystem(VectorC initial, VectorR potential, double step_x, double step_t,
					double h_bar, double mass, double x_size, double t_size) : 
					h_bar(h_bar), mass(mass), step_x(step_x), step_t(step_t),
					x_size(x_size), t_size(t_size) {
			this->N_space = initial.rows();
			this->N_time = (int)t_size/step_t;
			this->wavefunction = initial;
			if(potential.rows() != initial.rows()) throw; //dimensions mismatch
			this->potential = potential;
			this->hamilton = hamiltonian(potential);
}

MatrixC QMSystem::hamiltonian(VectorR pot) {
	/*	Updates the hamiltonian of the system	*/
	double constant = -(h_bar)*(h_bar)/(2*mass*(step_x*step_x));
	MatrixC ans = constant * utils::three_pt_diff(N_space); 	// second derivative part
	ans += pot.asDiagonal(); 					// add potential
	OUTPUT("Peek on Hamiltonian: ", std::endl << ans.block(0,0,5,5));
	return ans;
}

VectorC QMSystem::cranknicolson() {
	MatrixC identity = MatrixC::Identity(N_space,N_space);

	MatrixC factor_minus = identity - (0.5/h_bar * dcplx_i * step_t) * hamilton; //RHS factor
	MatrixC factor_plus = identity + (0.5/h_bar * dcplx_i * step_t) * hamilton;	//coefficients matrix

	VectorC ans(N_space), rhs(N_space);
	rhs = factor_minus*wavefunction;
	ans = utils::thomas_solver(factor_plus, rhs); //solve implicit scheme

	this->wavefunction = ans;
	return ans;
}