#ifndef QMSYSTEM_H
#define QMSYSTEM_H

#include "common.h"

namespace utils {
	VectorR linspace(double step_x, double max);
	MatrixR three_pt_diff(int size);
	VectorC linear_solver(MatrixC coefficients, MatrixC rhs);
	dcplx gauss_function(dcplx, double,  double, double);
	VectorC thomas_solver(MatrixC, VectorC);
	dcplx linear_function(dcplx, double,  double, double);
	dcplx heaviside_function(dcplx, double, double, short);
	dcplx step_function(dcplx, double, double);
}

class QMSystem {
	/* Whole QM System with the following parameters:
	space
	hamiltonian
	step x
	step t
	number of steps x
	number of steps t
	h bar
	mass
	*/
	public:
		QMSystem(VectorC, VectorR, double, double, double h_bar, double mass, double x_size, double t_size);
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
		double x_size, t_size;
};

#endif
