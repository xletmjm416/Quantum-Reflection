#ifndef QMSYSTEM_H
#define QMSYSTEM_H

#include "common.h"

namespace utils {
	VectorR linspace(double step_x, double max);
	MatrixR three_pt_diff(int size);
	VectorC linear_solver(MatrixC coefficients, MatrixC rhs);
	dcplx gauss_state(dcplx, double,  double, double);
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
		QMSystem(VectorC, VectorR, double, double, double h_bar = 1, double mass = 1, int t_size = 1);
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
