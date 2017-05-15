#include <iostream>
#include <cmath>
#include <complex>
#include "Eigen/Dense"

typedef std::complex<double> dcplx;

dcplx sin(dcplx x) {
	return std::sin(x);
}

Eigen::VectorXd linspace(double N_steps) {
	/*	Initialises the linearly divided space.
	*	N_steps		- number of total steps
	*	system_size	- size of the whole system (in case that pointer arithmetic fails miserably)
	*/
	Eigen::VectorXd ans;
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

int main()
{
	Eigen::VectorXcd space = linspace(10);
	Eigen::VectorXcd psi = map(sin,space);
	std::cout << psi << std::endl;
	
	Eigen::MatrixXd TPD = second_derv(10);
	std::cout << TPD << std::endl;
	
	Eigen::MatrixXcd multip = TPD*psi;
	std::cout << multip << std::endl;
	char c;
	std::cin >> c;
	return 0;
}