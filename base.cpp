#include <iostream>
#include <cmath>
#include "Eigen/Dense"

double sinh(double x) {
	return std::sin(x);
}

Eigen::VectorXd linspace(double (*func)(double x), double N_steps, double system_size = 1) {
	/*	Initialises the wavefunction (system vector) from function callback "func"
	*	func		- function f(x) which defines the wavefunction at t=0;
	*	N_steps		- number of total steps
	*	system_size	- size of the whole system (in case that pointer arithmetic fails miserably)
	*/
	Eigen::VectorXd ans;
	ans = Eigen::VectorXd::LinSpaced(N_steps, 0, system_size);
	ans = ans.unaryExpr(std::ptr_fun(func)).eval();
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
	Eigen::MatrixXd A(3,3);
	Eigen::VectorXd B(3);
	A << 1,2,3,4,6,54,3,7,3;
	B << 1,2,3;
	std::cout << A*B << std::endl;
	
	Eigen::VectorXd psi = linspace(sin, 10);
	std::cout << psi << std::endl;
	
	Eigen::MatrixXd TPD = second_derv(10);
	std::cout << TPD << std::endl;
	
	Eigen::MatrixXd multip = TPD*psi;
	std::cout << multip << std::endl;
	char c;
	std::cin >> c;
	return 0;
}