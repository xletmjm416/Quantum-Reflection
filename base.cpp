#include <iostream>
#include <cmath>
#include "Eigen/Dense"

double func(double x) {
	return std::sinh(x);
}

Eigen::VectorXd wavefunction(double (*func)(double x), double N_steps, double system_size = 1) {
	/*	Initialises the wavefunction (system vector) from function callback "initial"
	*	func		- function f(x) which defines the wavefunction at t=0;
	*	N_steps		- number of total steps
	*	system_size	- size of the whole system (in case that pointer arithmetic fails miserably)
	*/
	Eigen::VectorXd ans;
	ans = Eigen::VectorXd::LinSpaced(N_steps, 0, system_size);
	ans = ans.unaryExpr(std::ptr_fun(func)).eval();
	return ans;
}

int main()
{
	Eigen::MatrixXd A(3,3);
	Eigen::VectorXd B(3);
	A << 1,2,3,4,6,54,3,7,3;
	B << 1,2,3;
	
	Eigen::VectorXd psi = wavefunction(func, 10);
	
	std::cout << A*B;
	std::cout << psi;
	char c;
	std::cin >> c;
	return 0;
}