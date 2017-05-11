#include "matrix.h"

int main() {
	Matrix<double> m1 = Matrix<double>(2,2);
	Matrix<double> m2 = Matrix<double>(2,2);
	
	m1.set_element(1,0,0);
	m1.set_element(0,0,1);
	m1.set_element(0,1,0);
	m1.set_element(1,1,1);
	
	m2.set_element(1,0,0);
	m2.set_element(0,0,1);
	m2.set_element(-1,1,0);
	m2.set_element(2,1,1);

	Matrix<double> m3 = m1 + (-m2);
	Matrix<double> m4 = m2*m1;
	
	std::cout << m1(0,0) << " " << m1(0,1) << std::endl << m1(1,0) << " " << m1(1,1) << std::endl << std::endl;
	std::cout << m2(0,0) << " " << m2(0,1) << std::endl << m2(1,0) << " " << m2(1,1) << std::endl << std::endl;
	std::cout << m3(0,0) << " " << m3(0,1) << std::endl << m3(1,0) << " " << m3(1,1) << std::endl << std::endl;
	std::cout << m4(0,0) << " " << m4(0,1) << std::endl << m4(1,0) << " " << m4(1,1) << std::endl << std::endl;
	
	Matrix<double> m5 = identity<double>(2);
	std::cout << m5(0,0) << " " << m5(0,1) << std::endl << m5(1,0) << " " << m5(1,1) << std::endl << std::endl;
	//pause program
	char c;
	std::cin >> c;
	return 0;
} //main