#include "matrix.h"

int main() {
	Matrix<double> A = Matrix<double>(3,3);	
	A.set_element(1,0,0);
	A.set_element(7,0,1);
	A.set_element(3,0,2);
	A.set_element(6,1,0);
	A.set_element(-5,1,0);
	A.set_element(4,1,2);
	A.set_element(-1,2,0);
	A.set_element(2,2,1);
	A.set_element(5,2,2);

	Matrix<double> b = Matrix<double>(3,1);
	b.set_element(-1,0,0);
	b.set_element(5,1,0);
	b.set_element(2,2,0);
	
	
	Matrix<double> m3 = augment<double>(A,b);
	Matrix<double> m4 = gauss_elim<double>(A,b);
	
	A.show();
	std::cout << std::endl;
	b.show();
	std::cout << std::endl;
	m3.show();
	std::cout << std::endl;
	m4.show();
	
	//pause program
	char c;
	std::cin >> c;
	return 0;
} //main