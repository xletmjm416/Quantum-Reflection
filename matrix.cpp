#include "matrix.h"

int main() {
	Matrix<int> m1 = Matrix<int>(2,2);	
	m1.set_element(2,0,0);
	m1.set_element(0,0,1);
	m1.set_element(0,1,0);
	m1.set_element(1,1,1);
	
	
	Matrix<int> m2 = Matrix<int>(2,2);
	m2.set_element(5,0,0);
	m2.set_element(4,0,1);
	m2.set_element(-1,1,0);
	m2.set_element(2,1,1);

	Matrix<int> b = Matrix<int>(2,1);
	b.set_element(-1,0,0);
	b.set_element(5,1,0);
	
	Matrix<int> m3 = augment<int>(m1,m2);

	m1.show();
	std::cout << std::endl;
	m2.show();
	std::cout << std::endl;
	m3.show();
	std::cout << std::endl;
	
	std::cout << m3(0,0) << m3(0,1) << m3(0,2) << std::endl << m3(1,0) << m3(1,1) << m3(1,2) << std::endl;
	//pause program
	char c;
	std::cin >> c;
	return 0;
} //main