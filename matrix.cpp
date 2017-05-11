#include "matrix.h"

int main() {
	/*Matrix<double> A = Matrix<double>(3,3);	
	A.set_element(1,0,0);
	A.set_element(0,0,1);
	A.set_element(1,0,2);
	A.set_element(1,1,0);
	A.set_element(1,1,1);
	A.set_element(1,1,2);
	A.set_element(1,2,0);
	A.set_element(-1,2,1);
	A.set_element(1,2,2);

	Matrix<double> b = Matrix<double>(3,1);
	b.set_element(1,0,0);
	b.set_element(2,1,0);
	b.set_element(1,2,0);
	*/
	Matrix<double> A = random<double>(2,2,10);	
	Matrix<double> b = random<double>(2,1,10);
	
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