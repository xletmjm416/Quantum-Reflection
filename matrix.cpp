#include <iostream>
#include <cstdlib>

template<class T> //T - data type
class Matrix {
	/*	Matrix Class
	*	Basic matrix operations and solving linear equations.
	*/
	T *data; //entries of the matrix
	public:
	int rows, cols;
	
	Matrix(int,int);
	Matrix<T> operator+ (const Matrix& other);
	Matrix<T> operator- ();
	Matrix<T> operator* (const Matrix& other);
	T operator() (int row, int col); //retrieving elements
	bool set_element(T, int, int); //setting elements
};


template<class T>
Matrix<T>::Matrix(int r, int c) : rows(r), cols(c) {
	/*	Constructs a linearised vrsion of a matrix	*/
	data = new T[rows*cols];
	if(data==NULL) {
		std::cerr << "Error with matrix dynamic memory allocation." << std::endl;
		std::exit(1); //exit with error code 1
	}
}

template<class T>
T Matrix<T>::operator() (int row, int col) {
	/*	Retrieve element "(row,col)" from the matrix
		Indexing starts at 0!
	*/
	if(row > this->rows || col > this->cols || col < 0 || row < 0) {
		std::cerr << "Warning: reading out of bound of a matrix" << std::endl;
		return 0;
		}
	return *(this->data+col+row*rows);
}

template<class T>
bool Matrix<T>::set_element(T val, int row, int col) {
	/*Set element "(row,col)" to value "val" in the matrix.*/
	if(row > this->rows || col > this->cols || col < 0 || row < 0) return 1;
	else {
		T *ptr = ((this->data)+col+row*rows);
		*ptr = val;
		return 0;
	}
}

int main() {
	Matrix<double> m1 = Matrix<double>(3,3);
	m1.set_element(2.3,1,1);
	std::cout << m1(1,1) << std::endl;
	std::cout << m1(0,3) << std::endl;
	//pause program
	char c;
	std::cin >> c;
	return 0;
}