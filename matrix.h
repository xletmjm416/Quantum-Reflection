#ifndef MATRIX
#define MATRIX
#include <iostream>
#include <cstdlib>
#include <exception>

template<class T> //T - data type
class Matrix {
	/*	Matrix Class
	*	Basic matrix operations and solving linear equations.
	*/
	T *data; //entries of the matrix
	public:
	int rows, cols;
	
	//special memebers
	Matrix(int,int);
	~Matrix();
	//access
	T operator() (int row, int col) const; //retrieving elements
	void set_element(T, int, int); //setting elements
	
	//algebra
	Matrix<T> operator+ (const Matrix& other);
	Matrix<T> operator- ();
	Matrix<T> operator* (const Matrix& other);
}; //class Matrix

template<class T>
Matrix<T>::Matrix(int r, int c) : rows(r), cols(c) {
	/*	Constructs a linearised vrsion of a matrix	*/
	try {
		data = new T[rows*cols];
	}
	catch(std::exception& e) {
		throw;
	}
} //Matrix constructor

template<class T>
Matrix<T>::~Matrix() {
	/* Matrix destructor */
	delete this->data;
}

template<class T>
T Matrix<T>::operator() (int row, int col) const {
	/*	Retrieve element "(row,col)" from the matrix
		Indexing starts at 0!
	*/
	if(row > this->rows || col > this->cols || col < 0 || row < 0) {
		//reading outside of matrix
		throw std::exception();
		}
	return *(this->data+col+row*rows);
} //Matrix operator()

template<class T>
void Matrix<T>::set_element(T val, int row, int col) {
	/*Set element "(row,col)" to value "val" in the matrix.*/
	if(row > this->rows || col > this->cols || col < 0 || row < 0) throw std::exception();
	else {
		T *ptr = ((this->data)+col+row*rows);
		//elements are arranged for 3x3 matrix: 00 01 02 10 11 12 20 21 22
		*ptr = val;
		return;
	}
} //Matrix set_element

template<class T>
Matrix<T> Matrix<T>::operator+ (const Matrix<T>& other) {
	/* Matrix addition */
	if(rows != other.rows || cols != other.cols) {
		//dimensions mismatch
		throw std::exception();
	}
	
	int row, col;
	Matrix<T> ans(rows,cols);
	for(col = 0; col < cols; col++) {
		for(row = 0; row < rows; row++) {
			T summand1 = this->operator()(row, col);
			T summand2 = other(row,col);
			T val = summand1 + summand2;
			ans.set_element(val,row,col);
		}
	}
	return ans;
}

template<class T>
Matrix<T> Matrix<T>::operator- () {
	/* Negative of a matrix */
	Matrix<T> ans(rows, cols);
	int col, row;
	
	for(col = 0; col < cols; col++) {
		for(row = 0; row < rows; row++) {
			T val = -this->operator()(row, col);
			ans.set_element(val,row,col);
		}
	}
	return ans;
}

template<class T>
Matrix<T> Matrix<T>::operator* (const Matrix& other) {
	/* Matrix multiplication */
	if(cols != other.rows) {
		//dimensions mismatch
		throw std::exception();
	}
	
	Matrix<T> ans(rows, other.cols);
	int row, col, dummy;
	
	for(row = 0; row < rows; row++) {
		for(col = 0; col < other.cols; col++) {
			ans.set_element(0,row,col);
			for(dummy = 0; dummy < cols; dummy++) {
				T element1 = this->operator()(row, dummy);
				T element2 = other(dummy,col);
				T val = ans(row,col) + element1*element2;
				ans.set_element(val,row,col);
			}
		}
	}
	return ans;
}

#endif //header guard