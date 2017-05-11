#ifndef MATRIX
#define MATRIX
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <exception>

template<class T> //T - data type
class Matrix {
	/*	Matrix Class
	*	Basic matrix operations and solving linear equations.
	*/
	public:
	
	//special memebers
	Matrix(int,int);
	~Matrix();
	
	//access
	T operator() (int row, int col) const; //retrieving elements
	void set_element(T, int, int); //setting elements
	void set_element(T, int); //linear indexing
	int get_rows() const;
	int get_cols() const;
	
	//algebra
	Matrix<T> operator+ (const Matrix& other) const;
	Matrix<T> operator- () const;
	Matrix<T> operator* (const Matrix& other) const;
	
	//other
	void show();
	
	private:
	T *data; //entries of the matrix
	const int rows, cols;
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
	delete [] this->data;
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
	return *(this->data+col+row*cols);
} //Matrix operator()

template<class T>
void Matrix<T>::set_element(T val, int row, int col) {
	/*Set element "(row,col)" to value "val" in the matrix.*/
	if(row > this->rows || col > this->cols || col < 0 || row < 0) throw std::exception();
	else {
		T *ptr = &( (this->data)[col+row*(this->cols)] );
		//elements are arranged for 3x3 matrix: 00 01 02 10 11 12 20 21 22
		*ptr = val;
		return;
	}
} //Matrix set_element

template<class T>
void Matrix<T>::set_element(T val, int idx) {
	/*Set element "(idx)" to value "val" in the matrix when it is treated as a linear collection of numbers.*/
	if(idx > rows*cols || idx < 0) throw std::exception();
	else {
		T *ptr = ((this->data)+idx);
		*ptr = val;
		return;
	}
} //Matrix set_element

template<class T>
int Matrix<T>::get_rows() const {
	return this->rows;
}

template<class T>
int Matrix<T>::get_cols() const {
	return this->cols;
}

template<class T>
Matrix<T> Matrix<T>::operator+ (const Matrix<T>& other) const {
	/* Matrix addition.
	Direct computation by definition.
	Operation cost: O(n^2). (for A_(nk)+B_(nk), Cost=O(nk))
	*/
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
} //Matrix addition

template<class T>
void Matrix<T>::show() {
	/* Print the matrix on the console
	*/
	for(int row = 0; row < this->get_rows(); row++) {
		for(int col = 0; col < this->get_cols(); col++) {
			std::cout << this->operator()(row,col) << " ";
		}
		std::cout << std::endl;
	}
} //matrix show

template<class T>
Matrix<T> Matrix<T>::operator- () const {
	/* Negative of a matrix.
	Direct computation by definition.
	Operation cost: O(n^2). (for A_(nk), Cost=O(nk))
	*/
	Matrix<T> ans(rows, cols);
	int col, row;
	
	for(col = 0; col < cols; col++) {
		for(row = 0; row < rows; row++) {
			T val = -this->operator()(row, col);
			ans.set_element(val,row,col);
		}
	}
	return ans;
} //Matrix negative

template<class T>
Matrix<T> Matrix<T>::operator* (const Matrix& other) const {
	/* Matrix multiplication.
	Direct computation by definition.
	Operation cost: O(n^3). (for A_(nk)*B_(kj), Cost=O(nkj))
	*/
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
} //Matrix multiplication

template<class T>
Matrix<T> identity(int n, T val = 1) {
	/*Creates an identity matrix, or, if argument val is passed, a matrix with val on diagonal.*/
	Matrix<T> ans(n,n);
	int row, col;
	for(col = 0; col < n; col++) {
		for(row = 0; row < n; row++) {
			if (row == col) ans.set_element(val,row,col);
			else ans.set_element(0,row,col);
		}
	}
	return ans;
}

template<class T>
Matrix<T> augment(const Matrix<T> &A, const Matrix<T> &B) {
	/*	Produces augmented matrix (A|B) */
	if(A.get_rows() != B.get_rows()) throw std::exception();
	int rows = A.get_rows();
	int colsA = A.get_cols();
	int colsB = B.get_cols();
	T val; //dummy vars
	
	Matrix<T> *aug_p = new Matrix<T>(rows,colsA+colsB); //augmented matrix
	
	for(int row=0; row < rows; row++) {
		for(int col=0; col < colsA; col++) {
			val = A(row,col);
			aug_p->set_element(val,row,col);
		}
		for(int col=0; col < colsB; col++) {
			val = B(row,col);
			aug_p->set_element(val,row,colsA+col);
		}
	}	
	return *aug_p;
}	//augment

template<class T>
Matrix<T> gauss_elim(const Matrix<T> &coeffs, const Matrix<T> &rhs) {
	/*	Solve for x in coeffs*x=rhs.
	*/
	if(coeffs.get_rows() != coeffs.get_cols() || rhs.get_cols() != 1 || rhs.get_rows() != coeffs.get_rows()) throw std::exception();
	
	Matrix<T> aug = augment(coeffs, rhs);
	int rows = aug.get_rows();
	int cols = aug.get_cols();

	//forward elimination
	
	#define PEEK(info, x) std::cout << info << x << std::endl;
	T pivot, factor, val;
	int piv_row = 0; //pivot's row
	for(int sweep=0; sweep < cols; sweep++) {
		//sweeps the row searching for pivot (leftmost non-zero term);
		val = aug(piv_row,sweep);
		if(val==0) continue;
		else { pivot = val; pivot_col=sweep; break; }
	}
	if(pivot == 0) throw std::exception();
	PEEK("pivot ", pivot)
	for(int fwd=0; fwd < rows; fwd++) {
		PEEK("fwd ", fwd)
		T factor = aug(fwd, piv_col) / pivot; //TODO if not zero
		PEEK("factor ", factor)
		for(int col=0; col < cols; col++) {
			//update consecutive rows column by column
			//eliminate with the pivot
			if(fwd == piv_row) {
			//normalise the pivot's row
				val = aug(fwd, col) / pivot;
				aug.set_element(val, fwd, col);
			}
			else {
				val = aug(fwd, col) - factor * aug(piv_row,col);
				PEEK("val ", val)
				aug.set_element(val, fwd, col);
			}
		}
	}

	return aug;
} //gauss elim

#endif //header guard