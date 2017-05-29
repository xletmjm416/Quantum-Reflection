#ifndef FFT_H
#define FFT_H
#include "common.h"

namespace FFT {
	unsigned int bitrev(unsigned int, unsigned int);
	VectorC collect(VectorC, dcplx, int);
	dcplx kth_elem(VectorC, int);
	VectorC FFT(VectorC);
}

#endif