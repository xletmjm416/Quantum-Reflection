#include "FFT.h"

namespace FFT {
	// From Numerical Recipes in C++, 2nd ed. p. 513
	
	unsigned int bitrev(unsigned int n, unsigned int bits) {
		// from http://www.katjaas.nl/bitreversal/bitreversal.html
		unsigned int nrev, N;
		unsigned int count;   
		N = 1<<bits;
		count = bits-1;   // initialize the count variable
		nrev = n;
		for(n>>=1; n; n>>=1)
		{
			nrev <<= 1;
			nrev |= n & 1;
			count--;
		}

		nrev <<= count;
		nrev &= N - 1;

		return nrev;
	} //bitrev

	VectorC collect(VectorC sample, dcplx root, int k) {
		//one step in the Danielson-Lanczos formula - reduces sample to half its size
		int N = sample.rows();
		VectorC ans(N/2);
		for(int i=0; i<N; i+=2) {
			ans(i/2) = sample(i) + std::pow(root,(double)k) * sample(i+1);
		}
		return ans;
	}
	
	dcplx kth_elem(VectorC sample, int k) {
		int N = sample.rows();
		VectorC ans(N);
		
		dcplx root = std::exp(dcplx_i*2.0*PI_CONST/(double)N);
		//obtain k-th elementh
		VectorC reduced = collect(sample, root, k);
		while(reduced.rows() != 1) {
			reduced = collect(reduced, root, k);
		}
		return reduced(0);
	} //kth-elem
	
	VectorC FFT(VectorC sample) {
		/* sample is of size N = 2^n */
		
		int N = sample.rows();
		VectorC rev(N);
		rev = sample;
		if((N & (N-1)) != 0) throw; //sample size is not a power of two
		int power = (int)(std::log(N)/std::log(2));
		
		//shuffle to get bit-reversed order
		for(int i = 0; i<N; i++) {
			int reversed = bitrev(i,power);	
			if(i != reversed){
				rev(i) = sample(reversed);
				rev(reversed) = sample(i);
			}
		}
		
		VectorC ans(N);
		for(int j=0; j<N; j++) {
			ans(j) = kth_elem(rev,j);
		}
		return ans;
	} //FFT

	VectorC DITFFT(VectorC sample) {
		int N = sample.rows();
		VectorC ans(N);
		ans = sample;
		
		if((N & (N-1)) != 0) throw; //sample size is not a power of two
		int power = (int)(std::log(N)/std::log(2));
		
		//shuffle to get bit-reversed order
		for(int i = 0; i<N; i++) {
			int reversed = bitrev(i,power);	
			if(i != reversed){
				ans(i) = sample(reversed);
				ans(reversed) = sample(i);
			}
		}
		
		// DIT FFT
		int step = 0;
		int elems = 1;
		
		for (int i=0; i < power; i++) {
			step = elems;
			elems *=2;	//increase number of elements taken every step
			
			for (int j=0; j < step; j++) {
				double theta = 2*PI_CONST/elems*j;
				dcplx twiddle = std::exp(theta); //twiddle factor
				
				for (int k=j; k < N; k += elems) {
					//butterfly
					dcplx temp = twiddle * ans(k+step);
					
					ans(k+step) = ans(k) - temp;					
					ans(k) = ans(k) + temp;
				}
			}
		}
											  
		return ans;
	}
} //namespace FFT
