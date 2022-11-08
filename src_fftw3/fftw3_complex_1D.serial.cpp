//
// FFTW3 1D complex FFT sample program in C++
//

#include <complex>
//	#include <stdio.h>
#include <fftw3.h>
//	#include <iostream>
//	#include <complex.h>
using namespace std;

int main(int argc, char **argv)
{
const int nx_dim=8;
fftw_plan plan, pinv;

//	fftw_complex *indata, *outdata;
//	complex *indata, *outdata;
complex<double> *indata, *outdata;

int i;

//	indata = new complex<double> [nx_dim];
//	outdata = new complex<double> [nx_dim];

indata = (complex<double> *)fftw_alloc_complex(nx_dim);
outdata = (complex<double> *)fftw_alloc_complex(nx_dim);

// set the input data 
for (i=0; i<nx_dim; ++i) {
	indata[i] = complex<double>(i, -i);
}

printf("initial input data.\n");
for (i=0; i<nx_dim; ++i) {
	printf("%20.15f, %20.15f \n", real(indata[i]), imag(indata[i]));
	}

// create plans for out of place FFT
	plan = fftw_plan_dft_1d(nx_dim, reinterpret_cast<fftw_complex*>(indata), reinterpret_cast<fftw_complex*>(outdata), FFTW_FORWARD, FFTW_ESTIMATE);
	pinv = fftw_plan_dft_1d(nx_dim, reinterpret_cast<fftw_complex*>(outdata), reinterpret_cast<fftw_complex*>(indata), FFTW_BACKWARD, FFTW_ESTIMATE);

// forward transformation
fftw_execute(plan);

printf("After the FWD trans.\n");
for (i=0; i<nx_dim; ++i) {
	printf("%20.15f, %20.15f \n", real(outdata[i]), imag(outdata[i]));
	}

// backward synthesis
fftw_execute(pinv);

printf("After the BACK trans.\n");
for (i=0; i<nx_dim; ++i) {
	printf("%20.15f, %20.15f \n", real(indata[i]), imag(indata[i]));
	}

// normalization
for (i=0; i<nx_dim; ++i) {
	indata[i] /= (double)(nx_dim);
}

printf("After the normalization.\n");
for (i=0; i<nx_dim; ++i) {
	printf("%20.15f, %20.15f \n", real(indata[i]), imag(indata[i]));
	}

fftw_destroy_plan(plan);
fftw_destroy_plan(pinv);

}


