//
// FFTW3 1D complex FFT sample program
//

//	#include <stdio.h>
#include <complex.h>
#include <fftw3.h>

int main(int argc, char **argv)
{
const ptrdiff_t nx_dim=8;
int debug_print=1;

fftw_plan plan, pinv;
fftw_complex *indata, *outdata;

//	double complex imaginary(0,-1);

int i;
printf("Testing C driver.\n");

indata = fftw_alloc_complex(nx_dim);
outdata = fftw_alloc_complex(nx_dim);

// set the input data 
for (i=0; i<nx_dim; ++i) {
	indata[i] = i - i*I;
}

printf("initial data.\n");
for (i=0; i<nx_dim; ++i) {
	printf("%20.15f, %20.15f \n", creal(indata[i]), cimag(indata[i]));
	}

// create plans for out of place FFT
	plan = fftw_plan_dft_1d(nx_dim, indata, outdata, FFTW_FORWARD, FFTW_ESTIMATE);
	pinv = fftw_plan_dft_1d(nx_dim, outdata, indata, FFTW_BACKWARD, FFTW_ESTIMATE);

// forward transformation
fftw_execute(plan);

printf("After the FWD trans.\n");
for (i=0; i<nx_dim; ++i) {
	printf("%20.15f, %20.15f \n", creal(outdata[i]), cimag(outdata[i]));
	}

// backward synthesis
fftw_execute(pinv);

printf("After the BACK trans.\n");
for (i=0; i<nx_dim; ++i) {
	printf("%20.15f, %20.15f \n", creal(indata[i]), cimag(indata[i]));
	}

// normalization
for (i=0; i<nx_dim; ++i) {
	indata[i] /= (double)(nx_dim);
}

printf("After the normalization.\n");
for (i=0; i<nx_dim; ++i) {
	printf("%20.15f, %20.15f \n", creal(indata[i]), cimag(indata[i]));
	}

fftw_destroy_plan(plan);
fftw_destroy_plan(pinv);

}


#include        <unistd.h>
#include        <sys/time.h>
#include        <sys/resource.h>
double second_()
{
  struct timeval s_val;
  gettimeofday(&s_val,0);
  return ((double) s_val.tv_sec + 0.000001*s_val.tv_usec);
}
//	double dptime_timef()
//	{
//		return timef_();
//	}
