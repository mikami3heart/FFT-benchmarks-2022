// FFTW3 3D complex FFT sample program
//
// The original data is 1D array[nx_dim * ny_dim * nz_dim]
// 	equivalent to 3D array array[nx_dim][ny_dim][nz_dim]
//

#include <stdio.h>
#include <complex.h>
#include <fftw3-mpi.h>
#include <mpi.h>
#include <omp.h>
double second_(void);
double dptime_timef(void);
double dptime_mpi(void);
int threads_ok;


int main(int argc, char **argv)
{
//	const ptrdiff_t nx_dim=4, ny_dim=2, nz_dim=4;
//	const ptrdiff_t nx_dim=64, ny_dim=64, nz_dim=64;
//	const ptrdiff_t nx_dim=128, ny_dim=128, nz_dim=128;
const ptrdiff_t nx_dim=256, ny_dim=256, nz_dim=256;
//	const ptrdiff_t nx_dim=512, ny_dim=512, nz_dim=512;
int debug_print=0;

fftw_plan plan, pinv;
fftw_complex *data;
ptrdiff_t alloc_local, my_dim, my_dim_start;

int i,j,k,ip;
int npes, my_id, nthreads=1;
double w1,w2,w3,w4,w5,w6,w7,w8;

// if using both MPI x Threads: p66 of FFTW user manual version 3.3.2, 28 April 2012
int provided;
//	MPI_Init(&argc, &argv);
MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
MPI_Comm_size(MPI_COMM_WORLD, &npes);

threads_ok = provided >= MPI_THREAD_FUNNELED;
if (threads_ok) {
	if (debug_print == 2) {
		printf("MPI_Init_thread() my_id=%d thread request accepted\n", my_id);
	}
	threads_ok = fftw_init_threads();
} else {
	printf("*** MPI_Init_thread() my_id=%d thread request failed.\n", my_id);
	MPI_Finalize(); return(0);
}

fftw_mpi_init();

nthreads=omp_get_max_threads();
fftw_plan_with_nthreads(nthreads);

if (debug_print == 2 && my_id == 0) {
	printf("3D complex FFT(%d,%d,%d) nthreads=%d\n", nx_dim,ny_dim,nz_dim, nthreads);
}


// set local size 
alloc_local = fftw_mpi_local_size_3d(nx_dim, ny_dim, nz_dim,  MPI_COMM_WORLD,
	&my_dim, &my_dim_start);

if (debug_print == 2) { MPI_Barrier(MPI_COMM_WORLD);
	printf("my_id=%d, my_dim=%d, my_dim_start=%d\n", my_id, my_dim, my_dim_start);
}

data = fftw_alloc_complex(alloc_local);

// create plans for in-place FFT
	w1=dptime_mpi();
plan = fftw_mpi_plan_dft_3d(nx_dim, ny_dim, nz_dim, data, data, MPI_COMM_WORLD,
		FFTW_FORWARD, FFTW_ESTIMATE);
	w2=dptime_mpi();
pinv = fftw_mpi_plan_dft_3d(nx_dim, ny_dim, nz_dim, data, data, MPI_COMM_WORLD,
		FFTW_BACKWARD, FFTW_ESTIMATE);
	w3=dptime_mpi();

// create plans for out of place FFT
//	plan = fftw_mpi_plan_dft_3d(nx_dim, ny_dim, nz_dim, in, out, MPI_COMM_WORLD,
//		FFTW_FORWARD, FFTW_ESTIMATE);
//	pinv = fftw_mpi_plan_dft_3d(nx_dim, ny_dim, nz_dim, out, in, MPI_COMM_WORLD,
//		FFTW_BACKWARD, FFTW_ESTIMATE);

// set the input data 
ip=0;
for (i=0; i<my_dim; ++i) {
for (j=0; j<ny_dim; ++j) {
for (k=0; k<nz_dim; ++k) {
	//	data[i*nz_dim*ny_dim + j*nz_dim +k ] = (10*j+k) + (i + 100*my_id)*I;
	data[ip] = ip + my_id*I;
	ip++;
}
}
}
// check print
if (debug_print == 2) { MPI_Barrier(MPI_COMM_WORLD);
for (i=0; i<my_dim; ++i) {
for (j=0; j<ny_dim; ++j) {
	printf("my_id=%d in  i=%d j=%d k: ",my_id,i,j);
for (k=0; k<nz_dim; ++k) {
	printf("(%5.1f,%5.1f) ",
	creal(data[i*nz_dim*ny_dim + j*nz_dim +k]),
	cimag(data[i*nz_dim*ny_dim + j*nz_dim +k]));
	}
	printf("\n");
	}
	}
}



// forward transformation
	w4=dptime_mpi();
fftw_execute(plan);
	w5=dptime_mpi();

// check print
if (debug_print == 2) { MPI_Barrier(MPI_COMM_WORLD);
for (i=0; i<my_dim; ++i) {
for (j=0; j<ny_dim; ++j) {
	printf("my_id=%d FWD i=%d j=%d k: ",my_id,i,j);
for (k=0; k<nz_dim; ++k) {
	printf("(%5.1f,%5.1f) ",
	creal(data[i*nz_dim*ny_dim + j*nz_dim +k]),
	cimag(data[i*nz_dim*ny_dim + j*nz_dim +k]));
	}
	printf("\n");
	}
	}
}

// backward synthesis
	w6=dptime_mpi();
fftw_execute(pinv);
	w7=dptime_mpi();

// normalization
for (i=0; i<my_dim; ++i) {
for (j=0; j<ny_dim; ++j) {
for (k=0; k<nz_dim; ++k) {
	data[i*nz_dim*ny_dim + j*nz_dim +k ] /= (double)(nx_dim*ny_dim*nz_dim);
}
}
}

// check print
if (debug_print == 2) { MPI_Barrier(MPI_COMM_WORLD);
for (i=0; i<my_dim; ++i) {
for (j=0; j<ny_dim; ++j) {
	printf("my_id=%d BACK i=%d j=%d k: ",my_id,i,j);
for (k=0; k<nz_dim; ++k) {
	printf("(%5.1f,%5.1f) ",
	creal(data[i*nz_dim*ny_dim + j*nz_dim +k]),
	cimag(data[i*nz_dim*ny_dim + j*nz_dim +k]));
	}
	printf("\n");
	}
	}
}

MPI_Barrier(MPI_COMM_WORLD);
if (my_id == 0) {
	printf("3D complex FFT(%d,%d,%d) nthreads=%d\n", nx_dim,ny_dim,nz_dim, nthreads);
	printf("planning time: FWD=%f, BACK=%f\n", w2-w1, w3-w2);
	printf("computing time: FWD=%f, BACK=%f\n", w5-w4, w7-w6);
}

fftw_destroy_plan(plan);
fftw_destroy_plan(pinv);
MPI_Finalize();
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
double dptime_mpi()
{
	return MPI_Wtime();
}
