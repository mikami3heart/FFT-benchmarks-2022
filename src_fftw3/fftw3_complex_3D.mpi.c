// FFTW3 3D complex FFT sample program
//
// The original data is 1D array[nx_dim * ny_dim * nz_dim]
// 	equivalent to 3D array array[nx_dim][ny_dim][nz_dim]
//

#include <stdio.h>
#include <complex.h>
#include <fftw3-mpi.h>
#include <mpi.h>
double second_(void);
double dptime_timef(void);
double dptime_mpi(void);

#define MY_PLAN FFTW_MEASURE
// choose from FFTW_ESTIMATE FFTW_MEASURE FFTW_PATIENT FFTW_EXHAUSTIVE


#undef DEBUG_PRINT

int main(int argc, char **argv)
{
//	const ptrdiff_t nx_dim=2, ny_dim=2, nz_dim=2;
//	const ptrdiff_t nx_dim=2, ny_dim=2, nz_dim=4;
//	const ptrdiff_t nx_dim=64, ny_dim=64, nz_dim=64;
//	const ptrdiff_t nx_dim=128, ny_dim=128, nz_dim=128;
//	const ptrdiff_t nx_dim=256, ny_dim=256, nz_dim=256;
//	const ptrdiff_t nx_dim=512, ny_dim=512, nz_dim=512;
const ptrdiff_t nx_dim=1024, ny_dim=1024, nz_dim=512;

#ifdef DEBUG_PRINT
int debug_print=0;	// 0: suppress debug info, 2: print debug info
#endif

fftw_plan plan, pinv;
fftw_complex *cdata;
fftw_complex *cftdata;
ptrdiff_t alloc_local, my_dim, my_dim_start;

int i,j,k,ip;
int npes, my_id, nthreads=1;
double w1,w2,w3,w4,w5,w6,w7,w8;

MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
MPI_Comm_size(MPI_COMM_WORLD, &npes);

fftw_mpi_init();

#ifdef DEBUG_PRINT
if (debug_print == 2 && my_id == 0) {
	printf("\n\n 3D complex FFT(%d,%d,%d) flat MPI\n", nx_dim,ny_dim,nz_dim);
}
#endif


// set local size 
alloc_local = fftw_mpi_local_size_3d(nx_dim, ny_dim, nz_dim,  MPI_COMM_WORLD,
	&my_dim, &my_dim_start);

#ifdef DEBUG_PRINT
if (debug_print == 2) { MPI_Barrier(MPI_COMM_WORLD);
	printf("my_id=%d, my_dim=%d, my_dim_start=%d\n", my_id, my_dim, my_dim_start);
}
#endif

cdata   = fftw_alloc_complex(alloc_local);
cftdata = fftw_alloc_complex(alloc_local);


// create plans for out of place FFT
w1=dptime_mpi();
plan = fftw_mpi_plan_dft_3d(nx_dim, ny_dim, nz_dim, cdata, cftdata,
	MPI_COMM_WORLD, FFTW_FORWARD, MY_PLAN);
w2=dptime_mpi();
pinv = fftw_mpi_plan_dft_3d(nx_dim, ny_dim, nz_dim, cftdata, cdata,
	MPI_COMM_WORLD, FFTW_BACKWARD, MY_PLAN);
w3=dptime_mpi();


// set the input data 
ip=0;
for (i=0; i<my_dim; ++i) {
for (j=0; j<ny_dim; ++j) {
for (k=0; k<nz_dim; ++k) {
	cdata[ip] = my_id + ip + 1/(double)(ip+1) - (ip + 1/(double)(ip+1))*I;
	ip++;
}
}
}
#ifdef DEBUG_PRINT
if (debug_print == 2) { MPI_Barrier(MPI_COMM_WORLD);
for (i=0; i<my_dim; ++i) {
for (j=0; j<ny_dim; ++j) {
	printf("my_id=%d in  i=%d j=%d k: ",my_id,i,j);
for (k=0; k<nz_dim; ++k) {
	printf("(%5.1f,%5.1f) ",
	creal(cdata[i*nz_dim*ny_dim + j*nz_dim +k]),
	cimag(cdata[i*nz_dim*ny_dim + j*nz_dim +k]));
	}
	printf("\n");
	}
	}
}
#endif



// forward transformation
w4=dptime_mpi();
fftw_execute(plan);
w5=dptime_mpi();

#ifdef DEBUG_PRINT
// check print
if (debug_print == 2) { MPI_Barrier(MPI_COMM_WORLD);
for (i=0; i<my_dim; ++i) {
for (j=0; j<ny_dim; ++j) {
	printf("my_id=%d FWD i=%d j=%d k: ",my_id,i,j);
for (k=0; k<nz_dim; ++k) {
	printf("(%5.1f,%5.1f) ",
	creal(cftdata[i*nz_dim*ny_dim + j*nz_dim +k]),
	cimag(cftdata[i*nz_dim*ny_dim + j*nz_dim +k]));
	}
	printf("\n");
	}
	}
}
#endif

// backward synthesis
w6=dptime_mpi();
fftw_execute(pinv);
w7=dptime_mpi();

// normalization
for (i=0; i<my_dim; ++i) {
for (j=0; j<ny_dim; ++j) {
for (k=0; k<nz_dim; ++k) {
	cdata[i*nz_dim*ny_dim + j*nz_dim +k ] /= (double)(nx_dim*ny_dim*nz_dim);
}
}
}

#ifdef DEBUG_PRINT
// check print
if (debug_print == 2) { MPI_Barrier(MPI_COMM_WORLD);
for (i=0; i<my_dim; ++i) {
for (j=0; j<ny_dim; ++j) {
	printf("my_id=%d BACK i=%d j=%d k: ",my_id,i,j);
for (k=0; k<nz_dim; ++k) {
	printf("(%5.1f,%5.1f) ",
	creal(cdata[i*nz_dim*ny_dim + j*nz_dim +k]),
	cimag(cdata[i*nz_dim*ny_dim + j*nz_dim +k]));
	}
	printf("\n");
	}
	}
}
#endif

MPI_Barrier(MPI_COMM_WORLD);
if (my_id == 0) {
	printf("3D complex FFT(%d,%d,%d) using MY_PLAN \n", nx_dim,ny_dim,nz_dim );
	printf("time[plan] %f ", (w2-w1)+(w3-w2));
	printf("time[transformation] %f \n", (w5-w4)+(w7-w6));
}

fftw_destroy_plan(plan);
fftw_destroy_plan(pinv);
MPI_Finalize();
}

double dptime_mpi()
{
	return MPI_Wtime();
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
