
program thread_3D_fftw3
use omp_lib

#undef ENABLE_MPI
#define ENABLE_THREADS true

#ifdef ENABLE_THREADS
use omp_lib
#endif
#ifdef ENABLE_MPI
#include "mpif.h"
#endif

#include "fftw3.f"

integer, parameter :: max2power=10
!cx integer, parameter :: max2power=2
integer, parameter :: max_dim=2**max2power+1
integer*8 plan, pinv
complex(kind=8), allocatable :: in(:,:,:), out(:,:,:)
integer :: nx_dim, ny_dim, nz_dim

real(kind=8) :: w1, w2, w3, w4, w5, w6, w7
real(kind=8) :: dptime_timef, dptime_mpi, dptime_omp
logical :: debug_print=.false.

integer :: nthreads=1
character*20 string, env_string
integer :: env_integer

!cx     #define MY_FFTW_PLAN FFTW_ESTIMATE
!cx an environment variable MY_FFTW_PLAN should be given
!cx : FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE
!cx  quick and easy < ----------------------------> long and optimized

write(string,'(a)') "MY_FFTW_PLAN"
call getenv(string,env_string)
if (env_string.eq."FFTW_ESTIMATE") then
    env_integer = FFTW_ESTIMATE
else if (env_string.eq."FFTW_MEASURE") then
    env_integer = FFTW_MEASURE
else if (env_string.eq."FFTW_PATIENT") then
    env_integer = FFTW_PATIENT
else if (env_string.eq."FFTW_EXHAUSTIVE") then
    env_integer = FFTW_EXHAUSTIVE
else
    write(*,'(a)') "*** environment variable MY_FFTW_PLAN is not defined. FFTW_ESTIMATE is set. ***"
    env_integer = FFTW_ESTIMATE
endif
write(*,'(a,a,a)') string, "=", env_string

my_id=0
#ifdef ENABLE_MPI
call mpi_init(ierr )
call mpi_comm_rank( MPI_COMM_WORLD, my_id, ierr )
call mpi_comm_size( MPI_COMM_WORLD, npes, ierr )
#endif ENABLE_MPI

#ifdef ENABLE_THREADS
call dfftw_init_threads(iret)
nthreads=omp_get_max_threads()
call dfftw_plan_with_nthreads(nthreads)
#endif ENABLE_THREADS


allocate(in(max_dim, max_dim, max_dim/2+2), stat=iaq1)
allocate(out(max_dim, max_dim, max_dim/2+2), stat=iaq2)
if ((iaq1+iaq2).ne.0) then ; write(*,*) "allocate failed." ; stop ; endif
if (debug_print) then
    write(*,'(a,a,i3,$)') "starting <complex 1D FFT> plan= ",  env_string, env_integer
    write(*,'(3x,a,i2)') "using nthreads=", nthreads
endif
write(*,'(a)') "starting the 1D-benchmark. The first warm-up loop should be discarded."


do loop=2,max2power
nx_dim=2**loop
ny_dim=2**loop
nz_dim=2**loop
if (nz_dim.ge.1000) nz_dim=512    !cx to avoid OOM

w1=dptime_timef()
call dfftw_plan_dft_3d(plan,nx_dim,ny_dim,nz_dim, in, out, FFTW_FORWARD, env_integer)
w2=dptime_timef()
call dfftw_plan_dft_3d(pinv,nx_dim,ny_dim,nz_dim, out, in, FFTW_FORWARD, env_integer)
w3=dptime_timef()


do k = 1,nz_dim
do j = 1,ny_dim
do i = 1,nx_dim
    in(i,j,k) = dcmplx(i,-i) + 0.1*dcmplx(i-1) *0.1*dcmplx(0,(j-1)) * 0.1*dcmplx(k-1)
    if(i.gt.nx_dim/2) then
    in(i,j,k) = dcmplx(i,nx_dim-i) + 0.1*dcmplx(0,(j-1)) * 0.1*dcmplx(k-1)
    endif
enddo
enddo
enddo

if (debug_print) then
	write(*,*) 'input vector:'
	do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	write(*,'(a,3(i2,1x),a,2f10.3)') ' in(',i,j,k,') = ',in(i,j,k)
	enddo; enddo; enddo
endif


w4=dptime_timef()
call dfftw_execute_dft (plan,in,out)
w5=dptime_timef()

if (debug_print) then
	write(*,*) 'after the forward transform:'
	do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	write(*,'(a,3(i2,1x),a,2f10.3)') ' out(',i,j,k,') = ',out(i,j,k)
	enddo; enddo; enddo
endif


w6=dptime_timef()
call dfftw_execute_dft ( pinv,out,in)
w7=dptime_timef()

if (debug_print) then
	write(*,*) 'after the backward synthesis:'
	do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	write(*,'(a,3(i2,1x),a,2f10.3)') ' in(',i,j,k,') = ',in(i,j,k)
	enddo; enddo; enddo
endif


scale = real(nx_dim) * real(ny_dim) * real(nz_dim)
do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	in(i,j,k) = in(i,j,k)/scale
enddo; enddo; enddo

if (debug_print) then
	write(*,*) 'after the normalization:'
	do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	write(*,'(a,3(i2,1x),a,2f10.3)') ' in(',i,j,k,') = ',in(i,j,k)
	enddo; enddo; enddo
endif

call dfftw_destroy_plan ( plan)
call dfftw_destroy_plan ( pinv)

if (loop.le.2) cycle
write(*,'(a,i3,1x,a,3i5,$)') &
& "loop", loop, ":", nx_dim,ny_dim,nz_dim
write(*,'(a,f15.9,3x,a,f15.9)')  &
& "[plan(F+B)]=", real(w2-w1)+real(w3-w2), &
& "[execute(F+B)]=", real(w7-w6)+real(w5-w4)

end do

deallocate(in)
deallocate(out)
#ifdef ENABLE_MPI
call MPI_Finalize( ierr )
#endif ENABLE_MPI

stop
end



real(kind=8) function dptime_timef()
integer(8) :: cv,cr
!cx call system_clock(count=cv,count_rate=cr)
!cx dptime_timef = real(cv,8)/real(cr,8)

!cx real(kind=8) timef
!cx	dptime_timef = timef()*0.001

!cx	use omp_lib
!cx	dptime_timef = omp_get_wtime()

call gettod(dptime_timef)
dptime_timef = dptime_timef * 1.0e-6

end

