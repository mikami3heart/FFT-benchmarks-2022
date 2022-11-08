
program thread_1D__fftw3
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

integer, parameter :: max2power=20
!cx integer, parameter :: max2power=5
!cx	integer, parameter :: max2power=2
integer*8 plan, pinv
complex(kind=8), allocatable :: in(:), out(:)
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

allocate(in(2**max2power),stat=iaq1)
allocate(out(2**max2power),stat=iaq2)
if ((iaq1+iaq2).ne.0) then ; write(*,*) "allocate failed. size=",2**max2power ; stop ; endif

#ifdef ENABLE_THREADS
call dfftw_init_threads(iret)
nthreads=omp_get_max_threads()
call dfftw_plan_with_nthreads(nthreads)
#endif ENABLE_THREADS

if (debug_print) then
    write(*,'(a,a,i3,$)') "starting <complex 1D FFT> plan= ",  env_string, env_integer
    write(*,'(3x,a,i2)') "using nthreads=", nthreads
endif


write(*,'(a)') "starting the 1D-benchmark. The first loop is a warm-up and should be discarded."


do loop=2,max2power
ndata=2**loop

w1=dptime_timef()
call dfftw_plan_dft_1d ( plan, ndata, in, out, FFTW_FORWARD, env_integer)
w2=dptime_timef()
call dfftw_plan_dft_1d ( pinv, ndata, out, in, FFTW_BACKWARD, env_integer)
w3=dptime_timef()

do i = 1,ndata
in(i) = dcmplx(real(i), -real(i))
enddo

if (debug_print) then
	write(*,*) 'Input array:'
	do i = 1,ndata
	write(*,'(a,i3,a,2f10.3)')  ' in(',i,') = ', in(i)
	enddo
endif


w4=dptime_timef()
call dfftw_execute ( plan)
w5=dptime_timef()

    if (debug_print) then
    	write(*,*) 'after forward FFT dft_1d:'
    	do i = 1,ndata
    	write(*,'(a,i3,a,2f10.3)')  ' out(',i,') = ',out(i)
    	enddo
    endif

w6=dptime_timef()
call dfftw_execute ( pinv)
w7=dptime_timef()

    if (debug_print) then
    	write(*,*) 'after backward FFT dft_1d:'
    	do i = 1,ndata
    	write(*,'(a,i3,a,2f10.3)')  ' in(',i,') = ', in(i)
    	enddo
    endif

do i = 1,ndata
in(i)=in(i)/real(ndata)
enddo

    if (debug_print) then
    	write(*,*) 'after the normalization:'
    	do i = 1,ndata
    	write(*,'(a,i3,a,2f10.3)')  ' in(',i,') = ', in(i)
    	enddo
    endif

call dfftw_destroy_plan ( plan)
call dfftw_destroy_plan ( pinv)

if (loop.le.2) cycle

write(*,'(a,i3,1x,a,i8, 2(2x,a,f15.9))') "loop", loop, ":", ndata, &
& "[plan(F+B)]=", real(w2-w1)+real(w3-w2), &
& "[execute(F+B)]=", real(w7-w6)+real(w5-w4)

end do

deallocate(in)
deallocate(out)
#ifdef ENABLE_MPI
call MPI_Finalize( ierr )
#endif

stop
end



real(kind=8) function dptime_timef()
integer(8) :: cv,cr
call gettod(dptime_timef)
dptime_timef = dptime_timef * 1.0e-6

!cx call system_clock(count=cv,count_rate=cr)
!cx dptime_timef = real(cv,8)/real(cr,8)

!cx real(kind=8) timef
!cx dptime_timef = timef()*0.001

!cx use omp_lib
!cx dptime_timef = omp_get_wtime()
end

