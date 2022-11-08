program serial_1D__fftw3

#undef ENABLE_MPI
#undef ENABLE_THREADS

#ifdef ENABLE_THREADS
use omp_lib
#endif
#ifdef ENABLE_MPI
#include "mpif.h"
#endif

#include "fftw3.f"

integer, parameter :: max2power=16
!cx	integer, parameter :: max2power=2
integer*8 plan, pinv
complex(kind=8), allocatable :: in(:), out(:)
real(kind=8) :: w1, w2, w3, w4, w5, w6, w7
real(kind=8) :: dptime_timef, dptime_mpi, dptime_omp
logical :: debug_print=.false.
integer :: nrepeat=100
integer :: nthreads=1

#ifdef ENABLE_MPI
call mpi_init(ierr )
call mpi_comm_rank( MPI_COMM_WORLD, my_id, ierr )
call mpi_comm_size( MPI_COMM_WORLD, npes, ierr )
#endif ENABLE_MPI

ninit=10
call f_pm_initialize (ninit)

#ifdef ENABLE_THREADS
call f_pm_start ("dfftw_init_threads")
call dfftw_init_threads(iret)
nthreads=omp_get_max_threads()
call dfftw_plan_with_nthreads(nthreads)
call f_pm_stop  ("dfftw_init_threads")
write(*,'(a,i10)') "<complex 1D FFT using Threads> nthreads=", nthreads
#endif ENABLE_THREADS

allocate(in(2**max2power),stat=iaq1)
allocate(out(2**max2power),stat=iaq2)
if ((iaq1+iaq2).ne.0) then ; write(*,*) "allocate failed. size=",2**max2power ; stop ; endif


!cx	do loop=1,max2power
do loop=max2power,max2power
ndata=2**loop

write(*,'(a,i8,i5)') "<complex 1D FFT> size=", ndata
call f_pm_start ("dfftw_plan_1d")
w1=dptime_timef()
call dfftw_plan_dft_1d ( plan, ndata, in, out, FFTW_FORWARD, FFTW_EXHAUSTIVE)
w2=dptime_timef()
call dfftw_plan_dft_1d ( pinv, ndata, out, in, FFTW_BACKWARD, FFTW_EXHAUSTIVE)
w3=dptime_timef()
call f_pm_stop  ("dfftw_plan_1d")
!cx flag_for_plan : FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE
!cx             quick and easy < ----------------------------> long and optimized
!cx	write(*,'(a,2f12.6)') "<FFTW_ESTIMATE 1D plan> forward,backward=", real(w2-w1), real(w3-w2)
write(*,'(a,2f12.6)') "<FFTW_EXHAUSTIVE 1D plan> forward,backward=", real(w2-w1), real(w3-w2)

do i = 1,ndata
in(i) = dcmplx(real(i), -real(i))
enddo

if (debug_print) then
	write(*,*) 'Input array:'
	do i = 1,ndata
	write(*,*) ' in(',i,') = ',in(i)
	enddo
endif

w4=0.0; w5=0.0; w6=0.0; w7=0.0

do l=1,nrepeat

call f_pm_start ("dfftw_execute(FWD)")
w4=dptime_timef()
call dfftw_execute ( plan)
w5=w5+(dptime_timef()-w4)
call f_pm_stop  ("dfftw_execute(FWD)")

if (debug_print) then
	write(*,*) 'after forward FFT dft_1d:'
	do i = 1,ndata
	write(*,*) ' out(',i,') = ',out(i)
	enddo
endif

call f_pm_start ("dfftw_execute(BACK)")
w6=dptime_timef()
call dfftw_execute ( pinv)
w7=w7+(dptime_timef()-w6)
call f_pm_stop  ("dfftw_execute(BACK)")

do i = 1,ndata
in(i)=in(i)/real(ndata)
enddo

end do

write(*,'(a,f12.9)') "forward_transform  ave. time=",real(w5)/real(nrepeat)
write(*,'(a,f12.9)') "backward synthesis ave. time=",real(w7)/real(nrepeat)
if (debug_print) then
	write(*,*) 'after the normalization:'
	do i = 1,ndata
	write(*,*) ' in(',i,') = ',in(i)
	enddo
endif

call dfftw_destroy_plan ( plan)
call dfftw_destroy_plan ( pinv)
end do

call f_pm_report ("")
deallocate(in)
deallocate(out)
#ifdef ENABLE_MPI
call MPI_Finalize( ierr )
#endif

stop
end



real(kind=8) function dptime_timef()
integer(8) :: cv,cr
call system_clock(count=cv,count_rate=cr)
dptime_timef = real(cv,8)/real(cr,8)

!cx real(kind=8) timef
!cx dptime_timef = timef()*0.001

!cx use omp_lib
!cx dptime_timef = omp_get_wtime()
end

