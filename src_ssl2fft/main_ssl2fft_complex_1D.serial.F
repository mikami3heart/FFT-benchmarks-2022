program SSL2FFT

integer, parameter :: max_p2=20
integer, parameter :: max_length=2**max_p2
complex(kind=8), allocatable :: in(:)
complex(kind=8), allocatable :: work(:)
integer isign

integer i,j,k,l
real(kind=8) :: w1=0.0, w2=0.0, w3=0.0, w4=0.0, w5=0.0, w6=0.0, w7=0.0
real(kind=8) :: dptime_timef, dptime_mpi, dptime_omp
integer :: nthreads=1
logical :: debug_print=.false.

write(*,'(a,i8,i5)') "starting SSL2 complex 1D FFT"

allocate(in(max_length),stat=iaq1)
allocate(work(max_length*2+70),stat=iaq2)
if ((iaq1+iaq2).ne.0) then ; write(*,*) "allocate failed." ; stop ; endif


!cx do loop=3,3      !cx for debugging the results
!cx do loop=5,max_p2
do loop=1,max_p2

nx_dim=2**loop

do i= 1,nx_dim
	!cx in(i) = dcmplx(i, (i+3)/4)
	in(i) = dcmplx(i, -i)
enddo

if (debug_print) then
	write(*,*) 'input vector:'
	do i= 1,nx_dim
	write(*,'(a,i2,a,2f8.2)') ' in(',i,') = ', in(i)
	enddo
endif

!cx     write(*,*) 'first call with isw=1. setting the trigonometric table'
!cx     write(*,*) 'isw is set as 0, and reused for succeeding calls'
!cx     write(*,*) 'starting the forward_transformation:'
w1=dptime_timef()
      isw=1
      isign=1
      call dvcfm1(in,nx_dim,isw,isign,work,ierror)
w2=dptime_timef()


!cx calling forward transformation twice will break the data. don't do that.
    !cx w3=dptime_timef()
      !cx   isw=0
      !cx   isign=1
      !cx   call dvcfm1(in,nx_dim,isw,isign,work,ierror)
    !cx w4=dptime_timef()

if (debug_print) then
	write(*,*) 'after the forward_transform:'
	do i= 1,nx_dim
	write(*,'(a,i2,a,2f8.2)') ' in(',i,') = ', in(i)
	enddo
endif

!cx     write(*,*) 'starting the backward_synthesis:'
      isign=-1
w3=dptime_timef()
      call dvcfm1(in,nx_dim,isw,isign,work,ierror)
w4=dptime_timef()

if (debug_print) then
	write(*,*) 'after the backward_synthesis: ierror=', ierror
	do i= 1,nx_dim
	write(*,'(a,i2,a,2f8.2)') ' in(',i,') = ', in(i)
	enddo
endif

do i= 1,nx_dim
	in(i) = in(i) / real(nx_dim)
enddo

if (debug_print) then
	write(*,*) 'after the normalization:'
	do i= 1,nx_dim
	write(*,'(a,i2,a,2f8.2)') ' in(',i,') = ', in(i)
	enddo
endif

if (loop.le.2) cycle    !cx skip the cold first calls
write(*,'(a,i3,1x,a,i8, 2x,a,3f15.9)') "loop", loop, ":", nx_dim, &
&   "time[F,B,(F+B)]=", real(w2-w1), real(w4-w3), real(w2-w1)+real(w4-w3)

end do

deallocate(in)
deallocate(work)

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

