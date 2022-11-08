program SSL2FFT

integer, parameter :: max_p2=20
integer, parameter :: max_length=2**max_p2
complex(kind=8), allocatable :: in(:)
complex(kind=8), allocatable :: out(:)
integer isign

integer i,j,k,l
real(kind=8) :: w1=0.0, w2=0.0, w3=0.0, w4=0.0, w5=0.0, w6=0.0, w7=0.0
real(kind=8) :: dptime_timef, dptime_mpi, dptime_omp
integer :: nthreads=1
logical :: debug_print=.false.

write(*,'(a,i8,i5)') "starting SSL2 complex 3D FFT thread version"

allocate(in(max_length),stat=iaq1)
allocate(out(max_length),stat=iaq2)
if ((iaq1+iaq2).ne.0) then ; write(*,*) "allocate failed." ; stop ; endif


!cx do loop=3,3      !cx for debugging the results
!cx do loop=5,max_p2
do loop=1,max_p2

nx_dim=2**loop

do i= 1,nx_dim
	!cx in(i) = dcmplx(i, -i)
	in(i) = dcmplx(i, (i+3)/4)
	out(i) = dcmplx(0,0)
enddo

if (debug_print) then
	write(*,*) 'input vector:'
	do i= 1,nx_dim
	write(*,'(a,i2,a,2f8.2)') ' in(',i,') = ', in(i)
	enddo
endif

w1=dptime_timef()
      isw=1
      isign=1
      call dm_v1dcft2(in,nx_dim,out,isign,ierror)
w2=dptime_timef()



if (debug_print) then
	write(*,*) 'after the forward_transform:'
	do i= 1,nx_dim
	write(*,'(a,i2,a,2f8.2)') ' in(',i,') = ', in(i)
	enddo
endif

!cx     write(*,*) 'starting the backward_synthesis:'
      isign=-1
w3=dptime_timef()
      call dm_v1dcft2(out,nx_dim,in,isign,ierror)
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
deallocate(out)

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

