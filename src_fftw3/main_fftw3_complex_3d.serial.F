program fftw3_3D_serial
use omp_lib
implicit none
include 'fftw3.f'


integer, parameter :: max_dim_size=1024
!cx VASP benchmark size
!cx integer, parameter :: nx_dim=84, ny_dim=128, nz_dim=84
integer :: nx_dim, ny_dim, nz_dim
integer*8 plan, pinv

complex(kind=8), allocatable :: in(:,:,:), out(:,:,:)

integer i,j,k,loop
real(kind=8) :: w1=0.0, w2=0.0, w3=0.0, w4=0.0, w5=0.0, w6=0.0, w7=0.0
real(kind=8) :: dptime_timef, dptime_mpi, dptime_omp
logical :: first_loop=.true.
integer :: nthreads=1
integer :: iaq1, iaq2
logical :: debug_print=.false.

character*20 string, env_string
integer :: env_integer

!cx an environment variable MY_FFTW_PLAN should be given
!cx     FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE
!cx     quick and easy < ---------------------> long and optimized

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
if (debug_print) then
write(*,'(a,a,a)') string, "=", env_string
write(*,'(a,a)') "starting <complex 3D FFT> using FFTW plan: ", string
endif

write(*,'(a)') "starting the benchmark. The first loop is a warm-up and should be discarded."


do loop=4,10
nx_dim=2**loop

if (nx_dim.gt.max_dim_size) cycle
ny_dim=nx_dim
nz_dim=nx_dim
!cx if (nz_dim.ge.1000) nz_dim=nz_dim/10
if (nz_dim.ge.1000) nz_dim=512

write(*,'(a,i3,3x,a,3i5)') "loop", loop, "nx_dim, ny_dim, nz_dim=", nx_dim,ny_dim,nz_dim

allocate(in(nx_dim,ny_dim,nz_dim),stat=iaq1)
allocate(out(nx_dim,ny_dim,nz_dim),stat=iaq2)
if ((iaq1+iaq2).ne.0) then ; write(*,*) "allocate failed." ; stop ; endif

    w1=dptime_timef()
call dfftw_plan_dft_3d(plan,nx_dim,ny_dim,nz_dim, in, out, FFTW_FORWARD, env_integer)
    w2=dptime_timef()
call dfftw_plan_dft_3d(pinv,nx_dim,ny_dim,nz_dim, out, in, FFTW_BACKWARD, env_integer)
    w3=dptime_timef()

write(*,'(a,i5,2f12.6)') "[plan] size(nx), time(fwd), time(back) =", &
& nx_dim, real(w2-w1), real(w3-w2)


!cx	in  = (0.0, 0.0)
!cx	out = (0.0, 0.0)
do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	in(i,j,k) = dcmplx(i, 1 + (j-1)*nx_dim + (k-1)*nx_dim*ny_dim)
enddo; enddo; enddo

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
	write(*,*) 'after the forward_transform:'
	do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	write(*,'(a,3(i2,1x),a,2f10.3)') ' out(',i,j,k,') = ',out(i,j,k)
	enddo; enddo; enddo
endif

    w6=dptime_timef()
call dfftw_execute_dft ( pinv,out,in)
    w7=dptime_timef()

if (debug_print) then
	write(*,*) 'after the backward_synthesis:'
	do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	write(*,'(a,3(i2,1x),a,2f10.3)') ' in(',i,j,k,') = ',in(i,j,k)
	enddo; enddo; enddo
endif

do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	in(i,j,k) = in(i,j,k)/(real(nx_dim))
enddo; enddo; enddo

if (debug_print) then
	write(*,*) 'after the normalization:'
	do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	write(*,'(a,3(i2,1x),a,2f10.3)') ' in(',i,j,k,') = ',in(i,j,k)
	enddo; enddo; enddo
endif

write(*,'(a,i5,2f12.6)') "[execute] size(nx), time(fwd),time(back)=", &
& nx_dim, real(w5-w4), real(w7-w6)

call dfftw_destroy_plan ( plan)
call dfftw_destroy_plan ( pinv)
deallocate(in)
deallocate(out)
end do  !cx outer most loop

stop
end



real(kind=8) function dptime_timef()
integer(8) :: cv,cr
call system_clock(count=cv,count_rate=cr)
dptime_timef = real(cv,8)/real(cr,8)

!cx real(kind=8) timef
!cx	dptime_timef = timef()*0.001

!cx	use omp_lib
!cx	dptime_timef = omp_get_wtime()
end

