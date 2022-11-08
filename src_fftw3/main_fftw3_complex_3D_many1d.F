#define ENABLE_THREADS true

program fftw3_3D_many1d
#ifdef ENABLE_THREADS
use omp_lib
#endif
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
integer :: iret, iaq1, iaq2
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

#ifdef ENABLE_THREADS
call dfftw_init_threads(iret)
nthreads=omp_get_max_threads()
call dfftw_plan_with_nthreads(nthreads)
#endif ENABLE_THREADS

write(*,'(a,a)') "starting <complex 3D FFT> using FFTW plan: ", env_string


do loop=1,10
nx_dim=2**loop
if (nx_dim.gt.max_dim_size) cycle
ny_dim=nx_dim
nz_dim=nx_dim
!cx if (nz_dim.ge.1000) nz_dim=nz_dim/10
if (nz_dim.ge.1000) nz_dim=512

allocate(in(nx_dim,ny_dim,nz_dim),stat=iaq1)
allocate(out(nx_dim,ny_dim,nz_dim),stat=iaq2)
if ((iaq1+iaq2).ne.0) then ; write(*,*) "allocate failed." ; stop ; endif

    w1=dptime_timef()
    call dfftw_plan_many_dft &
    &   (plan, 1, nx_dim, ny_dim*nz_dim, in, nx_dim, 1, nx_dim, out, nx_dim, 1, nx_dim, FFTW_FORWARD, env_integer);

    w2=dptime_timef()
    call dfftw_plan_many_dft &
    &   (pinv, 1, nx_dim, ny_dim*nz_dim, out, nx_dim, 1, nx_dim, in, nx_dim, 1, nx_dim, FFTW_BACKWARD, env_integer);
    w3=dptime_timef()


do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	!cx in(i,j,k) = dcmplx(i, -i)
    in(i,j,k) = dcmplx(i,-i) + 0.1*dcmplx(i-1) *0.1*dcmplx(0,(j-1)) * 0.1*dcmplx(k-1)
    if(i.gt.nx_dim/2) then
    in(i,j,k) = dcmplx(i,nx_dim-i) + 0.1*dcmplx(0,(j-1)) * 0.1*dcmplx(k-1)
    endif
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

call dfftw_destroy_plan ( plan)
call dfftw_destroy_plan ( pinv)
deallocate(in)
deallocate(out)

if (loop.le.2) cycle
write(*,'(a,i3,1x,a,3i5,2(2x,a,f15.9))') "loop",loop, ":",nx_dim,ny_dim,nz_dim, &
& "[plan(F+B)]=", real(w2-w1)+real(w3-w2), &
& "[execute(F+B)]=", real(w5-w4)+real(w7-w6)

end do  !cx outer most loop

stop
end



real(kind=8) function dptime_timef()
!cx integer(8) :: cv,cr
!cx call system_clock(count=cv,count_rate=cr)
!cx dptime_timef = real(cv,8)/real(cr,8)

!cx real(kind=8) timef
!cx	dptime_timef = timef()*0.001

!cx	use omp_lib
!cx	dptime_timef = omp_get_wtime()

call gettod(dptime_timef)
dptime_timef = dptime_timef * 1.0e-6
end

