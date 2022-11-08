program serial_3D_fftw
use omp_lib
include 'fftw3.f'

!cx	integer, parameter :: nx_dim=2, ny_dim=2, nz_dim=2
!cx	integer, parameter :: nx_dim=4, ny_dim=4, nz_dim=4
!cx	integer, parameter :: nx_dim=128, ny_dim=128, nz_dim=128
!cx	integer, parameter :: nx_dim=256, ny_dim=256, nz_dim=256
!cx	integer, parameter :: nx_dim=42, ny_dim=64, nz_dim=42
integer, parameter :: nx_dim=84, ny_dim=128, nz_dim=84
integer*8 plan, pinv
complex(kind=8), allocatable :: in(:,:,:), out(:,:,:)

integer i,j,k,l
real(kind=8) :: w1=0.0, w2=0.0, w3=0.0, w4=0.0, w5=0.0, w6=0.0, w7=0.0
real(kind=8) :: dptime_timef, dptime_mpi, dptime_omp
logical :: debug_print=.false.
integer :: nrepeat=100
integer :: nthreads=1

write(*,'(a,\)') "<complex 3D FFT>  "
write(*,'(a,4i5)') "size, threads=", nx_dim,ny_dim,nz_dim, nthreads
write(*,'(a,2i5)') "repeats fwd/back:", nrepeat

allocate(in(nx_dim,ny_dim,nz_dim),stat=iaq1)
allocate(out(nx_dim,ny_dim,nz_dim),stat=iaq2)
if ((iaq1+iaq2).ne.0) then ; write(*,*) "allocate failed." ; stop ; endif

!cx	!$omp parallel do 
!cx	!$omp end parallel do 
do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	in(i,j,k) = dcmplx(i + (j-1)*nx_dim + (k-1)*nx_dim*ny_dim, k)
enddo; enddo; enddo

if (debug_print) then
	write(*,*) 'input vector:'
	do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	write(*,'(a,3(i2,1x),a,2f10.3)') ' in(',i,j,k,') = ',in(i,j,k)
	enddo; enddo; enddo
endif

ninit=10
call f_pm_initialize (ninit)

!cx flag_for_plan : FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE
!cx             quick and easy < ----------------------------> long and optimized

call f_pm_start ("dfftw_plan")
w1=dptime_timef()
call dfftw_plan_dft_3d(plan,nx_dim,ny_dim,nz_dim, in, out, FFTW_FORWARD, FFTW_EXHAUSTIVE)
call f_pm_stop  ("dfftw_plan")
call f_pm_start ("dfftw_plan")
w2=dptime_timef()
call dfftw_plan_dft_3d(pinv,nx_dim,ny_dim,nz_dim, out, in, FFTW_BACKWARD, FFTW_EXHAUSTIVE)
call f_pm_stop  ("dfftw_plan")
w3=dptime_timef()
!cx write(*,'(a,2f12.6)') "<FFTW_ESTIMATE 3D plan> forward,backward=", real(w2-w1), real(w3-w2)
write(*,'(a,2f12.6)') "<FFTW_EXHAUSTIVE 3D plan> forward,backward=", real(w2-w1), real(w3-w2)


do l=1,nrepeat

call f_pm_start ("dfftw_execute(FWD)")
w4=dptime_timef()
call dfftw_execute_dft (plan,in,out)
w5=w5+(dptime_timef()-w4)
call f_pm_stop  ("dfftw_execute(FWD)")

if (debug_print) then
	write(*,*) 'after the forward_transform:'
	do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	!cx	write(*,*) ' out(',i,j,k,') = ',out(i,j,k)
	write(*,'(a,3(i2,1x),a,2f10.3)') ' out(',i,j,k,') = ',out(i,j,k)
	enddo; enddo; enddo
endif
!cx	in  = (0.0, 0.0)

call f_pm_start ("dfftw_execute(BACK)")
w6=dptime_timef()
call dfftw_execute_dft ( pinv,out,in)
w7=w7+(dptime_timef()-w6)
call f_pm_stop  ("dfftw_execute(BACK)")

if (debug_print) then
	do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	write(*,'(a,3(i2,1x),a,2f10.3)') ' in(',i,j,k,') = ',in(i,j,k)
	enddo; enddo; enddo
endif

	!cx	in  = in/(real(nx_dim*ny_dim)*nz_dim)
do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	in(i,j,k) = in(i,j,k)/(real(nx_dim*ny_dim)*nz_dim)
enddo; enddo; enddo

end do

write(*,'(a,f12.9)') "forward_transform  ave. time=",real(w5)/real(nrepeat)
write(*,'(a,f12.9)') "backward synthesis ave. time=",real(w7)/real(nrepeat)
if (debug_print) then
	write(*,*) 'after the normalization:'
	do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	write(*,'(a,3(i2,1x),a,2f10.3)') ' in(',i,j,k,') = ',in(i,j,k)
	enddo; enddo; enddo
endif

call f_pm_report ("")

call dfftw_destroy_plan ( plan)
call dfftw_destroy_plan ( pinv)
deallocate(in)
deallocate(out)

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

