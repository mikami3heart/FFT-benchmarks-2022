program SSL2FFT_FULL_3D

!cx benchmark program to perform complex FFT using SSL2 FFT

integer, parameter :: ndim=3
integer isize(ndim),isign(ndim)
!cx integer, parameter :: max_dim_size=1024
integer, parameter :: max_dim_size=1024+1
integer :: nx_dim, ny_dim, nz_dim
!cx VASP benchmark size :: nx_dim=84, ny_dim=128, nz_dim=84

complex(kind=8), allocatable :: in(:,:,:)

integer i,j,k,l, n_many
real(kind=8) :: w1=0.0, w2=0.0, w3=0.0, w4=0.0, w5=0.0, w6=0.0, w7=0.0
real(kind=8) :: dptime_timef, dptime_mpi, dptime_omp
logical :: first_loop=.true.
integer :: nrepeat=1
integer :: nthreads=1
logical :: debug_print=.false.

write(*,'(a,i8,i5)') "starting SSL2 complex 3D FFT"

allocate(in(max_dim_size,max_dim_size,max_dim_size),stat=iaq1)
if (iaq1.ne.0) then ; stop "allocate failed." ; endif

do loop=2,10

nx_dim=2**loop
if (nx_dim.gt.max_dim_size) cycle
ny_dim=nx_dim
nz_dim=nx_dim
!cx if (nz_dim.ge.1000) nz_dim=nz_dim/10
if (nz_dim.ge.1000) nz_dim=512

!cx For DEBUG
!cx    nx_dim=8
!cx    ny_dim=4
!cx    nz_dim=2

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

      isize(1)=nx_dim
      isize(2)=ny_dim
      isize(3)=nz_dim
      do i=1,ndim
      isign(i)=1
      end do

w1=dptime_timef()
      call dm_v3dcft2 (in,max_dim_size,max_dim_size,nx_dim,ny_dim,nz_dim,isign,ierror)

w2=dptime_timef()

if (debug_print) then
	write(*,*) 'after the forward_transform:'
	do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	write(*,'(a,3(i2,1x),a,2f10.3)') ' in(',i,j,k,') = ',in(i,j,k)
	enddo; enddo; enddo
endif

w3=dptime_timef()
      do i=1,ndim
      isign(i)=-1
      end do
      call dm_v3dcft2 (in,max_dim_size,max_dim_size,nx_dim,ny_dim,nz_dim,isign,ierror)
w4=dptime_timef()

if (debug_print) then
	write(*,*) 'after the backward_synthesis:'
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

if (loop.le.2) cycle    !cx skip the cold first calls
write(*,'(a,i3,1x,a,3i5, 2x,a,3f15.9)') "loop", loop, ":", nx_dim,ny_dim,nz_dim, &
&   "time[F,B,(F+B)]=", real(w2-w1), real(w4-w3), real(w2-w1)+real(w4-w3)

end do

deallocate(in)

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

