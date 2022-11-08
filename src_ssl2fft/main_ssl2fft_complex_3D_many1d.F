program SSL2FFT

integer, parameter :: max_dim_size=1024
!cx integer, parameter :: max_dim_size=32
!cx VASP benchmark size
!cx integer, parameter :: nx_dim=84, ny_dim=128, nz_dim=84
integer :: nx_dim, ny_dim, nz_dim


complex(kind=8), allocatable :: in(:,:,:)
integer, allocatable :: isize(:),isign(:)

integer i,j,k,l, n_many
real(kind=8) :: w1=0.0, w2=0.0, w3=0.0, w4=0.0, w5=0.0, w6=0.0, w7=0.0
real(kind=8) :: dptime_timef, dptime_mpi, dptime_omp
logical :: first_loop=.true.
integer :: nrepeat=1
integer :: nthreads=1
logical :: debug_print=.false.

	ndim=3
	allocate(isize(ndim),stat=iaq1)
	allocate(isign(ndim),stat=iaq1)

write(*,'(a,i8,i5)') "starting SSL2 complex 3D FFT per x-direction"

do loop=5,10

nx_dim=2**loop
if (nx_dim.gt.max_dim_size) cycle
ny_dim=nx_dim
nz_dim=nx_dim
!cx if (nz_dim.ge.1000) nz_dim=nz_dim/10
if (nz_dim.ge.1000) nz_dim=512

write(*,'(a,i3,3x,a,3i5)') "loop", loop, "nx_dim, ny_dim, nz_dim=", nx_dim,ny_dim,nz_dim

allocate(in(nx_dim,ny_dim,nz_dim),stat=iaq1)
if (iaq1.ne.0) then ; write(*,*) "allocate failed." ; stop ; endif

do k = 1,nz_dim
do j = 1,ny_dim
do i = 1,nx_dim/2
	in(i,j,k) = dcmplx(i,1) + dcmplx(0,(j-2)*10) + dcmplx(k-1)*100
enddo
do i = nx_dim/2+1, nx_dim
	in(i,j,k) = dcmplx(i,2) + dcmplx(0,(j-2)*10) + dcmplx(k-1)*100
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
      isign(1)=0
      isign(2)=0
      isign(3)=0

    !cx a dummy call to practice and discard the first overhead
    if(first_loop) then
        isign(1)=1
        call dvmcf2(in,isize,ndim,isign,ierror)
        isign(1)=-1
        call dvmcf2(in,isize,ndim,isign,ierror)
        first_loop=.false.
    endif
    !cx end of the dummy call


w1=dptime_timef()
      isign(1)=1
      call dvmcf2(in,isize,ndim,isign,ierror)
w2=dptime_timef()

if (debug_print) then
	write(*,*) 'after the forward_transform:'
	do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	write(*,'(a,3(i2,1x),a,2f10.3)') ' in(',i,j,k,') = ',in(i,j,k)
	enddo; enddo; enddo
endif

w3=dptime_timef()
      isign(1)=-1
      call dvmcf2(in,isize,ndim,isign,ierror)
w4=dptime_timef()

if (debug_print) then
	write(*,*) 'after the backward_synthesis:'
	do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	write(*,'(a,3(i2,1x),a,2f10.3)') ' in(',i,j,k,') = ',in(i,j,k)
	enddo; enddo; enddo
endif

    !cx scale = isize(1)* abs(isign(1)) + isize(2)* abs(isign(2)) + isize(3)* abs(isign(3))
    scale = real(nx_dim)

do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	in(i,j,k) = in(i,j,k)/scale
enddo; enddo; enddo

if (debug_print) then
	write(*,*) 'after the normalization:'
	do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	write(*,'(a,3(i2,1x),a,2f10.3)') ' in(',i,j,k,') = ',in(i,j,k)
	enddo; enddo; enddo
endif

write(*,'(a,i5,i8,2f12.6)') "size(nx), many(ny*nz), time(fwd),time(back)=" &
& , nx_dim, ny_dim*nz_dim,  real(w2-w1), real(w4-w3)

deallocate(in)

end do

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

