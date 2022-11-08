program SSL2FFT

integer, parameter :: nx_dim=8, ny_dim=2, nz_dim=2

!cx complex(kind=8), allocatable :: in(:,:,:), out(:,:,:)
real(kind=8), allocatable :: in(:,:,:), out(:,:,:)

integer, allocatable :: isize(:),isign(:)

integer i,j,k,l, n_many
real(kind=8) :: w1=0.0, w2=0.0, w3=0.0, w4=0.0, w5=0.0, w6=0.0, w7=0.0
real(kind=8) :: dptime_timef, dptime_mpi, dptime_omp
logical :: debug_print=.true.
integer :: nrepeat=1
integer :: nthreads=1

write(*,'(a,i8,i5)') "starting SSL2 real to complex FFT"
write(*,'(a,4i5)') "nx_dim,ny_dim,nz_dim=", nx_dim,ny_dim,nz_dim

allocate(in(nx_dim+2,ny_dim,nz_dim),stat=iaq1)
allocate(out(nx_dim+2,ny_dim,nz_dim),stat=iaq2)
if ((iaq1+iaq2).ne.0) then ; write(*,*) "allocate failed." ; stop ; endif

!cx REAL    case      ld1 = 2*(n1/2+1)+npad1
!cx COMPLEX case      ld1 = n1+npad1

!cx	in  = (0.0, 0.0)
!cx	out = (0.0, 0.0)
do k = 1,nz_dim
do j = 1,ny_dim;
do i = 1,nx_dim
	in(i,j,k) = real(i) + real((j-1)*10) + real(k-1)*100
	out(i,j,k) = 1.0
enddo
enddo
enddo

if (debug_print) then
	write(*,*) 'input vector:'
	do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	write(*,'(a,3(i2,1x),a,2f10.3)') ' in(',i,j,k,') = ',in(i,j,k)
	enddo; enddo; enddo
endif

	ndim=3
	allocate(isize(ndim),stat=iaq1)
	allocate(isign(ndim),stat=iaq1)

      isize(1)=nx_dim
      isize(2)=ny_dim
      isize(3)=nz_dim
      isign(1)=1
      isign(2)=0
      isign(3)=0
      isn=isign(1)

    !cx a dummy call to discard the first overhead
      call dvmrf2(out,isize,ndim,isign,isn,ierror)

	write(*,*) 'starting the forward_transformation:'
w1=dptime_timef()
      call dvmrf2(in,isize,ndim,isign,isn,ierror)
w2=dptime_timef()

if (debug_print) then
	write(*,*) 'after the forward_transform:'
	do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	write(*,'(a,3(i2,1x),a,2f10.3)') ' in(',i,j,k,') = ',in(i,j,k)
	enddo; enddo; enddo
endif

      !cx   isign(1)=1
      isign(1)=-1
      isign(2)=0
      isign(3)=0
      isn=isign(1)

	write(*,*) 'starting the backward_synthesis:'
w3=dptime_timef()
      call dvmrf2(in,isize,ndim,isign,isn,ierror)
w4=dptime_timef()

if (debug_print) then
	write(*,*) 'after the backward_synthesis:'
	do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	write(*,'(a,3(i2,1x),a,2f10.3)') ' in(',i,j,k,') = ',in(i,j,k)
	enddo; enddo; enddo
endif

    scale = isize(1)* abs(isign(1)) + isize(2)* abs(isign(2)) + isize(3)* abs(isign(3))

do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	!cx in(i,j,k) = in(i,j,k)/(real(nx_dim*ny_dim)*nz_dim)
	in(i,j,k) = in(i,j,k)/scale
enddo; enddo; enddo

if (debug_print) then
	write(*,*) 'after the normalization:'
	do k = 1,nz_dim; do j = 1,ny_dim; do i = 1,nx_dim
	write(*,'(a,3(i2,1x),a,2f10.3)') ' in(',i,j,k,') = ',in(i,j,k)
	enddo; enddo; enddo
endif

write(*,'(a,2f12.6)') "ssl2 forward,backward=", real(w2-w1), real(w4-w3)

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

