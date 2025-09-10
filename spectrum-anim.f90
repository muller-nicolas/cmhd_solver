!***********************************************************
Program RadialSpectrum
use parameters
implicit none
double precision spec1d(Nh), spec2d(Nh,N)
integer i, j, k, countk(Nh), Nmax, ii, istore
character (len=21) :: anim2D='out_spectrumEU-2D-'
character (len=21) :: anim1D='out_spectrumEU-1D-'

! print *,istore_sp
! open(30, file = 'out_parameter', status = 'old',form='formatted')
! read(30,*) deltaT, ndeltaT, inrj, kinj, ispec, ifields
! close(30)

istore = 100
Nmax = int(ndeltaT/ispec)

! TODO: loops in spec are wrong, but it works...
do ii = 1, Nmax
write(anim2D(19:21),'(i3)') istore
open(30, file = anim2D, status = 'old',form='formatted')
do i = 1, Nh
do j = 1, N
read(30,*) spec2d(i,j)
end do
end do
close(30)
spec1d = 0.
countk = 0
do i = 1, N
do j = 1, Nh
k = ceiling(sqrt(real((i-1)*(i-1)+(j-1)*(j-1)))-0.5)
if (k .le. Nh) then
spec1d(k) = spec2d(j,i) + spec1d(k)
countk(k) = countk(k) + 1
end if
end do
end do
spec1d = spec1d/(real(countk)+1.d-40)
write(anim1D(19:21),'(i3)') istore
open(31, file = anim1D, status = 'new',form='formatted')
write(31,*) spec1d(:)
close(31)
istore = istore + 1
end do

end program RadialSpectrum
