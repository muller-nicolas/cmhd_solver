!***********************************************************
Program RadialSpectrum
implicit none
integer, parameter :: N = 256
integer, parameter :: Nh = N/2+1
integer, parameter :: Na = N/3  ! partie entière
double precision spec1d(Nh), spec2d(Nh,N), deltaT, kinj
integer i, j, k, countk(Nh), Nmax, istore, ii
integer ndeltaT, inrj, ispec
character (len=22) :: anim2D='out_spectrumrho-2D-'
character (len=22) :: anim1D='out_spectrumrho-1D-'

open(30, file = 'out_parameter', status = 'old',form='formatted')
read(30,*) deltaT, ndeltaT, inrj, kinj, ispec
close(30)

Nmax = int(ndeltaT/ispec)
istore = 100

do ii = 1, Nmax
write(anim2D(20:22),'(i3)') istore
open(30, file = anim2D, status = 'old',form='formatted')
do i = 1, Nh
do j = 1, N
read(30,*) spec2d(i,j)
end do
end do
close(30)
spec1d = 0.
countk = 0
do i = 1, Nh
do j = 1, Nh
k = ceiling(sqrt(real((i-1)*(i-1)+(j-1)*(j-1)))-0.5)
if (k .le. Nh) then
spec1d(k) = spec2d(i,j) + spec1d(k)
countk(k) = countk(k) + 1
end if
end do
end do
spec1d = spec1d/(real(countk)+1.d-40)
write(anim1D(20:22),'(i3)') istore
open(31, file = anim1D, status = 'new',form='formatted')
write(31,*) spec1d(:)
close(31)
istore = istore + 1
end do

end program RadialSpectrum
