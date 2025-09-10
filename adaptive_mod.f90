MODULE adaptive_mod
use, intrinsic :: iso_c_binding
use parameters
use fftw_mod
implicit none

public 
double precision :: deltaT, nu
save

contains
!*****************************************************************
!     Adaptive timestep
!*****************************************************************
SUBROUTINE adaptiveT(ak,bk,ta)
double precision a(N,N), b(N,N), ta, max1, max2, max3, max4, max
double complex ak(Nh,N), bk(Nh,N)
integer ii, jj

call FFT_SP(ak,a)
call FFT_SP(bk,b)

max1 = dabs(a(1,1))
do ii = 1, N
do jj = 1, N
max2 = dabs(a(jj,ii))
if (max2 .gt. max1) then
max1=max2
end if
end do
end do

max3 = dabs(b(1,1))
do ii = 1, N
do jj = 1, N
max4 = dabs(b(jj,ii))
if (max4 .gt. max3) then
max3=max4
end if
end do
end do

max=max1
if (max3 .gt. max1) then
max=max3
end if
if (max .lt. 2.d0) then
max=1.d0
end if
ta=1.d0/(N*max)

RETURN
END SUBROUTINE adaptiveT

END MODULE adaptive_mod