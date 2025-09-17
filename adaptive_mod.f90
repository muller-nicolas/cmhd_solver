MODULE adaptive_mod
use, intrinsic :: iso_c_binding
! use parameters
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
double complex ak(Nh,N), bk(Nh,N), tmpk1(Nh,N), tmpk2(Nh,N)
integer i, j

tmpk1 = ak
tmpk2 = bk
call FFT_SP(tmpk1,a)
call FFT_SP(tmpk2,b)

max1 = dabs(a(1,1))
do i = 1, N
    do j = 1, N
        max2 = dabs(a(j,i))
        if (max2 .gt. max1) then
            max1=max2
        end if
    end do
end do

max3 = dabs(b(1,1))
do i = 1, N
    do j = 1, N
        max4 = dabs(b(j,i))
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