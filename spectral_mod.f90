MODULE spectral_mod

use, intrinsic :: iso_c_binding
use fftw_mod
implicit none

contains

SUBROUTINE Initk(kkx,kky,kkkd,ddk,Nh,N)
! Initialization of the wavenumbers
double precision kkx(N), kky(Nh), kkkd(Nh,N), ddk
integer Nh, N, ii

kky=(/(dfloat(ii-1)*ddk,ii=1,Nh,1)/)
kkx(1:N/2)=(/(dfloat(ii-1)*ddk,ii=1,N/2,1)/)
kkx(N/2+1:N)=(/(dfloat(ii-1-N)*ddk,ii=N/2+1,N,1)/)
do ii = 1, N
kkkd(:,ii) = kkx(ii)*kkx(ii) + kky(:)*kky(:)
end do

RETURN
END SUBROUTINE Initk

SUBROUTINE derivex(aa,bb,kkx,Nh,N)
! Computation of the x-derivative
complex*16, parameter :: imag = (0.0d0,1.0d0)
double precision aa(N,N), bb(N,N), kkx(N)
double complex :: ctmp(Nh,N)
integer Nh, N, ii

call dfftw_execute_dft_r2c(plan_for,aa,ctmp)
do ii = 1, N
ctmp(:,ii) = imag*ctmp(:,ii)*kkx(ii)
end do
call dfftw_execute_dft_c2r(plan_back,ctmp,bb)
bb=bb/real(N*N)

RETURN
END SUBROUTINE derivex

!*****************************************************************
SUBROUTINE derivey(aa,bb,kky,Nh,N)
! Computation of the y-derivative
complex*16, parameter :: imag = (0.0d0,1.0d0)
double precision aa(N,N), bb(N,N), kky(Nh)
double complex :: ctmp(Nh,N)
integer Nh, N, jj

call dfftw_execute_dft_r2c(plan_for,aa,ctmp)
do jj = 1, Nh
ctmp(jj,:) = imag*ctmp(jj,:)*kky(jj)
end do
call dfftw_execute_dft_c2r(plan_back,ctmp,bb)
bb=bb/real(N*N)

RETURN
END SUBROUTINE derivey

END MODULE spectral_mod
