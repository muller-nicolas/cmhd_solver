MODULE spectral_mod

use, intrinsic :: iso_c_binding
use fftw_mod
implicit none

double precision, allocatable, dimension(:,:) :: kill, kd
double precision, allocatable, dimension(:) :: kx, ky
save
contains

SUBROUTINE Initk(ddk,Nh,N)
! Initialization of the wavenumbers
double precision ddk, ks
integer Nh, N, ii, jj, kmax

allocate(kill(Nh,N))
allocate(kx(N))
allocate(ky(Nh))
allocate(kd(Nh,N))

ky=(/(dfloat(ii-1)*ddk,ii=1,Nh,1)/)
kx(1:N/2)=(/(dfloat(ii-1)*ddk,ii=1,N/2,1)/)
kx(N/2+1:N)=(/(dfloat(ii-1-N)*ddk,ii=N/2+1,N,1)/)
do ii = 1, N
kd(:,ii) = kx(ii)*kx(ii) + ky(:)*ky(:)
end do

kmax = int(N/3)
kill = 1.
do ii=1,N
    do jj=1,Nh
        ks = sqrt(kd(jj,ii))
        if (ks .gt. kmax) then
            kill(ii,jj) = 0.
        end if
    end do
end do

RETURN
END SUBROUTINE Initk

SUBROUTINE derivex(aa,bb,Nh,N)
! Computation of the x-derivative
complex*16, parameter :: imag = (0.0d0,1.0d0)
double precision aa(N,N), bb(N,N)
double complex :: ctmp(Nh,N)
integer Nh, N, ii

call dfftw_execute_dft_r2c(plan_for,aa,ctmp)
do ii = 1, N
ctmp(:,ii) = imag*ctmp(:,ii)*kx(ii)
end do
call dfftw_execute_dft_c2r(plan_back,ctmp,bb)
bb=bb/real(N*N)

RETURN
END SUBROUTINE derivex

!*****************************************************************
SUBROUTINE derivey(aa,bb,Nh,N)
! Computation of the y-derivative
complex*16, parameter :: imag = (0.0d0,1.0d0)
double precision aa(N,N), bb(N,N)
double complex :: ctmp(Nh,N)
integer Nh, N, jj

call dfftw_execute_dft_r2c(plan_for,aa,ctmp)
do jj = 1, Nh
ctmp(jj,:) = imag*ctmp(jj,:)*ky(jj)
end do
call dfftw_execute_dft_c2r(plan_back,ctmp,bb)
bb=bb/real(N*N)

RETURN
END SUBROUTINE derivey

!*****************************************************************
SUBROUTINE spectrum(ux,uy,Ek,Nh,N)
!***********compute the 2D spectrum
double precision ux(N,N), uy(N,N), Ek(Nh,N)
double complex :: ukx(Nh,N), uky(Nh,N), Ek1(Nh,N)
integer Nh, N, iia, iib

call dfftw_execute_dft_r2c(plan_for,ux,ukx)
call dfftw_execute_dft_r2c(plan_for,uy,uky)

Ek1 = abs(ukx)**2 + abs(uky)**2
do iia = 1, Nh
do iib = 1, N
Ek(iia,iib) = real(Ek1(iia,iib))
end do
end do

RETURN
END SUBROUTINE spectrum

!*****************************************************************
SUBROUTINE spectrumrho(rho,Ek,Nh,N)
!***********compute the 2D spectrum
double precision rho(N,N), Ek(Nh,N)
double complex :: rhok(Nh,N), Ek1(Nh,N)
integer Nh, N, iia, iib

call dfftw_execute_dft_r2c(plan_for,rho,rhok)
Ek1 = abs(rhok)**2
do iia = 1, Nh
do iib = 1, N
Ek(iia,iib) = real(Ek1(iia,iib))
end do
end do

RETURN
END SUBROUTINE spectrumrho

END MODULE spectral_mod
