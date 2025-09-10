MODULE spectral_mod

use, intrinsic :: iso_c_binding
use parameters
use fftw_mod
implicit none

public
double precision, allocatable, dimension(:,:) :: kill, kd
double precision, allocatable, dimension(:) :: kx, ky
save
contains

SUBROUTINE Initk
! Initialization of the wavenumbers
double precision ks, kmax
integer ii, jj

allocate(kill(Nh,N))
allocate(kx(N))
allocate(ky(Nh))
allocate(kd(Nh,N))

ky=(/(dfloat(ii-1)*dk,ii=1,Nh,1)/)
kx(1:N/2)=(/(dfloat(ii-1)*dk,ii=1,N/2,1)/)
kx(N/2+1:N)=(/(dfloat(ii-1-N)*dk,ii=N/2+1,N,1)/)
do ii = 1, N
kd(:,ii) = kx(ii)*kx(ii) + ky(:)*ky(:)
end do

kmax = real(N)/3*dk
kill(:,:) = 1.
do ii=1,N
    do jj=1,Nh
        ks = sqrt(kd(jj,ii))
        if (ks.gt.kmax) then
            kill(jj,ii) = 0.
        end if
    end do
end do

RETURN
END SUBROUTINE Initk

SUBROUTINE derivex(aa,bb)
! Computation of the x-derivative
! double precision, intent(in)  :: aa(N,N)
! double precision, intent(out) :: bb(N,N)
! double complex :: ctmp(Nh,N)
double complex, intent(in)  :: aa(Nh,N)
double complex, intent(out) :: bb(Nh,N)
integer ii

! call dfftw_execute_dft_r2c(plan_for,aa,ctmp)
do ii = 1, N
! ctmp(:,ii) = imag*ctmp(:,ii)*kx(ii)
bb(:,ii) = imag*aa(:,ii)*kx(ii)
end do
! call dfftw_execute_dft_c2r(plan_back,ctmp,bb)
! bb=bb/real(N*N)

RETURN
END SUBROUTINE derivex

!*****************************************************************
SUBROUTINE derivey(aa,bb)
! Computation of the y-derivative
double complex, intent(in)  :: aa(Nh,N)
double complex, intent(out) :: bb(Nh,N)
! double precision, intent(in)  :: aa(N,N)
! double precision, intent(out) :: bb(N,N)
! double complex :: ctmp(Nh,N)
integer jj

! call dfftw_execute_dft_r2c(plan_for,aa,ctmp)
do jj = 1, Nh
! ctmp(jj,:) = imag*ctmp(jj,:)*ky(jj)
bb(jj,:) = imag*aa(jj,:)*ky(jj)
end do
! call dfftw_execute_dft_c2r(plan_back,ctmp,bb)
! bb=bb/real(N*N)

RETURN
END SUBROUTINE derivey

!*****************************************************************
SUBROUTINE spectrum(ukx,uky,Ek)
!***********compute the 2D spectrum
! double precision, intent(in)  :: ukx(Nh,N), uky(Nh,N)
double precision, intent(out) :: Ek(Nh,N)
double complex :: ukx(Nh,N), uky(Nh,N), Ek1(Nh,N)
integer iia, iib

! call dfftw_execute_dft_r2c(plan_for,ux,ukx)
! call dfftw_execute_dft_r2c(plan_for,uy,uky)

Ek1 = abs(ukx)**2 + abs(uky)**2
do iia = 1, Nh
do iib = 1, N
Ek(iia,iib) = real(Ek1(iia,iib))
end do
end do

RETURN
END SUBROUTINE spectrum

!*****************************************************************
SUBROUTINE spectrumrho(rhok,Ek)
!***********compute the 2D spectrum
! double precision, intent(in)  :: rho(N,N)
double precision, intent(out) :: Ek(Nh,N)
double complex :: rhok(Nh,N), Ek1(Nh,N)
integer iia, iib

! call dfftw_execute_dft_r2c(plan_for,rho,rhok)
Ek1 = abs(rhok)**2
do iia = 1, Nh
do iib = 1, N
Ek(iia,iib) = real(Ek1(iia,iib))
end do
end do

RETURN
END SUBROUTINE spectrumrho

SUBROUTINE WriteSpatioTemporalSpectrum(rhok, ukx, uky, bkx, bky, time)
double complex :: rhok(Nh,N), ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double precision :: time, Euxx(Nh,N), Euyy(Nh,N), Ebxx(Nh,N), Ebyy(Nh,N)
integer :: uSTS
logical :: exist

character (len=11) :: sts_time='STS_time'
character (len=15) :: sts_Euxx_x='STS_Euxx_x'
character (len=15) :: sts_Euxx_y='STS_Euxx_y'
character (len=15) :: sts_Euyy_x='STS_Euyy_x'
character (len=15) :: sts_Euyy_y='STS_Euyy_y'
character (len=15) :: sts_Ebxx_x='STS_Ebxx_x'
character (len=15) :: sts_Ebxx_y='STS_Ebxx_y'
character (len=15) :: sts_Ebyy_x='STS_Ebyy_x'
character (len=15) :: sts_Ebyy_y='STS_Ebyy_y'
character (len=11) :: sts_uxkx='STS_uxkx'
character (len=11) :: sts_uxky='STS_uxky'
character (len=11) :: sts_uykx='STS_uykx'
character (len=11) :: sts_uyky='STS_uyky'
character (len=11) :: sts_bxkx='STS_bxkx'
character (len=11) :: sts_bxky='STS_bxky'
character (len=11) :: sts_bykx='STS_bykx'
character (len=11) :: sts_byky='STS_byky'
character (len=12) :: sts_rhokx='STS_rhokx'
character (len=12) :: sts_rhoky='STS_rhoky'

Euxx(:,:) = 0.5*(abs(ukx)**2)
Euyy(:,:) = 0.5*(abs(uky)**2)
Ebxx(:,:) = 0.5*(abs(bkx)**2)
Ebyy(:,:) = 0.5*(abs(bky)**2)
uSTS = 80
inquire(file=sts_time, exist=exist)

! if (exist) then
!     OPEN(30,file=sts_time,status="old",position='append',form='formatted')
! else
!     OPEN(30,file=sts_time,status="new",form='formatted')
! endif
OPEN(uSTS,file=sts_time,position='append',form='formatted')
write(uSTS,*) time
close(uSTS)

! if (exist) then
!     OPEN(uSTS,file=sts_uxkx,status="old",position='append',form='formatted')
! else
!     OPEN(uSTS,file=sts_uxkx,status="new",form='formatted')
! endif

! TODO: Expand complex conjugates when kx=0

OPEN (uSTS, file=STS_Euxx_x, access='stream', position='append',form='unformatted')
write(uSTS) Euxx(1,:) ! ky=0
close(uSTS)
OPEN (uSTS, file=STS_Euxx_y, access='stream', position='append',form='unformatted')
write(uSTS) Euxx(:,1) ! ky=0
close(uSTS)
OPEN (uSTS, file=STS_Euyy_x, access='stream', position='append',form='unformatted')
write(uSTS) Euyy(1,:) ! ky=0
close(uSTS)
OPEN (uSTS, file=STS_Euyy_y, access='stream', position='append',form='unformatted')
write(uSTS) Euyy(:,1) ! ky=0
close(uSTS)
OPEN (uSTS, file=STS_Ebxx_x, access='stream', position='append',form='unformatted')
write(uSTS) Ebxx(1,:) ! ky=0
close(uSTS)
OPEN (uSTS, file=STS_Ebxx_y, access='stream', position='append',form='unformatted')
write(uSTS) Ebxx(:,1) ! ky=0
close(uSTS)
OPEN (uSTS, file=STS_Ebyy_x, access='stream', position='append',form='unformatted')
write(uSTS) Ebyy(1,:) ! ky=0
close(uSTS)
OPEN (uSTS, file=STS_Ebyy_y, access='stream', position='append',form='unformatted')
write(uSTS) Ebyy(:,1) ! ky=0
close(uSTS)

! OPEN (uSTS, file=sts_uxkx, access='stream', position='append',form='unformatted')
! write(uSTS) ukx(1,:) ! ky=0
! close(uSTS)
! OPEN (uSTS, file=sts_uxky, access='stream', position='append',form='unformatted')
! write(uSTS) ukx(:,1) ! kx=0
! close(uSTS)
! OPEN (uSTS, file=sts_uykx, access='stream', position='append',form='unformatted')
! write(uSTS) uky(1,:) !ky=0
! close(uSTS)
! OPEN (uSTS, file=sts_uyky, access='stream', position='append',form='unformatted')
! write(uSTS) uky(:,1) !kx=0
! close(uSTS)
! OPEN (uSTS, file=sts_bxkx, access='stream', position='append',form='unformatted')
! write(uSTS) bkx(1,:) ! ky=0
! close(uSTS)
! OPEN (uSTS, file=sts_bxky, access='stream', position='append',form='unformatted')
! write(uSTS) bkx(:,1) ! kx=0
! close(uSTS)
! OPEN (uSTS, file=sts_bykx, access='stream', position='append',form='unformatted')
! write(uSTS) bky(1,:) !ky=0
! close(uSTS)
! OPEN (uSTS, file=sts_byky, access='stream', position='append',form='unformatted')
! write(uSTS) bky(:,1) !kx=0
! close(uSTS)
! OPEN (uSTS, file=sts_rhokx, access='stream', position='append',form='unformatted')
! write(uSTS) rhok(1,:) !ky=0
! close(uSTS)
! OPEN (uSTS, file=sts_rhoky, access='stream', position='append',form='unformatted')
! write(uSTS) rhok(:,1) !kx=0
! close(uSTS)

END SUBROUTINE WriteSpatioTemporalSpectrum

END MODULE spectral_mod

