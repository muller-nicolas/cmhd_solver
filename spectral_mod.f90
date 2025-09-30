MODULE spectral_mod
use, intrinsic :: iso_c_binding
use omp_lib
use parameters
use fftw_mod
implicit none

public
double precision, allocatable, dimension(:,:) :: kill, kd
double precision, allocatable, dimension(:) :: kx, ky
double precision :: kmax
save
contains

!*****************************************************************
SUBROUTINE Initk
! Initialization of the wavenumbers
double precision ks
integer i, j

allocate(kill(Nh,N), kd(Nh,N))
allocate(kx(N), ky(Nh))

ky=(/(dfloat(i-1)*dk,i=1,Nh,1)/)
kx(1:N/2)=(/(dfloat(i-1)*dk,i=1,N/2,1)/)
kx(N/2+1:N)=(/(dfloat(i-1-N)*dk,i=N/2+1,N,1)/)
do i = 1, N
kd(:,i) = kx(i)*kx(i) + ky(:)*ky(:)
end do

kmax = real(N)/3.*dk
kill(:,:) = 1.
!$omp parallel do private(i, j, ks) 
do i=1,N
    do j=1,Nh
        ks = sqrt(kd(j,i))
        if (ks.gt.kmax) then
            kill(j,i) = 0.
        end if
    end do
end do

RETURN
END SUBROUTINE Initk

!*****************************************************************
SUBROUTINE derivex(aa,bb)
! Computation of the x-derivative
double complex, intent(in)  :: aa(Nh,N)
double complex, intent(out) :: bb(Nh,N)
integer i,j

!$omp parallel do private(i, j) 
do i = 1, N
    do j = 1, Nh
        bb(j,i) = imag*aa(j,i)*kx(i)
    end do
end do

RETURN
END SUBROUTINE derivex

!*****************************************************************
SUBROUTINE derivey(aa,bb)
! Computation of the y-derivative
double complex, intent(in)  :: aa(Nh,N)
double complex, intent(out) :: bb(Nh,N)
integer i,j

!$omp parallel do private(i, j)
do i=1,N
    do j = 1, Nh
        bb(j,i) = imag*aa(j,i)*ky(j)
    end do
end do

RETURN
END SUBROUTINE derivey

!*****************************************************************
SUBROUTINE divergence(Akx,Aky,divk)
double complex, intent(in)  :: Akx(Nh,N), Aky(Nh,N)
double complex, intent(out) :: divk(Nh,N)
double complex  :: tmpk(Nh,N)

call derivex(Akx,tmpk)
divk = tmpk
call derivey(Aky,tmpk)
divk = divk + tmpk

RETURN
END SUBROUTINE divergence

!*****************************************************************
SUBROUTINE curl(Akx,Aky,curlk)
double complex, intent(in)  :: Akx(Nh,N), Aky(Nh,N)
double complex, intent(out) :: curlk(Nh,N)
double complex  :: tmpk(Nh,N)

call derivex(Aky,tmpk)
curlk = tmpk
call derivey(Akx,tmpk)
curlk = curlk - tmpk

RETURN
END SUBROUTINE curl

!*****************************************************************
SUBROUTINE cross_product(Ax,Ay,Bx,By,C)
double precision, dimension(N,N), intent(in) :: Ax, Ay, Bx, By
double precision, dimension(N,N), intent(out) :: C
C = Ax*By - Ay*Bx
RETURN
END SUBROUTINE cross_product

!*****************************************************************
SUBROUTINE spectrum1D(Akx,Aky,spec1d)
! Computes the 1D spectrum of the 2D fields Akx and Aky. Works for velocity and magnetic fields.
double precision :: spec1d(Nh) 
double complex, intent(in) :: Akx(Nh,N), Aky(Nh,N)
integer i, j, k 

spec1d = 0.
!$omp parallel do private(i,j,k) reduction(+:spec1d)
do i = 1, N
    k = int(sqrt(kd(1,i))/dk+0.5)+1
    if (k .le. Nh) then
        spec1d(k) = spec1d(k) + (abs(Akx(1,i))**2 + abs(Aky(1,i))**2)
    end if
    do j = 2, Nh-1
        k = int(sqrt(kd(j,i))/dk+0.5)+1
        if (k .le. Nh) then
            spec1d(k) = spec1d(k) + 2*(abs(Akx(j,i))**2 + abs(Aky(j,i))**2)
        end if
    end do
    k = int(sqrt(kd(Nh,i))/dk+0.5)+1
    if (k .le. Nh) then
        spec1d(k) = spec1d(k) + (abs(Akx(Nh,i))**2 + abs(Aky(Nh,i))**2)
    end if
end do
spec1d = spec1d/real(N,kind=8)**4

RETURN
END SUBROUTINE spectrum1D

!*****************************************************************
SUBROUTINE spectrumrho1D(rhok,spec1d)
! Computes the 1D spectrum of the 2D field rhok.
double precision :: spec1d(Nh) 
double complex, intent(in) :: rhok(Nh,N)
integer i, j, k 

spec1d = 0.
!$omp parallel do private(i,j,k) reduction(+:spec1d)
do i = 1, N
    k = int(sqrt(kd(1,i))/dk+0.5)+1
    if (k .le. Nh) then
        spec1d(k) = spec1d(k) + abs(rhok(1,i))**2
    end if
    do j = 2, Nh-1
        k = int(sqrt(kd(j,i))/dk+0.5)+1
        if (k .le. Nh) then
            spec1d(k) = spec1d(k) + 2*abs(rhok(j,i))**2
        end if
    end do
    k = int(sqrt(kd(Nh,i))/dk+0.5)
    if (k .le. Nh) then
        spec1d(k) = spec1d(k) + abs(rhok(Nh,i))**2
    end if
end do
spec1d = spec1d/real(N,kind=8)**4

RETURN
END SUBROUTINE spectrumrho1D

!*****************************************************************
SUBROUTINE spectrumAk(Akx,Aky,specE,specH)
! Computes the 1D spectrum of the 2D fields Akx and Aky. Works for velocity and magnetic fields.
double complex, intent(in) :: Akx(Nh,N), Aky(Nh,N)
double precision, intent(out) :: specE(Nh), specH(Nh)
integer i, j, k 

specE = 0.
specH = 0.
!$omp parallel do private(i,j,k) reduction(+:specE,specH)
do i = 1, N
    k = int(sqrt(kd(1,i))/dk+0.5)+1
    if (k .le. Nh) then
        specE(k) = specE(k) + 0.5*(abs(Akx(1,i))**2 + abs(Aky(1,i))**2)
        specH(k) = specH(k) + 0.5*(abs(Akx(1,i))**2 - abs(Aky(1,i))**2) * kx(i)/sqrt(kd(1,i))
    end if
    do j = 2, Nh-1
        k = int(sqrt(kd(j,i))/dk+0.5)+1
        if (k .le. Nh) then
            specE(k) = specE(k) + (abs(Akx(j,i))**2 + abs(Aky(j,i))**2)
            specH(k) = specH(k) + (abs(Akx(j,i))**2 - abs(Aky(j,i))**2) * kx(i)/sqrt(kd(j,i))
        end if
    end do
    k = int(sqrt(kd(Nh,i))/dk+0.5)+1
    if (k .le. Nh) then
        specE(k) = specE(k) + 0.5*(abs(Akx(Nh,i))**2 + abs(Aky(Nh,i))**2)
        specH(k) = specH(k) + 0.5*(abs(Akx(Nh,i))**2 - abs(Aky(Nh,i))**2) * kx(i)/sqrt(kd(j,i))
    end if
end do
specH(1) = 0.
specE = specE/real(N,kind=8)**4
specH = specH/real(N,kind=8)**4

RETURN
END SUBROUTINE spectrumAk

!*****************************************************************
SUBROUTINE canonical_variables(ukx,uky,bkx,bky,Ak1,Ak2)
double complex, intent(in) :: ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double complex, intent(out) :: Ak1(Nh,N), Ak2(Nh,N)
double complex :: tmpk1, tmpk2
double precision :: eps = 1.d-15
integer :: i,j

do i=1,N
    do j=1,Nh
        ! Ak1(j,i) = ky(j)/kx(i)*ukx(j,i) - uky(j,i) - sqrt(kd(j,i))/ky(j)*bkx(j,i)
        ! Ak2(j,i) = ky(j)/kx(i)*ukx(j,i) - uky(j,i) + sqrt(kd(j,i))/ky(j)*bkx(j,i)
        tmpk1 = ky(j)/kx(i)*ukx(j,i) - uky(j,i)
        tmpk2 = (kx(i)*bky(j,i) - ky(j)*bkx(j,i)) / (sqrt(kd(j,i))+eps)
        ! tmpk1 = ky(j)/kx(i)*ukx(j,i) - uky(j,i)
        ! tmpk2 = - bkx(j,i)*sqrt(kd(j,i)) / ky(j)
        Ak1(j,i) = tmpk1 + tmpk2
        Ak2(j,i) = tmpk1 - tmpk2
    end do
end do
Ak1(:,1) = 0.
Ak2(:,1) = 0.
! Ak1(1,:) = 0.
! Ak2(1,:) = 0.

RETURN
END SUBROUTINE canonical_variables

!*****************************************************************
SUBROUTINE Ek2D_and_Hk2D(Ak1,Ak2,Ek,Hk)
double complex, intent(in) :: Ak1(Nh,N), Ak2(Nh,N)
double complex, intent(out) :: Ek(Nh,N), Hk(Nh,N)
integer :: i,j

!$omp parallel do private(i,j)
do i=1,N
    do j=1,Nh
        Ek(j,i) = 0.25*(abs(Ak1(j,i))**2 + abs(Ak2(j,i))**2)
        Hk(j,i) =-0.25*(abs(Ak1(j,i))**2 - abs(Ak2(j,i))**2) * kx(i)/sqrt(kd(j,i))
    end do
end do

RETURN
END SUBROUTINE Ek2D_and_Hk2D

!*****************************************************************
SUBROUTINE WriteSpatioTemporalSpectrum(rhok, ukx, uky, bkx, bky, time)
double complex, intent(in) :: rhok(Nh,N), ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double complex :: Ak1(Nh,N), Ak2(Nh,N)
double precision :: time, tmpk1(Nh,N)
integer :: uSTS=80

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
character (len=12) :: sts_Ak1kx='STS_Ak1kx'
character (len=12) :: sts_Ak1ky='STS_Ak1ky'
character (len=12) :: sts_Ak2kx='STS_Ak2kx'
character (len=12) :: sts_Ak2ky='STS_Ak2ky'

OPEN(uSTS,file=sts_time,position='append',form='formatted')
write(uSTS,*) time
close(uSTS)

tmpk1(:,:) = 0.5*(abs(ukx)**2)
OPEN (uSTS, file=STS_Euxx_x, access='stream', position='append',form='unformatted')
write(uSTS) tmpk1(1,:) ! ky=0
close(uSTS)
OPEN (uSTS, file=STS_Euxx_y, access='stream', position='append',form='unformatted')
write(uSTS) tmpk1(:,1) ! kx=0
close(uSTS)

tmpk1(:,:) = 0.5*(abs(uky)**2)
OPEN (uSTS, file=STS_Euyy_x, access='stream', position='append',form='unformatted')
write(uSTS) tmpk1(1,:) ! ky=0
close(uSTS)
OPEN (uSTS, file=STS_Euyy_y, access='stream', position='append',form='unformatted')
write(uSTS) tmpk1(:,1) ! kx=0
close(uSTS)

tmpk1(:,:) = 0.5*(abs(bkx)**2)
OPEN (uSTS, file=STS_Ebxx_x, access='stream', position='append',form='unformatted')
write(uSTS) tmpk1(1,:) ! ky=0
close(uSTS)
OPEN (uSTS, file=STS_Ebxx_y, access='stream', position='append',form='unformatted')
write(uSTS) tmpk1(:,1) ! kx=0
close(uSTS)

tmpk1(:,:) = 0.5*(abs(bky)**2)
OPEN (uSTS, file=STS_Ebyy_x, access='stream', position='append',form='unformatted')
write(uSTS) tmpk1(1,:) ! ky=0
close(uSTS)
OPEN (uSTS, file=STS_Ebyy_y, access='stream', position='append',form='unformatted')
write(uSTS) tmpk1(:,1) ! kx=0
close(uSTS)

OPEN (uSTS, file=sts_uxkx, access='stream', position='append',form='unformatted')
write(uSTS) ukx(1,:) ! ky=0
close(uSTS)
OPEN (uSTS, file=sts_uxky, access='stream', position='append',form='unformatted')
write(uSTS) ukx(:,1) ! kx=0
close(uSTS)
OPEN (uSTS, file=sts_uykx, access='stream', position='append',form='unformatted')
write(uSTS) uky(1,:) !ky=0
close(uSTS)
OPEN (uSTS, file=sts_uyky, access='stream', position='append',form='unformatted')
write(uSTS) uky(:,1) !kx=0
close(uSTS)
OPEN (uSTS, file=sts_bxkx, access='stream', position='append',form='unformatted')
write(uSTS) bkx(1,:) ! ky=0
close(uSTS)
OPEN (uSTS, file=sts_bxky, access='stream', position='append',form='unformatted')
write(uSTS) bkx(:,1) ! kx=0
close(uSTS)
OPEN (uSTS, file=sts_bykx, access='stream', position='append',form='unformatted')
write(uSTS) bky(1,:) !ky=0
close(uSTS)
OPEN (uSTS, file=sts_byky, access='stream', position='append',form='unformatted')
write(uSTS) bky(:,1) !kx=0
close(uSTS)
OPEN (uSTS, file=sts_rhokx, access='stream', position='append',form='unformatted')
write(uSTS) rhok(1,:) !ky=0
close(uSTS)
OPEN (uSTS, file=sts_rhoky, access='stream', position='append',form='unformatted')
write(uSTS) rhok(:,1) !kx=0
close(uSTS)

! Canonical variables
call canonical_variables(ukx,uky,bkx,bky,Ak1,Ak2)
! Ak+
OPEN (uSTS, file=sts_Ak1kx, access='stream', position='append',form='unformatted')
write(uSTS) Ak1(2,:) !ky=0
close(uSTS)
OPEN (uSTS, file=sts_Ak1ky, access='stream', position='append',form='unformatted')
write(uSTS) Ak1(:,2) !kx=0
close(uSTS)
! Ak-
OPEN (uSTS, file=sts_Ak2kx, access='stream', position='append',form='unformatted')
write(uSTS) Ak2(2,:) !ky=0
close(uSTS)
OPEN (uSTS, file=sts_Ak2ky, access='stream', position='append',form='unformatted')
write(uSTS) Ak2(:,2) !kx=0
close(uSTS) 

END SUBROUTINE WriteSpatioTemporalSpectrum

END MODULE spectral_mod

