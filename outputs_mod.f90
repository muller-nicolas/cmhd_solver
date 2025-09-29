MODULE outputs
use, intrinsic :: iso_c_binding
use omp_lib
use parameters
use fftw_mod
use adaptive_mod
use spectral_mod
implicit none

double complex, dimension(:,:), allocatable :: tmpk1, tmpk2, tmpk3
double precision, dimension(:,:), allocatable :: tmp1, tmp2, tmp3 
double precision, allocatable :: spec1d(:)

contains

SUBROUTINE save_energy(rhok,ukx,uky,bkx,bky)
double complex, intent(in) :: rhok(Nh,N), ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double precision :: EU, EB, Erho, Eint, divu, divb

call compute_energy(rhok,ukx,uky,bkx,bky,EU,EB,Erho,Eint,divu,divb)

open(40, file = 'out_EU', position='append',form='formatted')
open(41, file = 'out_EB', position='append',form='formatted')
open(42, file = 'out_Erho', position='append',form='formatted')
open(43, file = 'out_divu', position='append',form='formatted')
open(44, file = 'out_divb', position='append',form='formatted')
open(45, file = 'out_Eint', position='append',form='formatted')
write(40,*) EU
write(41,*) EB
write(42,*) Erho
write(43,*) divu
write(44,*) divb
write(45,*) Eint
close(40)
close(41)
close(42)
close(43)
close(44)
close(45)

END SUBROUTINE save_energy

!*****************************************************************
SUBROUTINE energy(Akx,Aky,Eout)
!***********compute energies
double complex, intent(in) :: Akx(Nh,N), Aky(Nh,N)
double precision, intent(out) :: Eout
double precision norm 
integer i, j

norm = 1.0d0/real(N,kind=8)**4
Eout = 0.0d0

!$omp parallel do private(i,j) reduction(+:Eout)
do i = 1, N
    Eout = Eout + 0.5*(abs(Akx(1,i))**2 + abs(Aky(1,i))**2)
    do j = 2, Nh-1
        Eout = Eout + (abs(Akx(j,i))**2 + abs(Aky(j,i))**2)
    end do
    Eout = Eout + 0.5*(abs(Akx(Nh,i))**2 + abs(Aky(Nh,i))**2)
end do
Eout = Eout*norm

RETURN
END SUBROUTINE energy

!*****************************************************************
SUBROUTINE compute_energy(rhok,ukx,uky,bkx,bky,EU,EB,Erho,Eint,divu,divb)
!***********compute energies
double complex, intent(in) :: rhok(Nh,N), ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double precision, intent(out) :: EU, EB, Eint, Erho, divu, divb
double precision norm, norm2
integer i, j

allocate (tmpk1(Nh,N), tmpk2(Nh,N), tmpk3(Nh,N))
allocate (tmp1(N,N), tmp2(N,N), tmp3(N,N))

norm = 1./real(N*N)
norm2 = norm*norm
EU = 0.
EB = 0.
Erho = 0.
divu = 0.
divb = 0.

! TODO: Erho is not properly computed
!$omp parallel do private(i,j) reduction(+:EU,EB)
do i = 1, N
    EU = EU + 0.5*(abs(ukx(1,i))**2 + abs(uky(1,i))**2)
    EB = EB + 0.5*(abs(bkx(1,i))**2 + abs(bky(1,i))**2)
    do j = 2, Nh-1
        EU = EU + (abs(ukx(j,i))**2 + abs(uky(j,i))**2)
        EB = EB + (abs(bkx(j,i))**2 + abs(bky(j,i))**2)
    end do
    EU = EU + 0.5*(abs(ukx(Nh,i))**2 + abs(uky(Nh,i))**2)
    EB = EB + 0.5*(abs(bkx(Nh,i))**2 + abs(bky(Nh,i))**2)
end do

EU = EU*norm2
EB = EB*norm2

tmpk1 = rhok
tmpk2 = ukx
tmpk3 = uky
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
call FFT_SP(tmpk3,tmp3)

!$omp parallel do private(i,j) reduction(+:Erho,Eint)
do i = 1, N
    do j = 1, N
        Erho = Erho + 0.5*(1.d0+tmp1(j,i))*(tmp2(j,i)**2 + tmp3(j,i)**2)
        Eint = Eint + cspeed**2/(gamma*(gamma-1.d0))*((tmp1(j,i)+1.d0)**gamma - 1.d0) ! Extract constant value
    end do
end do

Erho = Erho*norm
Eint = Eint*norm

call divergence(ukx,uky,tmpk1)
call divergence(bkx,bky,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)

! Divergences in real space
!$omp parallel do private(i,j) reduction(+:divu,divb)
do i = 1,N
    do j = 1,N
        divu = divu + tmp1(j,i)
        divb = divb + tmp2(j,i)
    end do
end do

divu = divu*norm
divb = divb*norm

deallocate (tmpk1, tmpk2, tmpk3, tmp1, tmp2, tmp3)

RETURN
END SUBROUTINE compute_energy

!*****************************************************************
SUBROUTINE energyF(ux,uy,EU)
!***********compute energy
double precision, intent(in) :: ux(N,N), uy(N,N)
double precision, intent(out) :: EU
integer i, j

EU = 0.
!$omp parallel do private(i,j) reduction(+:EU)
do i = 1, N
    do j = 1, N
        EU = EU + 0.5*(ux(j,i)**2 + uy(j,i)**2)
    end do
end do
EU = EU/real(N*N)

RETURN
END SUBROUTINE energyF

SUBROUTINE save_spectra(rhok,ukx,uky,bkx,bky,istore_sp)
double complex, intent(in) :: rhok(Nh,N), ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double complex :: Ak1(Nh,N), Ak2(Nh,N)
integer, intent(inout) :: istore_sp
character (len=21) :: animU1D='out_spectrumEU-1D-'
character (len=21) :: animB1D='out_spectrumEB-1D-'
character (len=22) :: animrho1D='out_spectrumrho-1D-'
character (len=22) :: animAk1='out_spectrumAk1-1D-'
character (len=22) :: animAk2='out_spectrumAk2-1D-'

allocate (spec1d(Nh))

! Save kinetic energy spectrum
write(animU1D(19:21),'(i3)') istore_sp
call spectrum1D(ukx,uky,spec1d)
open(31, file=animU1D, status = 'new',form='formatted')
write(31,'(1000(1X,E25.18))') spec1d(:)
close(31)

! Save magnetic energy spectrum
write(animB1D(19:21),'(i3)') istore_sp
call spectrum1D(bkx,bky,spec1d)
open(31, file=animB1D, status='new',form='formatted')
write(31,'(1000(1X,E25.18))') spec1d(:)
close(31)

! Save density spectrum
call spectrumrho1D(rhok,spec1d)
write(animrho1D(20:22),'(i3)') istore_sp
open(31, file=animrho1D, status='new',form='formatted')
write(31,'(1000(1X,E25.18))') spec1d(:)
close(31)

! Canonical variables
call canonical_variable(ukx,uky,bkx,bky,Ak1,Ak2)
! Ak+
call spectrumrho1D(Ak1,spec1d)
write(animAk1(20:22),'(i3)') istore_sp
open(31, file=animAk1, status='new',form='formatted')
write(31,'(1000(1X,E25.18))') spec1d(:)
close(31)
! Ak-
call spectrumrho1D(Ak2,spec1d)
write(animAk2(20:22),'(i3)') istore_sp
open(31, file=animAk2, status='new',form='formatted')
write(31,'(1000(1X,E25.18))') spec1d(:)
close(31)

istore_sp = istore_sp + 1

deallocate (spec1d)

END SUBROUTINE save_spectra

SUBROUTINE save_fields(rhok,ukx,uky,bkx,bky,istore_fields)
double complex, intent(in) :: rhok(Nh,N), ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
integer, intent(inout) :: istore_fields

character (len=14) :: animO='out_rho-2D-'
character (len=13) :: animUx='out_ux-2D-'
character (len=13) :: animUy='out_uy-2D-'
character (len=13) :: animBx='out_bx-2D-'
character (len=13) :: animBy='out_by-2D-'

allocate( tmpk1(Nh,N), tmp1(N,N))

tmpk1 = rhok
write(animO(12:14),'(i3)') istore_fields
call FFT_SP(tmpk1,tmp1)
open(30, file = animO, status='replace', form='unformatted', access='stream')
write(30) tmp1(:,:)
close(30)
!
tmpk1 = ukx
write(animUx(11:13),'(i3)') istore_fields
call FFT_SP(tmpk1,tmp1)
open(30, file = animUx, status='replace', form='unformatted', access='stream')
write(30) tmp1(:,:)
close(30)
!
tmpk1 = uky
write(animUy(11:13),'(i3)') istore_fields
call FFT_SP(tmpk1,tmp1)
open(30, file = animUy, status='replace', form='unformatted', access='stream')
write(30) tmp1(:,:)
close(30)
!
tmpk1 = bkx
write(animBx(11:13),'(i3)') istore_fields
call FFT_SP(tmpk1,tmp1)
open(30, file = animBx, status='replace', form='unformatted', access='stream')
write(30) tmp1(:,:)
close(30)
!
tmpk1 = bky
write(animBy(11:13),'(i3)') istore_fields
call FFT_SP(tmpk1,tmp1)
open(30, file = animBy, status='replace', form='unformatted', access='stream')
write(30) tmp1(:,:)
close(30)
!
istore_fields = istore_fields + 1

deallocate (tmpk1, tmp1)

END SUBROUTINE save_fields

END MODULE outputs