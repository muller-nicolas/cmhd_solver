MODULE outputs
use, intrinsic :: iso_c_binding
use parameters
use fftw_mod
use adaptive_mod
use spectral_mod
implicit none

contains

SUBROUTINE save_energy(rhok,ukx,uky,bkx,bky)
double complex, intent(in) :: rhok(Nh,N), ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double precision :: EU, EB, Erho, Erho2, divu, divb

call compute_energy(rhok,ukx,uky,bkx,bky,EU,EB,Erho,Erho2,divu,divb)

open(40, file = 'out_EU', position='append',form='formatted')
open(41, file = 'out_EB', position='append',form='formatted')
open(42, file = 'out_Erho', position='append',form='formatted')
open(43, file = 'out_divu', position='append',form='formatted')
open(44, file = 'out_divb', position='append',form='formatted')
open(45, file = 'out_Erho2', position='append',form='formatted')
write(40,*) EU
write(41,*) EB
write(42,*) Erho
write(43,*) divu
write(44,*) divb
write(45,*) Erho2
close(40)
close(41)
close(42)
close(43)
close(44)
close(45)

END SUBROUTINE save_energy

!*****************************************************************
SUBROUTINE compute_energy(rhok,ukx,uky,bkx,bky,EU,EB,Erho,Erho2,divu,divb)
!***********compute energies
use parameters
use fftw_mod
implicit none
double complex, intent(in) :: rhok(Nh,N), ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double precision, intent(out) :: EU, EB, Erho, Erho2, divu, divb
double complex tmpk1(Nh,N), tmpk2(Nh,N) 
double precision tmp1(N,N), tmp2(N,N) 
double precision norm, norm2
integer i, j

norm = 1./real(N*N)
norm2 = norm*norm
EU = 0.
EB = 0.
Erho = 0.
Erho2 = 0.
divu = 0.
divb = 0.

! TODO: Erho is not properly computed
do i = 1, N
    EU = EU + 0.5*(abs(ukx(1,i))**2 + abs(uky(1,i))**2)
    EB = EB + 0.5*(abs(bkx(1,i))**2 + abs(bky(1,i))**2)
    Erho = Erho + 0.5*(1.d0+abs(rhok(1,i)))*(abs(ukx(1,i))**2 + abs(uky(1,i))**2) ! CHECK
    Erho2 = Erho2 + abs(rhok(1,i))
    do j = 2, Nh
        EU = EU + (abs(ukx(j,i))**2 + abs(uky(j,i))**2)
        EB = EB + (abs(bkx(j,i))**2 + abs(bky(j,i))**2)
        Erho = Erho + (1.d0+abs(rhok(j,i)))*(abs(ukx(j,i))**2 + abs(uky(j,i))**2)
        Erho2 = Erho2 + abs(rhok(j,i))
    end do
end do

EU = EU*norm2
EB = EB*norm2
Erho = Erho*norm2
Erho2 = Erho2*norm2

call divergence(ukx,uky,tmpk1)
call divergence(bkx,bky,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)

! Divergences in real space
do i = 1,N
    do j = 1,N
        divu = divu + tmp1(j,i)
        divb = divb + tmp2(j,i)
    end do
end do

divu = divu*norm
divb = divb*norm

RETURN
END SUBROUTINE compute_energy

!*****************************************************************
SUBROUTINE energyF(ux,uy,EU)
use parameters
!***********compute energy
implicit none
double precision EU, ux(N,N), uy(N,N)
integer i, j

EU = 0.
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
double precision spec1d(Nh)
integer, intent(inout) :: istore_sp
character (len=21) :: animU1D='out_spectrumEU-1D-'
character (len=21) :: animB1D='out_spectrumEB-1D-'
character (len=22) :: animrho1D='out_spectrumrho-1D-'

! Save kinetic energy spectrum
write(animU1D(19:21),'(i3)') istore_sp
call spectrum1D(ukx,uky,spec1d)
open(31, file=animU1D, status = 'new',form='formatted')
write(31,*) spec1d(:)
close(31)

! Save magnetic energy spectrum
write(animB1D(19:21),'(i3)') istore_sp
call spectrum1D(bkx,bky,spec1d)
open(31, file=animB1D, status='new',form='formatted')
write(31,*) spec1d(:)
close(31)

! Save density spectrum
call spectrumrho1D(rhok,spec1d)
write(animrho1D(20:22),'(i3)') istore_sp
open(31, file=animrho1D, status='new',form='formatted')
write(31,*) spec1d(:)
close(31)

istore_sp = istore_sp + 1

END SUBROUTINE save_spectra

SUBROUTINE save_fields(rhok,ukx,uky,bkx,bky,istore_fields)
double complex, intent(in) :: rhok(Nh,N), ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
integer, intent(inout) :: istore_fields
double complex :: tmpk(Nh,N)
double precision :: tmp(N,N)

character (len=14) :: animO='out_rho-2D-'
character (len=13) :: animW='out_wz-2D-'
character (len=13) :: animJ='out_jz-2D-'
character (len=15) :: animdiv='out_divb-2D-'
character (len=15) :: animdivu='out_divu-2D-'

write(animO(12:14),'(i3)') istore_fields
call FFT_SP(rhok,tmp)
open(30, file = animO, status='replace', form='unformatted', access='stream')
write(30) tmp(:,:)
close(30)
!
write(animW(11:13),'(i3)') istore_fields
call curl(ukx,uky,tmpk)
call FFT_SP(tmpk,tmp)
open(30, file = animW, status='replace', form='unformatted', access='stream')
write(30) tmp(:,:)
close(30)
!
call curl(bkx,bky,tmpk)
call FFT_SP(tmpk,tmp)
write(animJ(11:13),'(i3)') istore_fields
open(30, file = animJ, status='replace', form='unformatted', access='stream')
write(30) tmp(:,:)
close(30)
!
write(animdiv(13:15),'(i3)') istore_fields
call divergence(bkx,bky,tmpk)
call FFT_SP(tmpk,tmp)
open(30, file = animdiv, status='replace', form='unformatted', access='stream')
write(30) tmp(:,:)
close(30)
!
write(animdivu(13:15),'(i3)') istore_fields
call divergence(ukx,uky,tmpk)
call FFT_SP(tmpk,tmp)
open(30, file = animdivu, status='replace', form='unformatted', access='stream')
write(30) tmp(:,:)
close(30)
istore_fields = istore_fields + 1

END SUBROUTINE save_fields

END MODULE outputs