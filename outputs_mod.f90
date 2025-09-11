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
SUBROUTINE compute_energy(rho,ux,uy,bx,by,EU,EB,Erho,Erho2,divu,divb)
!***********compute energies
use parameters
use fftw_mod
implicit none
double complex, intent(in) :: rho(Nh,N), ux(Nh,N), uy(Nh,N), bx(Nh,N), by(Nh,N)
double precision, intent(out) :: EU, EB, Erho, Erho2, divu, divb
double complex ukxdx(Nh,N), ukydy(Nh,N), bkxdx(Nh,N), bkydy(Nh,N)
double precision uxdx(N,N), uydy(N,N), bxdx(N,N), bydy(N,N)
double precision norm, norm2
integer ii, jj

norm = 1./real(N*N)
norm2 = norm*norm
EU = 0.
EB = 0.
Erho = 0.
Erho2 = 0.
divu = 0.
divb = 0.

call derivex(ux, ukxdx)
call derivey(uy, ukydy)
call derivex(bx, bkxdx)
call derivey(bx, bkydy)

call FFT_SP(ukxdx,uxdx)
call FFT_SP(ukydy,uydy)
call FFT_SP(bkxdx,bxdx)
call FFT_SP(bkydy,bydy)

! TODO: Erho is not properly computed
do ii = 1, N
EU = EU + 0.5*(abs(ux(1,ii))**2 + abs(uy(1,ii))**2)
EB = EB + 0.5*(abs(bx(1,ii))**2 + abs(by(1,ii))**2)
Erho = Erho + 0.5*(1.d0+abs(rho(1,ii)))*(abs(ux(1,ii))**2 + abs(uy(1,ii))**2) ! CHECK
Erho2 = Erho2 + abs(rho(1,ii))
do jj = 2, Nh
EU = EU + (abs(ux(jj,ii))**2 + abs(uy(jj,ii))**2)
EB = EB + (abs(bx(jj,ii))**2 + abs(by(jj,ii))**2)
Erho = Erho + (1.d0+abs(rho(jj,ii)))*(abs(ux(jj,ii))**2 + abs(uy(jj,ii))**2)
Erho2 = Erho2 + abs(rho(jj,ii))
end do
end do

EU = EU*norm2
EB = EB*norm2
Erho = Erho*norm2
Erho2 = Erho2*norm2

! Divergences in real space
do ii = 1,N
do jj = 1,N
divu = divu + uxdx(jj,ii) + uydy(jj,ii)
divb = divb + bxdx(jj,ii) + bydy(jj,ii)
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
integer ii, jj

EU = 0.
do ii = 1, N
do jj = 1, N
EU = EU + (ux(jj,ii)**2 + uy(jj,ii)**2)
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
double complex :: ukydx(Nh,N), ukxdy(Nh,N), bkydx(Nh,N), bkxdy(Nh,N) 
double complex :: ukxdx(Nh,N), ukydy(Nh,N), bkxdx(Nh,N), bkydy(Nh,N) 
double precision :: uydx(N,N), uxdy(N,N), bydx(N,N), bxdy(N,N), rho(N,N)
double precision :: uxdx(N,N), uydy(N,N), bxdx(N,N), bydy(N,N) 

character (len=14) :: animO='out_rho-2D-'
character (len=13) :: animW='out_wz-2D-'
character (len=13) :: animJ='out_jz-2D-'
character (len=15) :: animdiv='out_divb-2D-'
character (len=15) :: animdivu='out_divu-2D-'

write(animO(12:14),'(i3)') istore_fields
call FFT_SP(rhok,rho)
open(30, file = animO, status='replace', form='unformatted', access='stream')
write(30) rho(:,:)
close(30)
!
write(animW(11:13),'(i3)') istore_fields
call derivex(uky,ukydx)
call derivey(ukx,ukxdy)
call FFT_SP(ukydx,uydx)
call FFT_SP(ukxdy,uxdy)
open(30, file = animW, status='replace', form='unformatted', access='stream')
write(30) uydx(:,:)-uxdy(:,:)
close(30)
!
call derivex(bky,bkydx)
call derivey(bkx,bkxdy)
call FFT_SP(bkydx,bydx)
call FFT_SP(bkxdy,bxdy)
write(animJ(11:13),'(i3)') istore_fields
open(30, file = animJ, status='replace', form='unformatted', access='stream')
write(30) bydx(:,:)-bxdy(:,:)
close(30)
!
write(animdiv(13:15),'(i3)') istore_fields
call derivex(bkx,bkxdx)
call derivey(bky,bkydy)
call FFT_SP(bkxdx,bxdx)
call FFT_SP(bkydy,bydy)
open(30, file = animdiv, status='replace', form='unformatted', access='stream')
write(30) bxdx(:,:)+bydy(:,:)
close(30)
!
write(animdivu(13:15),'(i3)') istore_fields
call derivex(ukx,ukxdx)
call derivey(uky,ukydy)
call FFT_SP(ukxdx,uxdx)
call FFT_SP(ukydy,uydy)
open(30, file = animdivu, status='replace', form='unformatted', access='stream')
write(30) uxdx(:,:)+uydy(:,:)
close(30)
istore_fields = istore_fields + 1

END SUBROUTINE save_fields

END MODULE outputs