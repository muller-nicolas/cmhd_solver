!*********************************************************************
!  Numerical integration of 2D CMHD equations for F-wave turbulence
!  Sébastien Galtier - LPP - Version 8 (July 2025)
!*********************************************************************

! TODO: Reduce number of fields (better scalability)
! TODO: Pointers for fields in memory (Reduce memory usage)
! TODO: Save fields as binary files (Reduce memory stockage)
! TODO: include spectrum-anim in the code (for practical purposes)
! TODO: handle precision globally (for practical purposes)
! TODO: Parallelisation: MPI, openMP or GPU
! TODO: RK4 (Numerical precision)

Program CMHD

use, intrinsic :: iso_c_binding
use parameters
use fftw_mod
use spectral_mod
use cMHD_mod

implicit none

double precision Ek(Nh,N)
double precision t1, t2
double precision divu, divb, ta, norm, time, timests

double precision :: rho0(N,N), ux0(N,N), uy0(N,N), bx0(N,N), by0(N,N)
double precision :: rho1(N,N), ux1(N,N), uy1(N,N), bx1(N,N), by1(N,N)
double precision :: uxdx(N,N), uxdy(N,N), uydx(N,N), uydy(N,N)
double precision :: bxdy(N,N), bydx(N,N), bxdx(N,N), bydy(N,N)
double precision :: nonlinrho0(N,N), nonlinux0(N,N), nonlinuy0(N,N), nonlinbx0(N,N), nonlinby0(N,N)
double precision :: nonlinrho1(N,N), nonlinux1(N,N), nonlinuy1(N,N), nonlinbx1(N,N), nonlinby1(N,N)

double precision EU, EB, Erho, Erho2 

double complex :: ukx0(Nh,N), uky0(Nh,N), bkx0(Nh,N), bky0(Nh,N), rhok0(Nh,N)
double complex :: ukx1(Nh,N), uky1(Nh,N), bkx1(Nh,N), bky1(Nh,N), rhok1(Nh,N)
double complex :: ukx2(Nh,N), uky2(Nh,N), bkx2(Nh,N), bky2(Nh,N), rhok2(Nh,N)
double complex :: nonlinrhok0(Nh,N), nonlinukx0(Nh,N), nonlinuky0(Nh,N), nonlinbkx0(Nh,N), nonlinbky0(Nh,N)
double complex :: nonlinrhok1(Nh,N), nonlinukx1(Nh,N), nonlinuky1(Nh,N), nonlinbkx1(Nh,N), nonlinbky1(Nh,N)

integer i, j, it 
character (len=11) :: animR='restart-'
character (len=21) :: animE='out_spectrumEU-2D-'
character (len=21) :: animB='out_spectrumEB-2D-'
character (len=22) :: animrho='out_spectrumrho-2D-'
character (len=14) :: animO='out_rho-2D-'
character (len=13) :: animW='out_wz-2D-'
character (len=13) :: animJ='out_jz-2D-'
! character (len=13) :: animux='out_ux-2D-'
! character (len=13) :: animuy='out_uy-2D-'
! character (len=13) :: animbx='out_bx-2D-'
! character (len=13) :: animby='out_by-2D-'
character (len=15) :: animdiv='out_divb-2D-'
character (len=15) :: animdivu='out_divu-2D-'

call cpu_time(time=t1)

!**************Initialization

time = 0.d0
timests = 0.d0

call init_fftw
call Initk

!***************** In case of no restart the code starts down here
if (nrestart .ne. 0) then

open(30, file='out_parameter', status='new', form='formatted')
write(30,*) deltaT, ndeltaT, inrj, kinj, ispec, ifields, N, dk
close(30)

open(40, file = 'out_EU', status = 'new',form='formatted')
open(41, file = 'out_EB', status = 'new',form='formatted')
open(42, file = 'out_Erho', status = 'new',form='formatted')
open(43, file = 'out_divu', status = 'new',form='formatted')
open(44, file = 'out_divb', status = 'new',form='formatted')
open(45, file = 'out_Erho2', status = 'new',form='formatted')
open(51, file = 'out_deltaT', status = 'new',form='formatted')
open(52, file = 'out_time', status = 'new',form='formatted')
open(53, file = 'out_nu', status = 'new',form='formatted')

! Initilize velocity field
call RandomInit(ukx0,uky0)

call RHS(rhok0,ukx0,uky0,bkx0,bky0,nonlinrhok0,nonlinukx0,nonlinuky0,nonlinbkx0,nonlinbky0)

rhok1 = rhok0 + deltaTi*nonlinrhok0
ukx1 = ukx0 + deltaTi*nonlinukx0
uky1 = uky0 + deltaTi*nonlinuky0
bkx1 = bkx0 + deltaTi*nonlinbkx0
bky1 = bky0 + deltaTi*nonlinbky0
end if

!****************** In case of restart the code starts below *************
if (nrestart .eq. 0) then
open(66, file = 'restart', status = 'old',form='unformatted')
read(66) rho1(:,:),ux1(:,:),uy1(:,:),bx1(:,:),by1(:,:),rho0(:,:),ux0(:,:),uy0(:,:),bx0(:,:),by0(:,:)
read(66) nonlinrho0(:,:),nonlinux0(:,:),nonlinuy0(:,:),nonlinbx0(:,:),nonlinby0(:,:)
read(66) nonlinrho1(:,:),nonlinux1(:,:),nonlinuy1(:,:),nonlinbx1(:,:),nonlinby1(:,:)
close(66)
call FFT_PS(rho1,rhok1)
call FFT_PS(ux1,ukx1)
call FFT_PS(uy1,uky1)
call FFT_PS(bx1,bkx1)
call FFT_PS(by1,bky1)
call FFT_PS(nonlinrho0,nonlinrhok0)
call FFT_PS(nonlinuy0,nonlinuky0)
call FFT_PS(nonlinuy0,nonlinuky0)
call FFT_PS(nonlinbx0,nonlinbkx0)
call FFT_PS(nonlinby0,nonlinbky0)
end if

!************************************************************************
!*****Time evolution: Adams-Bashforth
!*****************************  Main loop  ******************************
do it = 1, ndeltaT

! Random forcing in ux0 and uy0
! call RandomF(ux0,kinj,pi,it,Na,Nh,N)
! call RandomF(uy0,kinj,pi,it,Na,Nh,N)

! ! Normalization of the forcing
! call energyF(ux0,uy0,EU)
! ux0=amp*ux0/sqrt(EU)
! uy0=amp*uy0/sqrt(EU)

call RHS(rhok1,ukx1,uky1,bkx1,bky1,nonlinrhok1,nonlinukx1,nonlinuky1,nonlinbkx1,nonlinbky1)
call check_nan(rhok1)

! Adams-Bashford method
rhok2 = rhok1 + deltaT*(1.5*nonlinrhok1 - 0.5*nonlinrhok0)
ukx2  = ukx1  + deltaT*(1.5*nonlinukx1  - 0.5*nonlinukx0) + ukx0*off
uky2  = uky1  + deltaT*(1.5*nonlinuky1  - 0.5*nonlinuky0) + uky0*off
bkx2  = bkx1  + deltaT*(1.5*nonlinbkx1  - 0.5*nonlinbkx0)
bky2  = bky1  + deltaT*(1.5*nonlinbky1  - 0.5*nonlinbky0)

call dissipation(rhok2,ukx2,uky2,bkx2,bky2)

! Rename variables for saving and use them as initial values for next loop
rhok1=rhok2
ukx1=ukx2
uky1=uky2
bkx1=bkx2
bky1=bky2
nonlinrhok0=nonlinrhok1
nonlinukx0=nonlinukx1
nonlinuky0=nonlinuky1
nonlinbkx0=nonlinbkx1
nonlinbky0=nonlinbky1

! Compute and write energy
if ( (mod(it,inrj) .eq. 0) ) then
! call energy(rho1,ux1,uy1,bx1,by1,uxdx,uydy,bxdx,bydy,EU,EB,Erho,Erho2,divu,divb,N)
call derivex(ukx1,ukxdx)
call derivey(uky1,ukydy)
call derivex(bkx1,bkxdx)
call derivey(bky1,bkydy)
call energy(rhok1,ukx1,uky1,bkx1,bky1,ukxdx,ukydy,bkxdx,bkydy,EU,EB,Erho,Erho2,divu,divb)
write(40,*) EU
write(41,*) EB
write(42,*) Erho
write(43,*) divu
write(44,*) divb
write(45,*) Erho2

!****************Adaptive timestep
write(53,*) nu
time = time + dfloat(inrj)*deltaT
write(52,*) time
! call adaptiveT(ukx2,rhok2,ta)
! write(51,*) ta
! deltaT = ta*0.2d0 ! condition CFL
! nu = 1./ta*(1.d0/N)**4 ! condition CFL
! eta = nu

end if

! Compute and write spectra
if ( (mod(it,ispec) .eq. 0) ) then
print *, it
! TODO: Write spectra properly -> fix anim_spec routines
write(animE(19:21),'(i3)') istore_sp
call spectrum(ukx2,uky2,Ek)
open(30, file = animE, status = 'new',form='formatted')
do i = 1, Nh
    do j = 1, N
        write(30,*) Ek(i,j)
    end do
end do
close(30)
write(animB(19:21),'(i3)') istore_sp
call spectrum(bkx2,bky2,Ek)
open(30, file = animB, status = 'new',form='formatted')
do i = 1, Nh
    do j = 1, N
        write(30,*) Ek(i,j)
    end do
end do
close(30)
write(animrho(20:22),'(i3)') istore_sp
call spectrumrho(rhok2,Ek)
open(30, file = animrho, status = 'new',form='formatted')
do i = 1, Nh
    do j = 1, N
        write(30,*) Ek(i,j)
    end do
end do
close(30)
istore_sp = istore_sp + 1
endif

! Write fields
! rho, wz, jz, divu, divb
if (mod(it,ifields) .eq. 0) then
write(animO(12:14),'(i3)') istore_fields
call FFT_SP(rhok2,rho1)
open(30, file = animO, status = 'new',form='formatted')
write(30,*) rho1(:,:)
close(30)
!
write(animW(11:13),'(i3)') istore_fields
call derivex(uky2,ukydx)
call derivey(ukx2,ukxdy)
call FFT_SP(ukydx,uydx)
call FFT_SP(ukxdy,uxdy)
open(30, file = animW, status = 'new',form='formatted')
write(30,*) uydx(:,:)-uxdy(:,:)
close(30)
!
call derivex(bky2,bkydx)
call derivey(bkx2,bkxdy)
call FFT_SP(bkydx,bydx)
call FFT_SP(bkxdy,bxdy)
write(animJ(11:13),'(i3)') istore_fields
open(30, file = animJ, status = 'new',form='formatted')
write(30,*) bydx(:,:)-bxdy(:,:)
close(30)
!
write(animdiv(13:15),'(i3)') istore_fields
call derivex(bkx2,bkxdx)
call derivey(bky2,bkydy)
call FFT_SP(bkxdx,bxdx)
call FFT_SP(bkydy,bydy)
open(30, file = animdiv, status = 'new',form='formatted')
write(30,*) bxdx(:,:)+bydy(:,:)
close(30)
!
write(animdivu(13:15),'(i3)') istore_fields
call derivex(ukx2,ukxdx)
call derivey(uky2,ukydy)
call FFT_SP(ukxdx,uxdx)
call FFT_SP(ukydy,uydy)
open(30, file = animdivu, status = 'new',form='formatted')
write(30,*) uxdx(:,:)+uydy(:,:)
close(30)
istore_fields = istore_fields + 1
end if

! Save spatio-temporal spectra
if (sts .eq. 1) then
    if (mod(it,ists) .eq. 0) then
        timests = timests + dfloat(ists)*deltaT
        call WriteSpatioTemporalSpectrum(ukx1, uky1, bkx1, bky1, rhok1, timests)
    end if
end if


end do ! end of the temporal loop


! Save fields in real space for restart (maybe should I do it in spectral space?)
call FFT_SP(rhok1,rho1)
call FFT_SP(ukx1,ux1)
call FFT_SP(uky1,uy1)
call FFT_SP(bkx1,bx1)
call FFT_SP(bky1,by1)
call FFT_SP(rhok0,rho0)
call FFT_SP(ukx0,ux0)
call FFT_SP(uky0,uy0)
call FFT_SP(bkx0,bx0)
call FFT_SP(bky0,by0)
call FFT_SP(nonlinrhok0,nonlinrho0)
call FFT_SP(nonlinukx0,nonlinux0)
call FFT_SP(nonlinuky0,nonlinuy0)
call FFT_SP(nonlinbkx0,nonlinbx0)
call FFT_SP(nonlinbky0,nonlinby0)
call FFT_SP(nonlinrhok1,nonlinrho1)
call FFT_SP(nonlinukx1,nonlinux1)
call FFT_SP(nonlinuky1,nonlinuy1)
call FFT_SP(nonlinbkx1,nonlinbx1)
call FFT_SP(nonlinbky1,nonlinby1)

write(animR(9:11),'(i3)') istore_fields
open(30, file = animR, status = 'new',form='unformatted')
write(30) rho1(:,:),ux1(:,:),uy1(:,:),bx1(:,:),by1(:,:),rho0(:,:),ux0(:,:),uy0(:,:),bx0(:,:),by0(:,:)
write(30) nonlinrho0(:,:),nonlinux0(:,:),nonlinuy0(:,:),nonlinbx0(:,:),nonlinby0(:,:)
write(30) nonlinrho1(:,:),nonlinux1(:,:),nonlinuy1(:,:),nonlinbx1(:,:),nonlinby1(:,:)
close(30)

call end_fftw ! Deallocate plans

close(40)
close(41)
call cpu_time(time=t2)
write(*,*) "cpu time", t2-t1

print *, 'OK'
end program CMHD

!*****************************************************************
SUBROUTINE RandomInit(ukxi,ukyi)
use parameters
use fftw_mod
! Initial random field spectra
implicit none
double precision spectri(Na+1,Na+1), uxtmp(N,N), uytmp(N,N)
double precision theta, knc, EU
double complex :: ukxi(Nh,N), ukyi(Nh,N)
integer ii, jj

call srand(seed)
spectri=0.
ukxi=0. 
ukyi=0. 
do ii = 1, Na+1
do jj = 1, Na+1
knc = (kinj - real(ii-1) - real(jj-1))**4
spectri(jj,ii) = dexp(-knc*100.)
end do
end do
do ii = 1, Na+1
do jj = 1, Na+1
theta = rand()*2.*pi
ukxi(jj,ii) = spectri(jj,ii)*(cos(theta) + imag*sin(theta))
theta = rand()*2.*pi
ukyi(jj,ii) = spectri(jj,ii)*(cos(theta) + imag*sin(theta))
end do
end do
call FFT_SP(ukxi,uxtmp)
call FFT_SP(ukyi,uytmp)

! Normalization
call energyF(uxtmp,uytmp,EU)
uxtmp=amp*uxtmp/sqrt(EU)
uytmp=amp*uytmp/sqrt(EU)

call FFT_PS(uxtmp,ukxi)
call FFT_PS(uytmp,ukyi)

RETURN
END SUBROUTINE RandomInit

!*****************************************************************
! Random forcing spectrum
!*****************************************************************
SUBROUTINE RandomF(field1)
use parameters
use fftw_mod
implicit none
double precision spectri(Na+1,Na+1), field1(N,N)
double precision theta, knc
double complex :: spectric1(Nh,N)
integer ii, jj

call srand(seed)
spectri=0.
spectric1=0.
do ii = 1, Na+1
do jj = 1, Na+1
knc = (kinj - real(ii-1) - real(jj-1))**4
spectri(jj,ii) = dexp(-knc*100.)
end do
end do
do ii = 1, Na+1
do jj = 1, Na+1
theta = rand()*2.*pi
spectric1(jj,ii) = spectri(jj,ii)*cmplx(cos(theta),sin(theta))
end do
end do
! call dfftw_execute_dft_c2r(plan_back,spectric1,field1)
call FFT_SP(spectric1,field1)

RETURN
END SUBROUTINE RandomF

!*****************************************************************
SUBROUTINE energy(rho,ux,uy,bx,by,uxdx,uydy,bxdx,bydy,EU,EB,Erho,Erho2,divu,divb)
!***********compute energies
use parameters
use fftw_mod
implicit none
double precision EU, EB, Erho, Erho2, divu, divb
! double precision rho(N,N), ux(N,N), uy(N,N), bx(N,N), by(N,N)
! double precision uxdx(N,N), uydy(N,N), bxdx(N,N), bydy(N,N)
double complex rho(Nh,N), ux(Nh,N), uy(Nh,N), bx(Nh,N), by(Nh,N)
double complex uxdx(Nh,N), uydy(Nh,N), bxdx(Nh,N), bydy(Nh,N)
double precision uxdxr(N,N), uydyr(N,N), bxdxr(N,N), bydyr(N,N)
double precision norm, norm2
integer ii, jj

norm = 1./real(N*N*N)
norm2 = norm*norm
EU = 0.
EB = 0.
Erho = 0.
Erho2 = 0.
divu = 0.
divb = 0.
call FFT_SP(uxdx,uxdxr)
call FFT_SP(uydy,uydyr)
call FFT_SP(bxdx,bxdxr)
call FFT_SP(bydy,bydyr)

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
divu = divu + uxdxr(jj,ii) + uydyr(jj,ii)
divb = divb + bxdxr(jj,ii) + bydyr(jj,ii)
end do
end do

divu = divu*norm
divb = divb*norm

! do ii = 1, N
! do jj = 1, Nh
! EU = EU + 0.5d0*(abs(ux(jj,ii))**2 + abs(uy(jj,ii))**2)
! EB = EB + 0.5d0*(abs(bx(jj,ii))**2 + abs(by(jj,ii))**2)
! Erho = Erho + 0.5d0*(1.d0+rho(jj,ii))*(ux(jj,ii)**2 + uy(jj,ii)**2)
! Erho2 = Erho2 + rho(jj,ii)
! divu = divu + uxdxr(jj,ii) + uydyr(jj,ii)
! divb = divb + bxdxr(jj,ii) + bydyr(jj,ii)
! end do
! end do

RETURN
END SUBROUTINE energy

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

!*****************************************************************
!     Adaptive timestep
!*****************************************************************
! SUBROUTINE adaptiveT(a,b,ta)
! use parameters
! implicit none
! double precision a(N,N), b(N,N), ta, max1, max2, max3, max4, max
! integer ii, jj
! include "fftw3.f"

! max1 = dabs(a(1,1))
! do ii = 1, N
! do jj = 1, N
! max2 = dabs(a(jj,ii))
! if (max2 .gt. max1) then
! max1=max2
! end if
! end do
! end do

! max3 = dabs(b(1,1))
! do ii = 1, N
! do jj = 1, N
! max4 = dabs(b(jj,ii))
! if (max4 .gt. max3) then
! max3=max4
! end if
! end do
! end do

! max=max1
! if (max3 .gt. max1) then
! max=max3
! end if
! if (max .lt. 2.d0) then
! max=1.d0
! end if
! ta=1.d0/(N*max)

! RETURN
! END SUBROUTINE adaptiveT
!*****************************************************************




!*****************************************************************