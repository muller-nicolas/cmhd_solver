!*********************************************************************
!  Numerical integration of 2D CMHD equations for F-wave turbulence
!  Sébastien Galtier - LPP - Version 8 (July 2025)
!*********************************************************************

! TODO: outputs as functions
! TODO: Save in real time, and not in the end
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
use outputs

implicit none

double precision Ek(Nh,N)
double precision t1, t2
double precision divu, divb, ta, norm, time, timests

double precision :: rho0(N,N), ux0(N,N), uy0(N,N), bx0(N,N), by0(N,N)
double precision :: rho1(N,N), ux1(N,N), uy1(N,N), bx1(N,N), by1(N,N)
double precision :: uxdx(N,N), uxdy(N,N), uydx(N,N), uydy(N,N)
double precision :: bxdy(N,N), bydx(N,N), bxdx(N,N), bydy(N,N)
! double precision :: nonlinrho0(N,N), nonlinux0(N,N), nonlinuy0(N,N), nonlinbx0(N,N), nonlinby0(N,N)
! double precision :: nonlinrho1(N,N), nonlinux1(N,N), nonlinuy1(N,N), nonlinbx1(N,N), nonlinby1(N,N)

double precision EU, EB, Erho, Erho2

double complex :: ukx0(Nh,N), uky0(Nh,N), bkx0(Nh,N), bky0(Nh,N), rhok0(Nh,N)
double complex :: ukx1(Nh,N), uky1(Nh,N), bkx1(Nh,N), bky1(Nh,N), rhok1(Nh,N)
double complex :: ukx2(Nh,N), uky2(Nh,N), bkx2(Nh,N), bky2(Nh,N), rhok2(Nh,N)
double complex :: nonlinrhok0(Nh,N), nonlinukx0(Nh,N), nonlinuky0(Nh,N), nonlinbkx0(Nh,N), nonlinbky0(Nh,N)
double complex :: nonlinrhok1(Nh,N), nonlinukx1(Nh,N), nonlinuky1(Nh,N), nonlinbkx1(Nh,N), nonlinbky1(Nh,N)

integer i, j, it 
character (len=11) :: animR='restart-'
! character (len=21) :: animE='out_spectrumEU-2D-'
! character (len=21) :: animB='out_spectrumEB-2D-'
! character (len=22) :: animrho='out_spectrumrho-2D-'
! character (len=14) :: animO='out_rho-2D-'
! character (len=13) :: animW='out_wz-2D-'
! character (len=13) :: animJ='out_jz-2D-'
! character (len=15) :: animdiv='out_divb-2D-'
! character (len=15) :: animdivu='out_divu-2D-'

call cpu_time(time=t1)

!**************Initialization

deltaT = deltaT0
nu = nu0
time = 0.d0
timests = 0.d0

call init_fftw
call Initk

!***************** In case of no restart the code starts down here
if (nrestart .ne. 0) then

open(30, file='out_parameter', status='new', form='formatted')
write(30,*) deltaT, ndeltaT, inrj, kinj, ispec, ifields, N, dk
close(30)

open(51, file = 'out_deltaT', status = 'new',form='formatted')
! open(52, file = 'out_time', status = 'new',form='formatted')
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
read(66) rhok1(:,:),ukx1(:,:),uky1(:,:),bkx1(:,:),bky1(:,:)
read(66) rhok0(:,:),ukx0(:,:),uky0(:,:),bkx0(:,:),bky0(:,:)
read(66) nonlinrhok0(:,:),nonlinukx0(:,:),nonlinuky0(:,:),nonlinbkx0(:,:),nonlinbky0(:,:)
read(66) nonlinrhok1(:,:),nonlinukx1(:,:),nonlinuky1(:,:),nonlinbkx1(:,:),nonlinbky1(:,:)
close(66)
end if

!************************************************************************
!*****Time evolution: Adams-Bashforth
!*****************************  Main loop  ******************************
do it = 1, ndeltaT

! Random forcing in ux0 and uy0
call RandomF(ukx0)
call RandomF(uky0)

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

call save_energy(rhok2,ukx2,uky2,bkx2,bky2)

open(52, file='out_time', position='append',form='formatted')
time = time + dfloat(inrj)*deltaT
write(52,*) time
close(52)

!****************Adaptive timestep
! write(53,*) nu
! call adaptiveT(ukx2,rhok2,ta)
! write(51,*) ta
! deltaT = ta*0.2d0 ! condition CFL
! nu = 1./ta*(1.d0/N)**4 ! condition CFL
! eta = nu

end if

! Compute and write spectra
if ( (mod(it,ispec) .eq. 0) ) then
print *, it
call save_spectra(rhok2,ukx2,uky2,bkx2,bky2,istore_sp)
endif

! Write fields
! rho, wz, jz, divu, divb
if (mod(it,ifields) .eq. 0) then
call save_fields(rhok2,ukx2,uky2,bkx2,bky2,istore_fields)
end if

! Save spatio-temporal spectra
if (sts .eq. 1) then
    if (mod(it,ists) .eq. 0) then
        timests = timests + dfloat(ists)*deltaT
        call WriteSpatioTemporalSpectrum(ukx1, uky1, bkx1, bky1, rhok1, timests)
    end if
end if


end do ! end of the temporal loop


! Save fields in spectral space for restart
write(animR(9:11),'(i3)') istore_fields
open(30, file = animR, status = 'new',form='unformatted')
write(30) rhok1(:,:),ukx1(:,:),uky1(:,:),bkx1(:,:),bky1(:,:)
write(30) rhok0(:,:),ukx0(:,:),uky0(:,:),bkx0(:,:),bky0(:,:)
write(30) nonlinrhok0(:,:),nonlinukx0(:,:),nonlinuky0(:,:),nonlinbkx0(:,:),nonlinbky0(:,:)
write(30) nonlinrhok1(:,:),nonlinukx1(:,:),nonlinuky1(:,:),nonlinbkx1(:,:),nonlinbky1(:,:)
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
double precision spectri(Na+1,Na+1)
double precision theta, knc
double complex :: field1(Nh,N)
integer ii, jj

call srand(seed)
spectri=0.
field1=0.
do ii = 1, Na+1
do jj = 1, Na+1
knc = (kinj - real(ii-1) - real(jj-1))**4
spectri(jj,ii) = dexp(-knc*100.)
end do
end do
do ii = 1, Na+1
do jj = 1, Na+1
theta = rand()*2.*pi
field1(jj,ii) = spectri(jj,ii)*(cos(theta) + imag*sin(theta))
end do
end do

RETURN
END SUBROUTINE RandomF



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




!*****************************************************************