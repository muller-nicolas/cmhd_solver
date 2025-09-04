!*********************************************************************
!  Numerical integration of 2D CMHD equations for F-wave turbulence
!  Sébastien Galtier - LPP - Version 8 (July 2025)
!*********************************************************************

! TODO: Declare global plans for forward and backward fft
! TODO: Dealisaing?
! TODO: RK4
! TODO: Parallelisation: MPI, openMP or GPU
! TODO: Save fields as binary files

Program CMHD

use, intrinsic :: iso_c_binding
use fftw_mod

implicit none
integer, parameter :: N = 256
integer, parameter :: Nh = N/2+1
integer, parameter :: Na = N/3  ! partie entière
double precision pi, dk, deltaT, deltaTi, nu, eta, Ek(Nh,N)
double precision kx(N), ky(Nh), kd(Nh,N), kinj, t1, t2
double precision divu, divb, ta, time, disp, norm
double precision rho0(N,N), ux0(N,N), uy0(N,N), bx0(N,N), by0(N,N)
double precision rho1(N,N), ux1(N,N), uy1(N,N), bx1(N,N), by1(N,N)
double precision rho2(N,N), ux2(N,N), uy2(N,N), bx2(N,N), by2(N,N)
double precision uxdx0(N,N), uxdy0(N,N), uxdx1(N,N), uxdy1(N,N)
double precision uydx0(N,N), uydy0(N,N), bxdy0(N,N), bydx0(N,N)
double precision uydx1(N,N), uydy1(N,N), bxdx1(N,N), bxdy1(N,N)
double precision bydx1(N,N), bydy1(N,N), nonlinby(N,N), nonlinbydx(N,N)
double precision nonlinuxa(N,N), nonlinuxadx(N,N), nonlinuxb(N,N)
double precision nonlinuxbdy(N,N), nonlinuy(N,N), nonlinuydx(N,N)
double precision nonlinrho0(N,N), nonlinux0(N,N), nonlinuy0(N,N)
double precision nonlinbx0(N,N), nonlinbx(N,N), nonlinbxdy(N,N)
double precision rhoux(N,N), rhouy(N,N), rhouxdx(N,N), rhouydy(N,N)
double precision nonlinrho1(N,N), nonlinux1(N,N), nonlinuy1(N,N)
double precision nonlinbx1(N,N), EU, EB, Erho, Erho2, a, amp, off
double precision nonlinby0(N,N), nonlinby1(N,N), nonlinuxbdx(N,N)
double complex :: ukx2(Nh,N), uky2(Nh,N), bkx2(Nh,N), bky2(Nh,N)
double complex :: rhok(Nh,N)

include "fftw3.f"

integer (kind=8) plan_for, plan_back
integer i, j, ndeltaT, inrj, ispec, it, istore, irestart, nrestart
character (len=11) :: animR='restart-'
character (len=21) :: animE='out_spectrumEU-2D-'
character (len=21) :: animB='out_spectrumEB-2D-'
character (len=22) :: animrho='out_spectrumrho-2D-'
character (len=14) :: animO='out_rho-2D-'
character (len=13) :: animW='out_wz-2D-'
character (len=13) :: animJ='out_jz-2D-'
character (len=13) :: animux='out_ux-2D-'
character (len=13) :: animuy='out_uy-2D-'
character (len=13) :: animbx='out_bx-2D-'
character (len=13) :: animby='out_by-2D-'
character (len=15) :: animdiv='out_divb-2D-'
character (len=15) :: animdivu='out_divu-2D-'

call cpu_time(time=t1)

!**************Initialization
istore = 100
pi = 3.141592653589793238d0
deltaT = 1.d-4
deltaTi = deltaT/10.d0
ndeltaT = 1000
inrj = 1
ispec = 500  !*********must be a multiple of inrj
irestart = 1000
kinj = 3.
dk = 2.*pi
! nu = 1.d-6
nu = 6.d-8
eta = nu
disp = 3.d-5  ! without dispersion => 0.d0
a = 1.d0  !*********a=0. => linear equations; a=1. non-linear equations
amp = 1.d-3
off = 0.d0   !*********forcing => off = 1.d0
time = 0.d0
nrestart = 1    !*********for a restart => nrestart = 0


call init_fftw(plan_for, plan_back, ux0, ukx2, N, Nh)
call Initk(kx,ky,kd,dk,kinj,Nh,N)

!***************** In case of no restart the code starts down here
if (nrestart .ne. 0) then

open(30, file='out_parameter', status='new', form='formatted')
write(30,*) deltaT, ndeltaT, inrj, kinj, ispec, N, dk
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

call RandomInit(ux0,uy0,kinj,pi,Na,Nh,N)
call energyF(ux0,uy0,EU,N)
ux0=amp*ux0/sqrt(EU)
uy0=amp*uy0/sqrt(EU)
call derivex(ux0,uxdx0,kx,Nh,N)
call derivex(uy0,uydx0,kx,Nh,N)
call derivex(by0,bydx0,kx,Nh,N)
call derivey(ux0,uxdy0,ky,Nh,N)
call derivey(uy0,uydy0,ky,Nh,N)
call derivey(bx0,bxdy0,ky,Nh,N)

rhoux = rho0*ux0
rhouy = rho0*uy0
call derivex(rhoux,rhouxdx,kx,Nh,N)
call derivey(rhouy,rhouydy,ky,Nh,N)
nonlinbx = ux0*by0-bx0*uy0
call derivey(nonlinbx,nonlinbxdy,ky,Nh,N)
nonlinby = -nonlinbx
call derivex(nonlinby,nonlinbydx,kx,Nh,N)
nonlinuxa = (bx0*bx0-by0*by0-ux0*ux0)/2.
call derivex(nonlinuxa,nonlinuxadx,kx,Nh,N)
nonlinuxb = bx0*by0
call derivey(nonlinuxb,nonlinuxbdy,ky,Nh,N)
nonlinuy = -(uy0*uy0+bx0*bx0-by0*by0)/2.
call derivey(nonlinuy,nonlinuydx,ky,Nh,N)
call derivex(nonlinuxb,nonlinuxbdx,kx,Nh,N)

nonlinrho0 = -uxdx0-uydy0 - a*(rhouxdx+rhouydy)
nonlinux0 = a*(nonlinuxadx + nonlinuxbdy - uy0*uxdy0)
nonlinuy0 = bydx0 - bxdy0 + a*(nonlinuydx+nonlinuxbdx-rho0*(bydx0-bxdy0))
nonlinbx0 = -uydy0 + a*nonlinbxdy
nonlinby0 = uydx0 + a*nonlinbydx
ux1 = ux0 + deltaTi*nonlinux0
uy1 = uy0 + deltaTi*nonlinuy0
bx1 = bx0 + deltaTi*nonlinbx0
by1 = by0 + deltaTi*nonlinby0
end if

!****************** In case of restart the code starts below *************
if (nrestart .eq. 0) then
open(66, file = 'restart', status = 'old',form='unformatted')
read(66) rho1(:,:),ux1(:,:),uy1(:,:),bx1(:,:),by1(:,:),rho0(:,:),ux0(:,:),uy0(:,:),bx0(:,:),by0(:,:)
read(66) nonlinrho0(:,:),nonlinux0(:,:),nonlinuy0(:,:),nonlinbx0(:,:),nonlinby0(:,:)
read(66) nonlinrho1(:,:),nonlinux1(:,:),nonlinuy1(:,:),nonlinbx1(:,:),nonlinby1(:,:)
close(66)
end if

!************************************************************************
!*****Time evolution: Adams-Bashforth
!*****************************  Main loop  ******************************
do it = 1, ndeltaT

! Random forcing in ux0 and uy0
call RandomF(ux0,kinj,pi,it,Na,Nh,N)
call RandomF(uy0,kinj,pi,it,Na,Nh,N)

! Normalization of the forcing
call energyF(ux0,uy0,EU,N)
ux0=amp*ux0/sqrt(EU)
uy0=amp*uy0/sqrt(EU)

! d/dx
call derivex(ux1,uxdx1,kx,Nh,N)
call derivex(uy1,uydx1,kx,Nh,N)
call derivex(bx1,bxdx1,kx,Nh,N)
call derivex(by1,bydx1,kx,Nh,N)

! d/dy
call derivey(ux1,uxdy1,ky,Nh,N)
call derivey(uy1,uydy1,ky,Nh,N)
call derivey(bx1,bxdy1,ky,Nh,N)
call derivey(by1,bydy1,ky,Nh,N)

! Compute non-linear terms

! NL terms for Density equation
rhoux = rho1*ux1
rhouy = rho1*uy1
call derivex(rhoux,rhouxdx,kx,Nh,N)
call derivey(rhouy,rhouydy,ky,Nh,N)

! NL terms for Induction equation
nonlinbx = ux1*by1-bx1*uy1
call derivey(nonlinbx,nonlinbxdy,ky,Nh,N)
nonlinby = -nonlinbx
call derivex(nonlinby,nonlinbydx,kx,Nh,N)

! NL terms for Velocity equation
nonlinuxa = (bx1*bx1-by1*by1-ux1*ux1)/2.
call derivex(nonlinuxa,nonlinuxadx,kx,Nh,N)
nonlinuxb = bx1*by1
call derivey(nonlinuxb,nonlinuxbdy,ky,Nh,N)
nonlinuy = -(uy1*uy1+bx1*bx1-by1*by1)/2.
call derivey(nonlinuy,nonlinuydx,ky,Nh,N)
call derivex(nonlinuxb,nonlinuxbdx,kx,Nh,N)

! NL terms all together
nonlinrho1 = -uxdx1-uydy1 - a*(rhouxdx+rhouydy)
nonlinux1 = a*(nonlinuxadx + nonlinuxbdy - uy1*uxdy1)
nonlinuy1 = bydx1 - bxdy1 + a*(nonlinuydx+nonlinuxbdx-rho1*(bydx1-bxdy1))
nonlinbx1 = -uydy1 + a*nonlinbxdy
nonlinby1 = uydx1 + a*nonlinbydx

! Adams-Bashford method
rho2 = rho1 + deltaT*(1.5*nonlinrho1 - 0.5*nonlinrho0)
ux2 = ux1 + deltaT*(1.5*nonlinux1 - 0.5*nonlinux0) + ux0*off
uy2 = uy1 + deltaT*(1.5*nonlinuy1 - 0.5*nonlinuy0) + uy0*off
bx2 = bx1 + deltaT*(1.5*nonlinbx1 - 0.5*nonlinbx0)
by2 = by1 + deltaT*(1.5*nonlinby1 - 0.5*nonlinby0)

! Implicit method for dissipation term in rho

norm = 1.d0/real(N*N)
call dfftw_execute_dft_r2c(plan_for,rho2,rhok)
do i = 1, N
    do j = 1, Nh
        rhok(j,i) = rhok(j,i)*exp(-(kd(j,i)*nu*kd(j,i)+cmplx(0,1)*disp*(kx(i)**3+ky(j)**3))*deltaT)
    end do
end do
call dfftw_execute_dft_c2r(plan_back,rhok,rho2)
rho2 = rho2*norm

! Implicit method for dissipation term in ux
call dfftw_execute_dft_r2c(plan_for,ux2,ukx2)
do i = 1, N
    do j = 1, Nh
        ukx2(j,i) = ukx2(j,i)*exp(-(kd(j,i)*nu*kd(j,i)+cmplx(0,1)*disp*(kx(i)**3+ky(j)**3))*deltaT)
    end do
end do
!ukx2(1,1)=0.d0
call dfftw_execute_dft_c2r(plan_back,ukx2,ux2)
ux2 = ux2*norm

! Implicit method for dissipation term in uy
call dfftw_execute_dft_r2c(plan_for,uy2,uky2)
do i = 1, N
    do j = 1, Nh
        uky2(j,i) = uky2(j,i)*exp(-(kd(j,i)*nu*kd(j,i)+cmplx(0,1)*disp*(kx(i)**3+ky(j)**3))*deltaT)
    end do
end do
!uky2(1,1)=0.d0
call dfftw_execute_dft_c2r(plan_back,uky2,uy2)
uy2 = uy2*norm

! Implicit method for dissipation term in bx
call dfftw_execute_dft_r2c(plan_for,bx2,bkx2)
do i = 1, N
    do j = 1, Nh
        bkx2(j,i) = bkx2(j,i)*exp(-(kd(j,i)*nu*kd(j,i)+cmplx(0,1)*disp*(kx(i)**3+ky(j)**3))*deltaT)
    end do
end do
call dfftw_execute_dft_c2r(plan_back,bkx2,bx2)
bx2 = bx2*norm

! Implicit method for dissipation term in by
call dfftw_execute_dft_r2c(plan_for,by2,bky2)
do i = 1, N
    do j = 1, Nh
        bky2(j,i) = bky2(j,i)*exp(-(kd(j,i)*nu*kd(j,i)+cmplx(0,1)*disp*(kx(i)**3+ky(j)**3))*deltaT)
    end do
end do
call dfftw_execute_dft_c2r(plan_back,bky2,by2)
by2 = by2*norm

! Rename variables for saving and use them as initial values for next loop
rho1=rho2
ux1=ux2
uy1=uy2
bx1=bx2
by1=by2
nonlinrho0=nonlinrho1
nonlinux0=nonlinux1
nonlinuy0=nonlinuy1
nonlinbx0=nonlinbx1
nonlinby0=nonlinby1

! Compute and write energy
if ( (mod(it,inrj) .eq. 0) ) then
call energy(rho1,ux1,uy1,bx1,by1,uxdx1,uydy1,bxdx1,bydy1,EU,EB,Erho,Erho2,divu,divb,N)
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
call adaptiveT(ux2,rho2,ta,N)
write(51,*) ta
deltaT = ta*0.2d0 ! condition CFL
nu = 1./ta*(1.d0/N)**4 ! condition CFL
eta = nu

end if

! Compute and write spectra
if ( (mod(it,ispec) .eq. 0) ) then
print *, it
write(animE(19:21),'(i3)') istore
call spectrum(ux2,uy2,Ek,Na,Nh,N)
open(30, file = animE, status = 'new',form='formatted')
do i = 1, Nh
    do j = 1, N
        write(30,*) Ek(i,j)
    end do
end do
close(30)
write(animB(19:21),'(i3)') istore
call spectrum(bx2,by2,Ek,Na,Nh,N)
open(30, file = animB, status = 'new',form='formatted')
do i = 1, Nh
    do j = 1, N
        write(30,*) Ek(i,j)
    end do
end do
close(30)
write(animrho(20:22),'(i3)') istore
call spectrumrho(rho2,Ek,Na,Nh,N)
open(30, file = animrho, status = 'new',form='formatted')
do i = 1, Nh
    do j = 1, N
        write(30,*) Ek(i,j)
    end do
end do
close(30)
!
write(animO(12:14),'(i3)') istore
open(30, file = animO, status = 'new',form='formatted')
write(30,*) rho2(:,:)
close(30)
!write(animux(11:13),'(i3)') istore
!open(30, file = animux, status = 'new',form='formatted')
!write(30,*) ux2(:,:)
!close(30)
!write(animuy(11:13),'(i3)') istore
!open(30, file = animuy, status = 'new',form='formatted')
!write(30,*) uy2(:,:)
!close(30)
!write(animbx(11:13),'(i3)') istore
!open(30, file = animbx, status = 'new',form='formatted')
!write(30,*) bx2(:,:)
!close(30)
!write(animby(11:13),'(i3)') istore
!open(30, file = animby, status = 'new',form='formatted')
!write(30,*) by2(:,:)
!close(30)
!
write(animW(11:13),'(i3)') istore
open(30, file = animW, status = 'new',form='formatted')
write(30,*) uydx1(:,:)-uxdy1(:,:)
close(30)
write(animJ(11:13),'(i3)') istore
open(30, file = animJ, status = 'new',form='formatted')
write(30,*) bydx1(:,:)-bxdy1(:,:)
close(30)
!
write(animdiv(13:15),'(i3)') istore
open(30, file = animdiv, status = 'new',form='formatted')
write(30,*) bxdx1(:,:)+bydy1(:,:)
close(30)
write(animdivu(13:15),'(i3)') istore
open(30, file = animdivu, status = 'new',form='formatted')
write(30,*) uxdx1(:,:)+uydy1(:,:)
close(30)
istore = istore + 1
end if

end do ! end of the temporal loop

! call end_fftw(plan_for, plan_back)

write(animR(9:11),'(i3)') istore
open(30, file = animR, status = 'new',form='unformatted')
write(30) rho1(:,:),ux1(:,:),uy1(:,:),bx1(:,:),by1(:,:),rho0(:,:),ux0(:,:),uy0(:,:),bx0(:,:),by0(:,:)
write(30) nonlinrho0(:,:),nonlinux0(:,:),nonlinuy0(:,:),nonlinbx0(:,:),nonlinby0(:,:)
write(30) nonlinrho1(:,:),nonlinux1(:,:),nonlinuy1(:,:),nonlinbx1(:,:),nonlinby1(:,:)
close(30)

close(40)
close(41)
call cpu_time(time=t2)
write(*,*) "cpu time", t2-t1

print *, 'OK'
end program CMHD

!*****************************************************************
SUBROUTINE Initk(kkx,kky,kkkd,ddk,kinj,Nh,N)
! Initialization of the wavenumbers
implicit none
double precision kkx(N), kky(Nh), kkkd(Nh,N), ddk, kinj
integer Nh, N, ii

kky=(/(dfloat(ii-1)*ddk,ii=1,Nh,1)/)
kkx(1:N/2)=(/(dfloat(ii-1)*ddk,ii=1,N/2,1)/)
kkx(N/2+1:N)=(/(dfloat(ii-1-N)*ddk,ii=N/2+1,N,1)/)
do ii = 1, N
kkkd(:,ii) = kkx(ii)*kkx(ii) + kky(:)*kky(:)
end do

RETURN
END SUBROUTINE Initk

!*****************************************************************
SUBROUTINE RandomInit(uxi,uyi,kinj,pi,Na,Nh,N)
! Initial random field spectra
implicit none
integer,parameter :: seed = 800
double precision spectri(Na+1,Na+1), uxi(N,N), uyi(N,N)
double precision pi, theta, knc, kinj
double complex :: spectric1(Nh,N), spectric2(Nh,N)
integer Na, N, Nh, ii, jj
integer (kind=8) plan_for, plan_back
include "fftw3.f"

call srand(seed)
uxi=0.
uyi=0.
spectri=0.
spectric1=0.
spectric2=0.
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
theta = rand()*2.*pi
spectric2(jj,ii) = spectri(jj,ii)*cmplx(cos(theta),sin(theta))
end do
end do
call dfftw_plan_dft_c2r_2d_(plan_back,N,N,spectric1,uxi,FFTW_ESTIMATE)
call dfftw_execute_(plan_back)
! call dfftw_destroy_plan(plan_back)
uxi=uxi/real(N*N)
call dfftw_plan_dft_c2r_2d_(plan_back,N,N,spectric2,uyi,FFTW_ESTIMATE)
call dfftw_execute_(plan_back)
call dfftw_destroy_plan(plan_back)
uyi=uyi/real(N*N)

RETURN
END SUBROUTINE RandomInit

!*****************************************************************
! Random forcing spectrum
!*****************************************************************
SUBROUTINE RandomF(field1,kinj,pi,seed,Na,Nh,N)
implicit none
double precision spectri(Na+1,Na+1), field1(N,N)
double precision pi, theta, knc, kinj
double complex :: spectric1(Nh,N)
integer N, Nh, Na, ii, jj, seed
integer (kind=8) plan_for, plan_back
include "fftw3.f"

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
call dfftw_plan_dft_c2r_2d_(plan_back,N,N,spectric1,field1,FFTW_ESTIMATE)
call dfftw_execute_(plan_back)
call dfftw_destroy_plan(plan_back)

RETURN
END SUBROUTINE RandomF

!*****************************************************************
SUBROUTINE derivex(aa,bb,kkx,Nh,N)
! Computation of the x-derivative
implicit none
complex*16, parameter :: imag = (0.0d0,1.0d0)
double precision aa(N,N), bb(N,N), kkx(N)
double complex :: cc(Nh,N)
integer (kind=8) plan_for, plan_back
integer Nh, N, ii
include "fftw3.f"

call dfftw_plan_dft_r2c_2d_(plan_for,N,N,aa,cc,FFTW_ESTIMATE)
call dfftw_execute_(plan_for)
! call dfftw_execute_dft_r2c(plan_for,aa,cc)
do ii = 1, N
cc(:,ii) = imag*cc(:,ii)*kkx(ii)
end do
call dfftw_plan_dft_c2r_2d_(plan_back,N,N,cc,bb,FFTW_ESTIMATE)
! call dfftw_execute_dft_c2r(plan_back,cc,bb)
call dfftw_execute_(plan_back)
call dfftw_destroy_plan(plan_back)
call dfftw_destroy_plan(plan_for)
bb=bb/real(N*N)

RETURN
END SUBROUTINE derivex

!*****************************************************************
SUBROUTINE derivey(aa,bb,kky,Nh,N)
! Computation of the y-derivative
implicit none
complex*16, parameter :: imag = (0.0d0,1.0d0)
double precision aa(N,N), bb(N,N), kky(Nh)
double complex :: cc(Nh,N)
integer (kind=8) plan_for, plan_back
integer Nh, N, jj
include "fftw3.f"

call dfftw_plan_dft_r2c_2d_(plan_for,N,N,aa,cc,FFTW_ESTIMATE)
call dfftw_execute_(plan_for)
do jj = 1, Nh
cc(jj,:) = imag*cc(jj,:)*kky(jj)
end do
call dfftw_plan_dft_c2r_2d_(plan_back,N,N,cc,bb,FFTW_ESTIMATE)
call dfftw_execute_(plan_back)
call dfftw_destroy_plan(plan_back)
call dfftw_destroy_plan(plan_for)
bb=bb/real(N*N)

RETURN
END SUBROUTINE derivey

!*****************************************************************
SUBROUTINE energy(rho,ux,uy,bx,by,uxdx,uydy,bxdx,bydy,EU,EB,Erho,Erho2,divu,divb,N)
!***********compute energies
implicit none
double precision EU, EB, Erho, Erho2, divu, divb
double precision rho(N,N), ux(N,N), uy(N,N), bx(N,N), by(N,N)
double precision uxdx(N,N), uydy(N,N), bxdx(N,N), bydy(N,N)
double precision norm
integer N, iia, iib
include "fftw3.f"

EU = 0.
EB = 0.
Erho = 0.
Erho2 = 0.
divu = 0.
divb = 0.
do iia = 1, N
do iib = 1, N
EU = EU + (ux(iia,iib)**2 + uy(iia,iib)**2)/2.d0
EB = EB + (bx(iia,iib)**2 + by(iia,iib)**2)/2.d0
Erho = Erho + (1.d0+rho(iia,iib))*(ux(iia,iib)**2 + uy(iia,iib)**2)/2.d0
divu = divu + uxdx(iia,iib) + uydy(iia,iib)
divb = divb + bxdx(iia,iib) + bydy(iia,iib)
Erho2 = Erho2 + rho(iia,iib)
end do
end do
norm = 1.d0/real(N*N)
EU = EU*norm
EB = EB*norm
Erho = Erho*norm
Erho2 = Erho2*norm
divu = divu*norm
divb = divb*norm

RETURN
END SUBROUTINE energy

!*****************************************************************
SUBROUTINE energyF(ux,uy,EU,N)
!***********compute energy
implicit none
double precision EU, ux(N,N), uy(N,N)
integer N, ii, jj
include "fftw3.f"

EU = 0.
do ii = 1, N
do jj = 1, N
EU = EU + (ux(ii,jj)**2 + uy(ii,jj)**2)
end do
end do
EU = EU/real(N*N)

RETURN
END SUBROUTINE energyF

!*****************************************************************
SUBROUTINE spectrum(ux,uy,Ek,Na,Nh,N)
!***********compute the 2D spectrum
implicit none
double precision ux(N,N), uy(N,N), Ek(Nh,N)
double complex :: ukx(Nh,N), uky(Nh,N), Ek1(Nh,N)
integer (kind=8) plan_for, plan_back
integer Nh, N, Na, iia, iib
include "fftw3.f"

call dfftw_plan_dft_r2c_2d_(plan_for,N,N,ux,ukx,FFTW_ESTIMATE)
call dfftw_execute_(plan_for)
call dfftw_plan_dft_r2c_2d_(plan_for,N,N,uy,uky,FFTW_ESTIMATE)
call dfftw_execute_(plan_for)
Ek1 = abs(ukx)**2 + abs(uky)**2
do iia = 1, Nh
do iib = 1, N
Ek(iia,iib) = real(Ek1(iia,iib))
end do
end do
call dfftw_destroy_plan(plan_for)

RETURN
END SUBROUTINE spectrum

!*****************************************************************
SUBROUTINE spectrumrho(rho,Ek,Na,Nh,N)
!***********compute the 2D spectrum
implicit none
double precision rho(N,N), Ek(Nh,N)
double complex :: rhok(Nh,N), Ek1(Nh,N)
integer (kind=8) plan_for, plan_back
integer Nh, N, Na, iia, iib
include "fftw3.f"

call dfftw_plan_dft_r2c_2d_(plan_for,N,N,rho,rhok,FFTW_ESTIMATE)
call dfftw_execute_(plan_for)
Ek1 = abs(rhok)**2
do iia = 1, Nh
do iib = 1, N
Ek(iia,iib) = real(Ek1(iia,iib))
end do
end do
call dfftw_destroy_plan(plan_for)

RETURN
END SUBROUTINE spectrumrho

!*****************************************************************
!     Adaptive timestep
!*****************************************************************
SUBROUTINE adaptiveT(a,b,ta,N)
implicit none
double precision a(N,N), b(N,N), ta, max1, max2, max3, max4, max
integer N, ii, jj
include "fftw3.f"

max1 = dabs(a(1,1))
do ii = 1, N
do jj = 1, N
max2 = dabs(a(ii,jj))
if (max2 .gt. max1) then
max1=max2
end if
end do
end do

max3 = dabs(b(1,1))
do ii = 1, N
do jj = 1, N
max4 = dabs(b(ii,jj))
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
!*****************************************************************
