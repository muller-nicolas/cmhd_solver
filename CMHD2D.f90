!*********************************************************************
!  Numerical integration of 2D CMHD equations for F-wave turbulence
!  Sébastien Galtier - LPP - Version 8 (July 2025)
!*********************************************************************

! TODO: correct time at restart
! TODO: Computation of flux
! TODO: Compare implicit vs explicit dissipation
! TODO: Read parameters from file after compilation
! TODO: Pointers for fields in memory (Reduce memory usage)
! TODO: write also forcing parameters in file
! TODO: handle precision globally (for practical purposes)
! TODO: Parallelisation: GPU

Program CMHD

use, intrinsic :: iso_c_binding
use omp_lib
use parameters
use fftw_mod
use spectral_mod
use cMHD_mod
use outputs

implicit none

double precision t1, t2
double precision ta, time, timests, phase
character(len=256) :: line, lastline
integer :: ios

double precision, dimension(:,:), allocatable :: ux1, uy1, bx1, by1, rho1
double complex, dimension(:,:), allocatable :: ukx1, uky1, bkx1, bky1, rhok1
double complex, dimension(:,:), allocatable :: nonlinrhok0, nonlinukx0, nonlinuky0, nonlinbkx0, nonlinbky0
double complex, dimension(:,:), allocatable :: nonlinrhok1, nonlinukx1, nonlinuky1, nonlinbkx1, nonlinbky1
double complex, dimension(:,:), allocatable :: fukx, fuky, fbkx, fbky

! double complex, dimension(:,:), allocatable :: ukx2, uky2, bkx2, bky2, rhok2 ! RK2

! double complex, dimension(:,:), allocatable :: ukx3, uky3, bkx3, bky3, rhok3 ! RK4
! double complex, dimension(:,:), allocatable :: ukx4, uky4, bkx4, bky4, rhok4 ! RK4
! double complex, dimension(:,:), allocatable :: nonlinrhok2, nonlinukx2, nonlinuky2, nonlinbkx2, nonlinbky2 ! RK4
! double complex, dimension(:,:), allocatable :: nonlinrhok3, nonlinukx3, nonlinuky3, nonlinbkx3, nonlinbky3 ! RK4

character(len=30) :: fname
character(len=12) :: numstr
integer i, j, it, corr

allocate(ux1(N,N), uy1(N,N), bx1(N,N), by1(N,N), rho1(N,N))
allocate(ukx1(Nh,N), uky1(Nh,N), bkx1(Nh,N), bky1(Nh,N), rhok1(Nh,N))
allocate(nonlinrhok0(Nh,N), nonlinukx0(Nh,N), nonlinuky0(Nh,N), nonlinbkx0(Nh,N), nonlinbky0(Nh,N))
allocate(nonlinrhok1(Nh,N), nonlinukx1(Nh,N), nonlinuky1(Nh,N), nonlinbkx1(Nh,N), nonlinbky1(Nh,N))
allocate(fukx(Nh,N), fuky(Nh,N), fbkx(Nh,N), fbky(Nh,N))

! allocate(ukx2(Nh,N), uky2(Nh,N), bkx2(Nh,N), bky2(Nh,N), rhok2(Nh,N)) ! RK2

! allocate(ukx3(Nh,N), uky3(Nh,N), bkx3(Nh,N), bky3(Nh,N), rhok3(Nh,N)) ! RK4
! allocate(ukx4(Nh,N), uky4(Nh,N), bkx4(Nh,N), bky4(Nh,N), rhok4(Nh,N)) ! RK4
! allocate(nonlinrhok2(Nh,N), nonlinukx2(Nh,N), nonlinuky2(Nh,N), nonlinbkx2(Nh,N), nonlinbky2(Nh,N)) ! RK4
! allocate(nonlinrhok3(Nh,N), nonlinukx3(Nh,N), nonlinuky3(Nh,N), nonlinbkx3(Nh,N), nonlinbky3(Nh,N)) ! RK4

! nthreads = omp_get_max_threads()
call omp_set_num_threads(nthreads)
print *,"Using",nthreads,"threads"

t1 = omp_get_wtime()

!**************Initialization

deltaT = deltaT0
nu = nu0
time = 0.d0
timests = 0.d0

call init_fftw
call Initk

!***************** In case of no restart the code starts down here
if (nrestart .eq. 0) then

    open(30, file='out_parameter', status='new', form='formatted')
    write(30,'(E15.8,1X,I10,1X,I10,1X,E15.8,1X,I10,1X,I10,1X,I10,1X,E15.8)') &
       deltaT, ndeltaT, inrj, kinj, ispec, ifields, N, dk
    close(30)
    open(31, file='out_dissipation', status='new', form='formatted')
    write(31,'(E15.8,1X,E15.8,1X,E15.8,1X,E15.8)') &
       nu0, eta, alpha, disp
    close(31)

    ! Initial conditions
    ! Initilize velocity field
    call RandomInit(ukx1,uky1)
    ! call GaussianInit(ukx1,uky1)

!****************** In case of restart the code starts below *************
elseif (nrestart .ne. 0) then

    write(numstr, '(I5)') nrestart
    fname = 'field_rho-2D-' // trim(adjustl(numstr))
    open(66, file = fname, status = 'old',form='unformatted', access='stream')
    read(66) rho1(:,:)
    close(66)
    fname = 'field_ux-2D-' // trim(adjustl(numstr))
    open(66, file = fname, status = 'old',form='unformatted', access='stream')
    read(66) ux1(:,:)
    close(66)
    fname = 'field_uy-2D-' // trim(adjustl(numstr))
    open(66, file = fname, status = 'old',form='unformatted', access='stream')
    read(66) uy1(:,:)
    close(66)
    fname = 'field_bx-2D-' // trim(adjustl(numstr))
    open(66, file = fname, status = 'old',form='unformatted', access='stream')
    read(66) bx1(:,:)
    close(66)
    fname = 'field_by-2D-' // trim(adjustl(numstr))
    open(66, file = fname, status = 'old',form='unformatted', access='stream')
    read(66) by1(:,:)
    close(66)

    open(unit=66, file='out_time', status='old', action='read')
    do
        read(66, '(A)', iostat=ios) line
        if (ios /= 0) exit
        lastline = trim(line)
    end do
    close(66)
    read(lastline, *) time

    if (sts.eq.1) then
        open(unit=66, file='STS_time', status='old', action='read')
        do
            read(66, '(A)', iostat=ios) line
            if (ios /= 0) exit
            lastline = trim(line)
        end do
        close(66)
        read(lastline, *) timests
    endif

    call FFT_PS(rho1,rhok1)
    call FFT_PS(ux1,ukx1)
    call FFT_PS(uy1,uky1)
    call FFT_PS(bx1,bkx1)
    call FFT_PS(by1,bky1)

    rhok1 = kill*rhok1
    ukx1 = kill*ukx1
    uky1 = kill*uky1
    bkx1 = kill*bkx1
    bky1 = kill*bky1

    istore_sp = nrestart+1
    istore_fields = nrestart+1
end if

! Initialize forcing in ux0 and uy0
! call GaussianF(fukx,fuky)
! call RandomF(fukx,fuky)
call CompressibleRandomF(fukx,fuky,CompAmp)
! call RandomF_Ak(fukx,fuky,fbkx,fbky)
! call RandomF(fbkx,fbky)
! fbkx = 0.0
! fbky = 0.0
! call PoloidalRandomF(fukx,fuky)
corr = min(int(corr0/deltaT),1)

!! Do one iteration of the time-stepping for AB2
call RHS(rhok1,ukx1,uky1,bkx1,bky1,nonlinrhok0,nonlinukx0,nonlinuky0,nonlinbkx0,nonlinbky0)

nonlinukx0 = nonlinukx0 + fukx
nonlinuky0 = nonlinuky0 + fuky
nonlinbkx0 = nonlinbkx0 + fbkx 
nonlinbky0 = nonlinbky0 + fbky
rhok1 = rhok1 + deltaT*nonlinrhok0
ukx1  = ukx1  + deltaT*nonlinukx0
uky1  = uky1  + deltaT*nonlinuky0
bkx1  = bkx1  + deltaT*nonlinbkx0
bky1  = bky1  + deltaT*nonlinbky0


!************************************************************************
!*****Time evolution: Adams-Bashforth
!*****************************  Main loop  ******************************
do it = 1, ndeltaT

! Random phase for the forcing
if (mod(it,corr).eq.0) then
    !$omp parallel do private(i,j,phase)
    do i=1,N
        do j=1,Nh
            call random_number(phase)
            phase = 2*pi*phase
            fukx(j,i) = fukx(j,i) * (cos(phase) + imag*sin(phase))
            call random_number(phase)
            phase = 2*pi*phase
            fuky(j,i) = fuky(j,i) * (cos(phase) + imag*sin(phase))
            call random_number(phase)
            phase = 2*pi*phase
            fbkx(j,i) = fbkx(j,i) * (cos(phase) + imag*sin(phase))
            call random_number(phase)
            phase = 2*pi*phase
            fbky(j,i) = fbky(j,i) * (cos(phase) + imag*sin(phase))
        end do
    end do
end if

! Check evolution
call check_nan(rhok1) 

!!! Runge-Kutta 4 (RK4)
! call RHS(rhok0,ukx0,uky0,bkx0,bky0,nonlinrhok0,nonlinukx0,nonlinuky0,nonlinbkx0,nonlinbky0)
! rhok1 = rhok0 + deltaT*nonlinrhok0
! ukx1  = ukx0  + deltaT*nonlinukx0
! uky1  = uky0  + deltaT*nonlinuky0
! bkx1  = bkx0  + deltaT*nonlinbkx0
! bky1  = bky0  + deltaT*nonlinbky0

! call RHS(rhok1,ukx1,uky1,bkx1,bky1,nonlinrhok1,nonlinukx1,nonlinuky1,nonlinbkx1,nonlinbky1)
! rhok3 = rhok0 + 0.5*deltaT*(nonlinrhok1)
! ukx3  = ukx0  + 0.5*deltaT*(nonlinukx1)
! uky3  = uky0  + 0.5*deltaT*(nonlinuky1)
! bkx3  = bkx0  + 0.5*deltaT*(nonlinbkx1)
! bky3  = bky0  + 0.5*deltaT*(nonlinbky1)

! call RHS(rhok3,ukx3,uky3,bkx3,bky3,nonlinrhok3,nonlinukx3,nonlinuky3,nonlinbkx3,nonlinbky3)
! rhok4 = rhok0 + 0.5*deltaT*(nonlinrhok3)
! ukx4  = ukx0  + 0.5*deltaT*(nonlinukx3)
! uky4  = uky0  + 0.5*deltaT*(nonlinuky3)
! bkx4  = bkx0  + 0.5*deltaT*(nonlinbkx3)
! bky4  = bky0  + 0.5*deltaT*(nonlinbky3)

! call RHS(rhok4,ukx4,uky4,bkx4,bky4,nonlinrhok2,nonlinukx2,nonlinuky2,nonlinbkx2,nonlinbky2)
! rhok2 = rhok0 + deltaT/6.*(nonlinrhok0 + 2*nonlinrhok1 + 2*nonlinrhok3 + nonlinrhok2)
! ukx2  = ukx0  + deltaT/6.*( nonlinukx0 + 2*nonlinukx1 + 2*nonlinukx3 + nonlinukx2)
! uky2  = uky0  + deltaT/6.*( nonlinuky0 + 2*nonlinuky1 + 2*nonlinuky3 + nonlinuky2)
! bkx2  = bkx0  + deltaT/6.*( nonlinbkx0 + 2*nonlinbkx1 + 2*nonlinbkx3 + nonlinbkx2)
! bky2  = bky0  + deltaT/6.*( nonlinbky0 + 2*nonlinbky1 + 2*nonlinbky3 + nonlinbky2)

! ! Rename variables for saving and use them as initial values for next loop (RK4)
! rhok1=rhok2
! ukx1=ukx2
! uky1=uky2
! bkx1=bkx2
! bky1=bky2




! !!! Heun method (RK2)
! call RHS(rhok1,ukx1,uky1,bkx1,bky1,nonlinrhok0,nonlinukx0,nonlinuky0,nonlinbkx0,nonlinbky0)
! rhok2 = rhok1 + deltat*nonlinrhok0
! ukx2  = ukx1  + deltat*(nonlinukx0 + fukx)
! uky2  = uky1  + deltat*(nonlinuky0 + fuky)
! bkx2  = bkx1  + deltat*nonlinbkx0
! bky2  = bky1  + deltat*nonlinbky0
! ! 
! call rhs(rhok2,ukx2,uky2,bkx2,bky2,nonlinrhok1,nonlinukx1,nonlinuky1,nonlinbkx1,nonlinbky1)
! rhok2 = rhok1 + 0.5*deltat*(nonlinrhok0 + nonlinrhok1)
! ukx2  = ukx1  + 0.5*deltat*( nonlinukx0 + nonlinukx1 + fukx)
! uky2  = uky1  + 0.5*deltat*( nonlinuky0 + nonlinuky1 + fuky)
! bkx2  = bkx1  + 0.5*deltat*( nonlinbkx0 + nonlinbkx1)
! bky2  = bky1  + 0.5*deltat*( nonlinbky0 + nonlinbky1)

! ! ! rename variables for saving and use them as initial values for next loop (RK2)
! rhok1=rhok2
! ukx1=ukx2
! uky1=uky2
! bkx1=bkx2
! bky1=bky2



!! Adams-Bashforth method (AB2)
call RHS(rhok1,ukx1,uky1,bkx1,bky1,nonlinrhok1,nonlinukx1,nonlinuky1,nonlinbkx1,nonlinbky1)

! !$omp parallel do private(i,j)
! do i=1,N
!     do j=1,Nh
!         nonlinukx1(j,i) = nonlinukx1(j,i) + fukx(j,i)
!         nonlinuky1(j,i) = nonlinuky1(j,i) + fuky(j,i)
!         nonlinbkx1(j,i) = nonlinbkx1(j,i) + fbkx(j,i)
!         nonlinbky1(j,i) = nonlinbky1(j,i) + fbky(j,i)

!         rhok1(j,i) = rhok1(j,i) + deltaT*(1.5*nonlinrhok1(j,i) - 0.5*nonlinrhok0(j,i))
!         ukx1(j,i)  = ukx1(j,i)  + deltaT*(1.5*nonlinukx1(j,i)  - 0.5*nonlinukx0(j,i))
!         uky1(j,i)  = uky1(j,i)  + deltaT*(1.5*nonlinuky1(j,i)  - 0.5*nonlinuky0(j,i))
!         bkx1(j,i)  = bkx1(j,i)  + deltaT*(1.5*nonlinbkx1(j,i)  - 0.5*nonlinbkx0(j,i))
!         bky1(j,i)  = bky1(j,i)  + deltaT*(1.5*nonlinbky1(j,i)  - 0.5*nonlinbky0(j,i))

!         ! Rename variables for saving and use them as initial values for next loop (AB2)
!         nonlinrhok0(j,i)=nonlinrhok1(j,i)
!         nonlinukx0(j,i)=nonlinukx1(j,i)
!         nonlinuky0(j,i)=nonlinuky1(j,i)
!         nonlinbkx0(j,i)=nonlinbkx1(j,i)
!         nonlinbky0(j,i)=nonlinbky1(j,i)
!     end do
! end do

nonlinukx1 = nonlinukx1 + fukx
nonlinuky1 = nonlinuky1 + fuky
nonlinbkx1 = nonlinbkx1 + fbkx
nonlinbky1 = nonlinbky1 + fbky

rhok1 = rhok1 + deltaT*(1.5*nonlinrhok1 - 0.5*nonlinrhok0)
ukx1  = ukx1  + deltaT*(1.5*nonlinukx1  - 0.5*nonlinukx0)
uky1  = uky1  + deltaT*(1.5*nonlinuky1  - 0.5*nonlinuky0)
bkx1  = bkx1  + deltaT*(1.5*nonlinbkx1  - 0.5*nonlinbkx0)
bky1  = bky1  + deltaT*(1.5*nonlinbky1  - 0.5*nonlinbky0)

! Rename variables for saving and use them as initial values for next loop (AB2)
nonlinrhok0=nonlinrhok1
nonlinukx0=nonlinukx1
nonlinuky0=nonlinuky1
nonlinbkx0=nonlinbkx1
nonlinbky0=nonlinbky1



!**************** Write quantities
! Compute and write energy
if ( (mod(it,inrj) .eq. 0) ) then
    call save_energy(rhok1,ukx1,uky1,bkx1,bky1)
    time = time + dfloat(inrj)*deltaT
    open(52, file='out_time', position='append',form='formatted')
    write(52,*) time
close(52)

!****************Adaptive timestep
! open(53, file = 'out_nu', position='append', form='formatted')
! write(53,*) nu
! close(53)
! call adaptiveT(ukx1,rhok1,ta)
! open(51, file = 'out_deltaT', position='append', form='formatted')
! write(51,*) ta
! close(51)
! deltaT = ta*0.2d0 ! condition CFL
! nu = 1./ta*(1.d0/N)**4 ! condition CFL
! eta = nu
end if

! Compute and write spectra
if ( (mod(it,ispec) .eq. 0) ) then
    print *, it
    call save_spectra(rhok1,ukx1,uky1,bkx1,bky1,istore_sp)
endif

! Write fields
! rho, wz, jz, divu, divb
if (mod(it,ifields) .eq. 0) then
    call save_fields(rhok1,ukx1,uky1,bkx1,bky1,istore_fields)
end if

! Save spatio-temporal spectra
if (sts .eq. 1) then
    if (mod(it,ists) .eq. 0) then
        timests = timests + dfloat(ists)*deltaT
        call WriteSpatioTemporalSpectrum(ukx1, uky1, bkx1, bky1, rhok1, timests)
    end if
end if

end do ! end of the temporal loop

!****************End of the time loop

t2 = omp_get_wtime()
write(*,*) "cpu time", t2-t1
call end_fftw ! Deallocate plans

deallocate(ukx1, uky1, bkx1, bky1, rhok1)
! deallocate(ukx2, uky2, bkx2, bky2, rhok2) ! RK2
! deallocate(ukx3, uky3, bkx3, bky3, rhok3) ! RK4
! deallocate(ukx4, uky4, bkx4, bky4, rhok4) ! RK4
deallocate(nonlinrhok0, nonlinukx0, nonlinuky0, nonlinbkx0, nonlinbky0)
deallocate(nonlinrhok1, nonlinukx1, nonlinuky1, nonlinbkx1, nonlinbky1)
! deallocate(nonlinrhok2, nonlinukx2, nonlinuky2, nonlinbkx2, nonlinbky2) ! RK4
! deallocate(nonlinrhok3, nonlinukx3, nonlinuky3, nonlinbkx3, nonlinbky3) ! RK4
deallocate(fukx, fuky, fbkx, fbky)

print *, 'OK'
end program CMHD






!*****************************************************************
SUBROUTINE RandomInit(ukxi,ukyi)
use omp_lib
use parameters
use fftw_mod
use spectral_mod
use outputs
! Initial random field spectra
implicit none
double precision phase, EU, kmn, kmx !, spectri, knc
double complex :: ukxi(Nh,N), ukyi(Nh,N)
integer i, j

kmn = dk**2
kmx = (kinj*dk)**2
ukxi=0.
ukyi=0.
!$omp parallel do private(i,j,phase)
do i = 1, N
    if ((kd(1,i).le.kmx).and.(kd(1,i).ge.kmn)) then
        call random_number(phase)
        phase = 2*pi*phase
        ukxi(1,i) = (cos(phase) + imag*sin(phase)) / sqrt(kd(1,i))
        call random_number(phase)
        phase = 2*pi*phase
        ukyi(1,i) = (cos(phase) + imag*sin(phase)) / sqrt(kd(1,i))
    endif
    do j = 2, Nh-1
        if ((kd(j,i).le.kmx).and.(kd(j,i).ge.kmn)) then
            call random_number(phase)
            phase = 2*pi*phase
            ukxi(j,i) = 2*(cos(phase) + imag*sin(phase)) / sqrt(kd(j,i))
            call random_number(phase)
            phase = 2*pi*phase
            ukyi(j,i) = 2*(cos(phase) + imag*sin(phase)) / sqrt(kd(j,i))
        endif
    end do
    if ((kd(Nh,i).le.kmx).and.(kd(Nh,i).ge.kmn)) then
        call random_number(phase)
        phase = 2*pi*phase
        ukxi(Nh,i) = (cos(phase) + imag*sin(phase)) / sqrt(kd(Nh,i))
        call random_number(phase)
        phase = 2*pi*phase
        ukyi(Nh,i) = (cos(phase) + imag*sin(phase)) / sqrt(kd(Nh,i))
    endif
end do

ukxi = kill*ukxi ! Dealiasing
ukyi = kill*ukyi ! Dealiasing

call energy(ukxi,ukyi,EU)
ukxi = amp*ukxi/sqrt(EU)
ukyi = amp*ukyi/sqrt(EU)

RETURN
END SUBROUTINE RandomInit

SUBROUTINE GaussianInit(Akx,Aky)
use omp_lib
use parameters
use spectral_mod
use outputs
implicit none
double complex, intent(inout) :: Akx(Nh,N), Aky(Nh,N)
double precision phase, E
integer i, j

! TODO: set seed properly
!$omp parallel do private(i,j,phase)
do i = 1, N
    call random_number(phase)
    phase = 2*pi*phase
    Akx(1,i) = (cos(phase) + imag*sin(phase))*exp(-.5*((sqrt(kd(1,i))/dk-kinj)/(width))**2)
    call random_number(phase)
    phase = 2*pi*phase
    Aky(1,i) = (cos(phase) + imag*sin(phase))*exp(-.5*((sqrt(kd(1,i))/dk-kinj)/(width))**2)
    do j = 2, Nh-1
        call random_number(phase)
        phase = 2*pi*phase
        Akx(j,i) = 2*(cos(phase) + imag*sin(phase))*exp(-.5*((sqrt(kd(j,i))/dk-kinj)/(width))**2)
        call random_number(phase)
        phase = 2*pi*phase
        Aky(j,i) = 2*(cos(phase) + imag*sin(phase))*exp(-.5*((sqrt(kd(j,i))/dk-kinj)/(width))**2)
    end do
    call random_number(phase)
    phase = 2*pi*phase
    Akx(Nh,i) = (cos(phase) + imag*sin(phase))*exp(-.5*((sqrt(kd(Nh,i))/dk-kinj)/(width))**2)
    call random_number(phase)
    phase = 2*pi*phase
    Aky(Nh,i) = (cos(phase) + imag*sin(phase))*exp(-.5*((sqrt(kd(Nh,i))/dk-kinj)/(width))**2)
end do

Akx = kill*Akx ! Dealiasing
Aky = kill*Aky ! Dealiasing

call energy(Akx,Aky,E)
Akx = amp*Akx/sqrt(E)
Aky = amp*Aky/sqrt(E)

RETURN
END SUBROUTINE GaussianInit

!*****************************************************************
! Random forcing spectrum
!*****************************************************************
! SUBROUTINE RandomF(field)
! use parameters
! use fftw_mod
! implicit none
! double complex, intent(inout) :: field(Nh,N)
! double precision spectri
! double precision theta, knc, alpha0
! integer ii, jj

! call srand(seed)
! alpha0 = 100.
! field=0.
! do ii = 1, Na+1
!     do jj = 1, Na+1
!         knc = (kinj - real(ii-1) - real(jj-1))**4
!         spectri = dexp(-knc*alpha0)
!         theta = rand()*2.*pi
!         field(jj,ii) = spectri*(cos(theta) + imag*sin(theta))
!     end do
! end do

! RETURN
! END SUBROUTINE RandomF

!*****************************************************************

SUBROUTINE RandomF(Akx,Aky)
use omp_lib
use parameters
use spectral_mod
use outputs
implicit none
double complex, intent(inout) :: Akx(Nh,N), Aky(Nh,N)
double precision phase, kmn, kmx, E
integer i, j

kmn = dk**2
kmx = (kinj*dk)**2
Akx=0.
Aky=0.

!$omp parallel do private(i,j,phase)
do i = 1, N
    if ((kd(1,i).le.kmx).and.(kd(1,i).ge.kmn)) then
        call random_number(phase)
        phase = 2*pi*phase
        Akx(1,i) = (cos(phase) + imag*sin(phase)) / sqrt(kd(1,i)) ! CHECK: kd/dk?
        call random_number(phase)
        phase = 2*pi*phase
        Aky(1,i) = (cos(phase) + imag*sin(phase)) / sqrt(kd(1,i))
    endif
    do j = 2, Nh-1
        if ((kd(j,i).le.kmx).and.(kd(j,i).ge.kmn)) then
            call random_number(phase)
            phase = 2*pi*phase
            Akx(j,i) = 2*(cos(phase) + imag*sin(phase)) / sqrt(kd(j,i))
            call random_number(phase)
            phase = 2*pi*phase
            Aky(j,i) = 2*(cos(phase) + imag*sin(phase)) / sqrt(kd(j,i))
        endif
    end do
    if ((kd(Nh,i).le.kmx).and.(kd(Nh,i).ge.kmn)) then
        call random_number(phase)
        phase = 2*pi*phase
        Akx(Nh,i) = (cos(phase) + imag*sin(phase)) / sqrt(kd(Nh,i))
        call random_number(phase)
        phase = 2*pi*phase
        Aky(Nh,i) = (cos(phase) + imag*sin(phase)) / sqrt(kd(Nh,i))
    endif
end do

Akx = kill*Akx ! Dealiasing
Aky = kill*Aky ! Dealiasing

call incompressible_projection(Akx,Aky)
! call compressible_projection(Akx,Aky)

call energy(Akx,Aky,E)
Akx = famp*Akx/sqrt(E)
Aky = famp*Aky/sqrt(E)

RETURN
END SUBROUTINE RandomF

!*****************************************************************

SUBROUTINE CompressibleRandomF(Akx,Aky,a)
use omp_lib
use parameters
use spectral_mod
use outputs
implicit none
double complex, intent(inout) :: Akx(Nh,N), Aky(Nh,N)
double complex :: f1, f2
double precision a, phase, kmn, kmx, E 
integer i, j

! Mixed incompressible and compressible forcing, controlled by the parameter a
! a : fraction of compressible forcing
! F = i k x z \psi_1 - i k \psi_2
! \psi_1 = (1-a) f_1
! \psi_2 = a f_2

kmn = dk**2
kmx = (kinj*dk)**2
Akx=0.
Aky=0.
f1=0.
f2=0.

!$omp parallel do private(i,j,phase)
do i = 1, N
    f1 = 0.
    f2 = 0.
    if ((kd(1,i).le.kmx).and.(kd(1,i).ge.kmn)) then
        call random_number(phase)
        phase = 2*pi*phase
        f1 = (cos(phase) + imag*sin(phase)) / sqrt(kd(1,i)) ! CHECK: kd/dk?
        call random_number(phase)
        phase = 2*pi*phase
        f2 = (cos(phase) + imag*sin(phase)) / sqrt(kd(1,i))
    endif
    Akx(1,i) = imag*ky(1)*(1.-a)*f1 - imag*kx(i)*a*f2
    Aky(1,i) = -imag*kx(i)*(1.-a)*f1 - imag*ky(1)*a*f2
    do j = 2, Nh-1
        f1 = 0.
        f2 = 0.
        if ((kd(j,i).le.kmx).and.(kd(j,i).ge.kmn)) then
            call random_number(phase)
            phase = 2*pi*phase
            f1 = 2*(cos(phase) + imag*sin(phase)) / sqrt(kd(j,i))
            call random_number(phase)
            phase = 2*pi*phase
            f2 = 2*(cos(phase) + imag*sin(phase)) / sqrt(kd(j,i))
        endif
        Akx(j,i) = imag*ky(j)*(1.-a)*f1 - imag*kx(i)*a*f2
        Aky(j,i) = -imag*kx(i)*(1.-a)*f1 - imag*ky(j)*a*f2
    end do
    f1 = 0.
    f2 = 0.
    if ((kd(Nh,i).le.kmx).and.(kd(Nh,i).ge.kmn)) then
        call random_number(phase)
        phase = 2*pi*phase
        f1 = (cos(phase) + imag*sin(phase)) / sqrt(kd(Nh,i))
        call random_number(phase)
        phase = 2*pi*phase
        f2 = (cos(phase) + imag*sin(phase)) / sqrt(kd(Nh,i))
    endif
    Akx(Nh,i) = imag*ky(Nh)*(1.-a)*f1 - imag*kx(i)*a*f2
    Aky(Nh,i) = -imag*kx(i)*(1.-a)*f1 - imag*ky(Nh)*a*f2
end do

call energy(Akx,Aky,E)
Akx = famp*Akx/sqrt(E)
Aky = famp*Aky/sqrt(E)

RETURN
END SUBROUTINE CompressibleRandomF

!*****************************************************************

SUBROUTINE GaussianF(Akx,Aky)
use omp_lib
use parameters
use spectral_mod
use outputs
implicit none
double complex, intent(inout) :: Akx(Nh,N), Aky(Nh,N)
double precision phase, E
integer i, j

! TODO: set seed properly
!$omp parallel do private(i,j,phase)
do i = 1, N
    call random_number(phase)
    phase = 2*pi*phase
    Akx(1,i) = (cos(phase) + imag*sin(phase))*exp(-.5*((sqrt(kd(1,i))/dk-kinj)/(width))**2)
    call random_number(phase)
    phase = 2*pi*phase
    Aky(1,i) = (cos(phase) + imag*sin(phase))*exp(-.5*((sqrt(kd(1,i))/dk-kinj)/(width))**2)
    do j = 2, Nh-1
        call random_number(phase)
        phase = 2*pi*phase
        Akx(j,i) = 2*(cos(phase) + imag*sin(phase))*exp(-.5*((sqrt(kd(j,i))/dk-kinj)/(width))**2)
        call random_number(phase)
        phase = 2*pi*phase
        Aky(j,i) = 2*(cos(phase) + imag*sin(phase))*exp(-.5*((sqrt(kd(j,i))/dk-kinj)/(width))**2)
    end do
    call random_number(phase)
    phase = 2*pi*phase
    Akx(Nh,i) = (cos(phase) + imag*sin(phase))*exp(-.5*((sqrt(kd(Nh,i))/dk-kinj)/(width))**2)
    call random_number(phase)
    phase = 2*pi*phase
    Aky(Nh,i) = (cos(phase) + imag*sin(phase))*exp(-.5*((sqrt(kd(Nh,i))/dk-kinj)/(width))**2)
end do

Akx = kill*Akx ! Dealiasing
Aky = kill*Aky ! Dealiasing

call energy(Akx,Aky,E)
Akx = famp*Akx/sqrt(E)
Aky = famp*Aky/sqrt(E)

RETURN
END SUBROUTINE GaussianF

!*****************************************************************

SUBROUTINE PoloidalRandomF(Akx,Aky)
use omp_lib
use parameters
! use fftw_mod
use spectral_mod
use outputs
implicit none
double complex, intent(inout) :: Akx(Nh,N), Aky(Nh,N)
double precision phase, kmn, kmx, E
integer i, j

kmn = dk**2
kmx = (kinj*dk)**2
Akx=0.
Aky=0.

!$omp parallel do private(i,j,phase)
do i = 1, N
    if ((kd(1,i).le.kmx).and.(kd(1,i).ge.kmn)) then
        call random_number(phase)
        phase = 2*pi*phase
        Akx(1,i) = (cos(phase) + imag*sin(phase)) / sqrt(kd(1,i)) ! CHECK: kd/dk?
    endif
    do j = 2, Nh-1
        if ((kd(j,i).le.kmx).and.(kd(j,i).ge.kmn)) then
            call random_number(phase)
            phase = 2*pi*phase
            Akx(j,i) = 2*(cos(phase) + imag*sin(phase)) / sqrt(kd(j,i))
        endif
    end do
    if ((kd(Nh,i).le.kmx).and.(kd(Nh,i).ge.kmn)) then
        call random_number(phase)
        phase = 2*pi*phase
        Akx(Nh,i) = (cos(phase) + imag*sin(phase)) / sqrt(kd(Nh,i))
    endif
end do

!$omp parallel do private(i,j)
do i=1,N
    do j=1,Nh
        Aky(j,i) = -Akx(j,i)*kx(i)*ky(j)
        Akx(j,i) =  Akx(j,i)*ky(j)*ky(j)
    enddo
enddo

Akx = kill*Akx ! Dealiasing
Aky = kill*Aky ! Dealiasing

call energy(Akx,Aky,E)
Akx = famp*Akx/sqrt(E)
Aky = famp*Aky/sqrt(E)

RETURN
END SUBROUTINE PoloidalRandomF

!*****************************************************************

SUBROUTINE RandomF_Ak(fukx,fuky,fbkx,fbky)
use omp_lib
use parameters
use spectral_mod
use outputs
implicit none
double complex, intent(inout) :: fukx(Nh,N), fuky(Nh,N), fbkx(Nh,N), fbky(Nh,N)
double precision phase, kmn, kmx, E
integer i, j

kmn = dk**2 ! k=1
kmx = (kinj*dk)**2
fukx=0.
fuky=0.
fbkx=0.
fbky=0.

!$omp parallel do private(i,j,phase)
do i = 1, N
    if ((kd(1,i).le.kmx).and.(kd(1,i).ge.kmn)) then
        call random_number(phase)
        phase = 2*pi*phase
        fukx(1,i) = 0.  ! Because ky=0
        call random_number(phase)
        phase = 2*pi*phase
        fuky(1,i) = -(cos(phase) + imag*sin(phase)) / sqrt(kd(1,i))
        call random_number(phase)
        phase = 2*pi*phase
        fbkx(1,i) = -ky(1)/sqrt(kd(1,i))*(cos(phase) + imag*sin(phase)) / sqrt(kd(1,i)) 
        call random_number(phase)
        phase = 2*pi*phase
        fbky(1,i) = kx(i)/sqrt(kd(1,i))*(cos(phase) + imag*sin(phase)) / sqrt(kd(1,i))
    endif
    do j = 2, Nh
        if ((kd(j,i).le.kmx).and.(kd(j,i).ge.kmn)) then
            call random_number(phase)
            phase = 2*pi*phase
            fukx(j,i) = 2*kx(i)/ky(j)*(cos(phase) + imag*sin(phase)) / sqrt(kd(j,i))
            call random_number(phase)
            phase = 2*pi*phase
            fuky(j,i) = -2*(cos(phase) + imag*sin(phase)) / sqrt(kd(j,i))
            call random_number(phase)
            phase = 2*pi*phase
            fbkx(j,i) = -2*ky(j)/sqrt(kd(j,i))*(cos(phase) + imag*sin(phase)) / sqrt(kd(j,i)) 
            call random_number(phase)
            phase = 2*pi*phase
            fbky(j,i) = 2*kx(i)/sqrt(kd(j,i))*(cos(phase) + imag*sin(phase)) / sqrt(kd(j,i))
        endif
    end do
end do

fukx(1,1) = 0.0
fuky(1,1) = 0.0
fbkx(1,1) = 0.0
fbky(1,1) = 0.0

fukx = kill*fukx ! Dealiasing
fuky = kill*fuky ! Dealiasing
fbkx = kill*fbkx ! Dealiasing
fbky = kill*fbky ! Dealiasing

call incompressible_projection(fukx,fuky)
call incompressible_projection(fbkx,fbky)

call energy(fukx,fuky,E)
fukx = famp*fukx/sqrt(E)
fuky = famp*fuky/sqrt(E)
call energy(fbkx,fbky,E)
fbkx = famp*fbkx/sqrt(E)
fbky = famp*fbky/sqrt(E)

RETURN
END SUBROUTINE RandomF_Ak