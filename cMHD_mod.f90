MODULE cMHD_mod
use, intrinsic :: iso_c_binding
use omp_lib
use parameters
use fftw_mod
use adaptive_mod
use spectral_mod
implicit none

! real(c_double), allocatable, target :: buffer1(:,:), buffer2(:,:), buffer3(:,:)
! real(c_double), pointer :: tmp1(:,:), tmp2(:,:), tmp3(:,:)               ! (N,N)
! complex(c_double_complex), pointer :: tmpk1(:,:), tmpk2(:,:), tmpk3(:,:) ! (Nh,N)

! TODO: pointers
double complex  , dimension(:,:), allocatable :: tmpk1, tmpk2, tmpk3
double precision, dimension(:,:), allocatable :: tmp1, tmp2, tmp3
double precision :: eps = 0.d-15

contains

SUBROUTINE RHS(rhok,ukx,uky,bkx,bky,nonlinrhok,nonlinukx,nonlinuky,nonlinbkx,nonlinbky)
double complex, dimension(Nh,N) :: rhok, ukx, uky, bkx, bky
double complex, dimension(Nh,N) :: nonlinrhok, nonlinukx, nonlinuky, nonlinbkx, nonlinbky
double complex :: dissip_nu, dissip_eta
double precision :: kx3, ky3, k4, nu3
integer :: i,j
! double precision, dimension(Nh,N) :: k4, k3

! call init_rhs
allocate (tmp1(N,N), tmp2(N,N), tmp3(N,N))
allocate (tmpk1(Nh,N), tmpk2(Nh,N), tmpk3(Nh,N))

call RHS1(rhok,ukx,uky,nonlinrhok)
call RHS2(rhok,ukx,uky,bkx,bky,nonlinukx)
call RHS3(rhok,ukx,uky,bkx,bky,nonlinuky)
call RHS4(ukx,uky,bkx,bky,nonlinbkx)
call RHS5(ukx,uky,bkx,bky,nonlinbky)

nu3 = nu/3.d0
! dissipation + dealiasing
! Hyperviscosity (bilaplacian) + hypoviscosity + dispersion
!$omp parallel do private(i,j,kx3,ky3,k4,dissip_nu,dissip_eta) 
do i = 1, N
    kx3 = kx(i)**3
    do j = 1, Nh
        k4 = kd(j,i)*kd(j,i)
        ky3 = ky(j)**3
        dissip_nu = (k4*nu  + alpha + imag*disp*(kx3+ky3)) ! Hyperviscosity
        dissip_eta= (k4*eta + alpha - imag*disp*(kx3+ky3))
        ! dissip_nu = (kd(j,i)*nu  + alpha + imag*disp*(kx3+ky3)) ! Viscosity
        ! dissip_eta= (kd(j,i)*eta + alpha - imag*disp*(kx3+ky3))

        nonlinrhok(j,i) = kill(j,i)*(nonlinrhok(j,i)) ! - dissip_nu*rhok(j,i))
        nonlinukx (j,i) = kill(j,i)*(nonlinukx (j,i) - dissip_nu*ukx(j,i) - nu3 * kd(j,i)*kx(i)*(kx(i)*ukx(j,i) + ky(j)*uky(j,i)))
        nonlinuky (j,i) = kill(j,i)*(nonlinuky (j,i) - dissip_nu*uky(j,i) - nu3 * kd(j,i)*ky(j)*(kx(i)*ukx(j,i) + ky(j)*uky(j,i)))
        nonlinbkx (j,i) = kill(j,i)*(nonlinbkx (j,i) - dissip_eta*bkx(j,i))
        nonlinbky (j,i) = kill(j,i)*(nonlinbky (j,i) - dissip_eta*bky(j,i))
        
        ! No dealiasing
        ! nonlinrhok(j,i) = (nonlinrhok(j,i)) ! - dissip_nu*rhok(j,i))
        ! nonlinukx (j,i) = (nonlinukx (j,i) - dissip_nu*ukx(j,i) - nu/3. * kd(j,i)*kx(i)*(kx(i)*ukx(j,i) + ky(j)*uky(j,i)))
        ! nonlinuky (j,i) = (nonlinuky (j,i) - dissip_nu*uky(j,i) - nu/3. * kd(j,i)*ky(j)*(kx(i)*ukx(j,i) + ky(j)*uky(j,i)))
        ! nonlinbkx (j,i) = (nonlinbkx (j,i) - dissip_eta*bkx(j,i))
        ! nonlinbky (j,i) = (nonlinbky (j,i) - dissip_eta*bky(j,i))

    end do
end do

! call divergence(ukx,uky,tmpk1)
! call derivex(tmpk1,tmpk2)
! call derivey(tmpk1,tmpk3)

! k4 = kd*kd
! k3 = kd**(1.5)
! nonlinrhok = kill*nonlinrhok
! nonlinukx  = kill*(nonlinukx - (k4*nu + alpha + imag*disp*(k3))*ukx + nu/3. * kd*tmpk2)
! nonlinuky  = kill*(nonlinuky - (k4*nu + alpha + imag*disp*(k3))*uky + nu/3. * kd*tmpk3)
! nonlinbkx  = kill*(nonlinbkx - (k4*eta+ alpha + imag*disp*(k3))*bkx)
! nonlinbky  = kill*(nonlinbky - (k4*eta+ alpha + imag*disp*(k3))*bky)

deallocate(tmpk1, tmpk2, tmpk3, tmp1, tmp2, tmp3)
! call end_rhs

END SUBROUTINE RHS


SUBROUTINE RHS1(rhok,ukx,uky,nonlinrhok)
double complex, intent(in) :: rhok(Nh,N), ukx(Nh,N), uky(Nh,N)
double complex, intent(out) :: nonlinrhok(Nh,N)
! integer i,j
! rhot = -divu - rho*(divu) - ux*rhox - uy*rhoy

! Add linear terms
call divergence(ukx,uky,tmpk2)
nonlinrhok = -tmpk2

tmpk1 = rhok
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! -rho*(uxdx + uydy)
! !$omp parallel do private(i,j)
! do i=1,N
!     do j=1,N
!         tmp3(j,i) = -tmp2(j,i)*tmp1(j,i) 
!     end do
! end do
tmp3 = -tmp1*tmp2

tmpk1 = ukx
call derivex(rhok,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! -ux*rhodx
! !$omp parallel do private(i,j)
! do i=1,N
!     do j=1,N
!         tmp3(j,i) = tmp3(j,i) - tmp2(j,i)*tmp1(j,i) 
!     end do
! end do
tmp3 = tmp3 - tmp1*tmp2 

tmpk1 = uky
call derivey(rhok,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! -uy*rhody
! !$omp parallel do private(i,j)
! do i=1,N
!     do j=1,N
!         tmp3(j,i) = tmp3(j,i) - tmp2(j,i)*tmp1(j,i) 
!     end do
! end do
tmp3 = tmp3 - tmp1*tmp2 

call FFT_PS(tmp3,tmpk3)
nonlinrhok = nonlinrhok + tmpk3

END SUBROUTINE RHS1

SUBROUTINE RHS2(rhok,ukx,uky,bkx,bky,nonlinukx)
double complex :: rhok(Nh,N), ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double complex :: nonlinukx(Nh,N)
! integer :: i,j
! uxt = -ux*uxdx - uy*uxdy - by*(curl(b))/(1+rho)

! No linear terms

tmpk1 = rhok
tmpk2 = bky
call curl(bkx,bky,tmpk3)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
call FFT_SP(tmpk3,tmp3)
! -by*(bydx - bxdy)/(1+rho)
! !$omp parallel do private(i,j)
! do i=1,N
!     do j=1,N
!         tmp3(j,i) = -tmp2(j,i)*tmp3(j,i) / (1.d0 + tmp1(j,i) + eps)
!     end do
! end do
tmp3 = -tmp2*tmp3/(1.d0+tmp1+eps)

!Pressure term
tmp1 = 1.d0 + tmp1 ! 1 + rho = rho_tot ; rho_0 = 1
call derivex(rhok,tmpk2)
call FFT_SP(tmpk2,tmp2)
! !$omp parallel do private(i,j)
! do i=1,N
!     do j=1,N
!         tmp3(j,i) = tmp3(j,i) - cspeed**2 * tmp2(j,i)*tmp1(j,i)**(gamma-1.d0) / (tmp1(j,i)+eps)
!     end do
! end do
tmp3 = tmp3 - cspeed**2 * tmp2*tmp1**(gamma-1.d0) / (tmp1+eps)

tmpk1 = ukx
call derivex(ukx,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! -ux*uxdx
! !$omp parallel do private(i,j)
! do i=1,N
!     do j=1,N
!         tmp3(j,i) = tmp3(j,i) - tmp2(j,i)*tmp1(j,i) 
!     end do
! end do
tmp3 = tmp3 - tmp1*tmp2 

tmpk1 = uky 
call derivey(ukx,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! -uy*uxdy
! !$omp parallel do private(i,j)
! do i=1,N
!     do j=1,N
!         tmp3(j,i) = tmp3(j,i) - tmp2(j,i)*tmp1(j,i) 
!     end do
! end do
tmp3 = tmp3 - tmp1*tmp2 

call FFT_PS(tmp3,tmpk3)
nonlinukx = tmpk3

END SUBROUTINE RHS2

SUBROUTINE RHS3(rhok,ukx,uky,bkx,bky,nonlinuky)
double complex :: rhok(Nh,N), ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double complex :: nonlinuky(Nh,N)
! integer :: i,j
! uyt = - ux*uydx - uy*uydy + (bx+b0)*curl(b)/(1+rho) with b0=1

! Add linear terms
call curl(bkx,bky,tmpk2)
nonlinuky = tmpk2
! nonlinuky = 0.d0

! Nonlinear terms
tmpk1 = rhok
tmpk3 = bkx
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
call FFT_SP(tmpk3,tmp3)
! (bx+1)*(bxdy - bydx) / (1+rho)
! Correct expression: Not working
! tmp3 = (1.d0 + tmp3)*tmp2/(1.d0+tmp1+eps)
! Taylor expansion: Working
! !$omp parallel do private(i,j)
! do i=1,N
!     do j=1,N
!         tmp3(j,i) = tmp3(j,i)*tmp2(j,i)/( 1.d0 + tmp1(j,i) + eps)
!         tmp3(j,i) = tmp3(j,i) - tmp2(j,i)*tmp1(j,i) 
!         tmp3(j,i) = tmp3(j,i) + tmp2(j,i)*tmp1(j,i)**2
!     end do
! end do
tmp3 = tmp3*tmp2/(1.d0+tmp1+eps)
tmp3 = tmp3 - tmp1*tmp2
tmp3 = tmp3 + 0.5*tmp1**2*tmp2

! NOTE: The correct expression is the one not working, wihtout adding the linear term in nonlinuky, but that doesn't work.
! So I add the linear term and do a Taylor expansion only of that term. TODO 

!Pressure term
tmp1 = 1.d0 + tmp1 ! 1 + rho = rho_tot ; rho_0 = 1
call derivey(rhok,tmpk2)
call FFT_SP(tmpk2,tmp2)
! !$omp parallel do private(i,j)
! do i=1,N
!     do j=1,N
!         tmp3(j,i) = tmp3(j,i) - cspeed**2 * tmp2(j,i)*tmp1(j,i)**(gamma-1.d0) / (tmp1(j,i)+eps)
!     end do
! end do
tmp3 = tmp3 - cspeed**2 * tmp2*tmp1**(gamma-1.d0) / (tmp1+eps)

tmpk1 = ukx
call derivex(uky,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! -ux*uydx
! !$omp parallel do private(i,j)
! do i=1,N
!     do j=1,N
!         tmp3(j,i) = tmp3(j,i) - tmp2(j,i)*tmp1(j,i) 
!     end do
! end do
tmp3 = tmp3 - tmp1*tmp2 

tmpk1 = uky
call derivey(uky,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! -uy*uydy
! !$omp parallel do private(i,j)
! do i=1,N
!     do j=1,N
!         tmp3(j,i) = tmp3(j,i) - tmp2(j,i)*tmp1(j,i) 
!     end do
! end do
tmp3 = tmp3 - tmp1*tmp2 

call FFT_PS(tmp3,tmpk3)

nonlinuky = nonlinuky + tmpk3 

END SUBROUTINE RHS3

SUBROUTINE RHS4(ukx,uky,bkx,bky,nonlinbkx)
double complex :: ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double complex :: nonlinbkx(Nh,N)
! bxt = -dyuy - bx*uydy + by*uxdy - ux*bxdx - uy*bxdy 

! Add linear terms
call derivey(uky,tmpk2)
nonlinbkx = -tmpk2

tmpk1 = ukx
tmpk2 = bky
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! ux*by
tmp3 = tmp1*tmp2

tmpk1 = uky
tmpk2 = bkx
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! -uy*bx
tmp3 = tmp3 - tmp1*tmp2 

call FFT_PS(tmp3,tmpk3)
call derivey(tmpk3,tmpk2)
nonlinbkx = nonlinbkx + tmpk2

END SUBROUTINE RHS4

SUBROUTINE RHS5(ukx,uky,bkx,bky,nonlinbky)
double complex :: ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double complex :: nonlinbky(Nh,N)
! byt = uydx + bx*uydx - ux*bydx - uy*bydy - by*uxdx

! Add linear terms 
call derivex(uky,tmpk2)
nonlinbky = tmpk2

tmpk1 = uky
tmpk2 = bkx
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! bx*uy
tmp3 = tmp1*tmp2

tmpk1 = ukx
tmpk2 = bky
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! -by*ux
tmp3 = tmp3 - tmp1*tmp2 

call FFT_PS(tmp3,tmpk3)
call derivex(tmpk3,tmpk2)
nonlinbky = nonlinbky + tmpk2

END SUBROUTINE RHS5

SUBROUTINE dissipation(ukx,uky,bkx,bky)
double complex, intent(inout) :: ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double precision :: k4, kx3, ky3
double complex :: dissip_nu, dissip_eta
integer :: i,j

! Implicit method for dissipation term
! Hypoviscosity (bilaplacian) + dispersion
! TODO: Dissipation in rho? Do we need dissipation at large scale?
!$omp parallel do private(i,j,kx3,ky3,k4,dissip_nu,dissip_eta) 
do i = 1, N
    kx3 = kx(i)**3
    do j = 1, Nh
        k4 = kd(j,i)*kd(j,i)
        ky3 = ky(j)**3
        dissip_nu = exp(-(k4*nu  + alpha + imag*disp*(kx3+ky3))*deltaT)
        dissip_eta= exp(-(k4*eta + alpha + imag*disp*(kx3+ky3))*deltaT)
        ukx(j,i) = ukx(j,i)*dissip_nu
        uky(j,i) = uky(j,i)*dissip_nu
        bkx(j,i) = bkx(j,i)*dissip_eta
        bky(j,i) = bky(j,i)*dissip_eta
    end do
end do

end SUBROUTINE dissipation

SUBROUTINE check_nan(arr)
double complex, intent(in) :: arr(Nh, N)   ! double precision complex
double precision :: aa

aa = sum(abs(arr))
if (isnan(aa)) then
    error stop "ERROR: NaN detected"
endif

END SUBROUTINE check_nan

! SUBROUTINE init_rhs
! allocate(buffer1(2*Nh, N), buffer2(2*Nh,N), buffer3(2*Nh,N)) ! factor 2 since complex = 2 reals

! ! Associate views
! call c_f_pointer(c_loc(buffer1), tmp1, [N,N])
! call c_f_pointer(c_loc(buffer1), tmpk1, [Nh,N])
! call c_f_pointer(c_loc(buffer2), tmp2, [N,N])
! call c_f_pointer(c_loc(buffer2), tmpk2, [Nh,N])
! call c_f_pointer(c_loc(buffer3), tmp3, [N,N])
! call c_f_pointer(c_loc(buffer3), tmpk3, [Nh,N])

! END SUBROUTINE init_rhs

! SUBROUTINE end_rhs
! deallocate(buffer1, buffer2, buffer3)
! END SUBROUTINE end_rhs

END MODULE cMHD_mod