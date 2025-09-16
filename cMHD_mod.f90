MODULE cMHD_mod
use, intrinsic :: iso_c_binding
use omp_lib
use parameters
use fftw_mod
use adaptive_mod
use spectral_mod
implicit none

! double complex :: tmpk1(Nh,N) 
! double complex :: field(Nh,N) 
! double precision :: tmp1(N,N), tmp2(N,N) 

double complex  , allocatable :: tmpk1(:,:), tmpk2(:,:)
double precision, allocatable :: tmp1(:,:), tmp2(:,:) 

contains

SUBROUTINE RHS(rhok,ukx,uky,bkx,bky,nonlinrhok,nonlinukx,nonlinuky,nonlinbkx,nonlinbky)
double complex :: rhok(Nh,N), ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double complex :: nonlinrhok(Nh,N), nonlinukx(Nh,N), nonlinuky(Nh,N) 
double complex :: nonlinbkx(Nh,N), nonlinbky(Nh,N)

allocate(tmpk1(Nh,N), tmpk2(Nh,N), tmp1(N,N), tmp2(N,N))

call RHS1(rhok,ukx,uky,nonlinrhok)
call RHS2(ukx,uky,bkx,bky,nonlinukx)
call RHS3(rhok,ukx,uky,bkx,bky,nonlinuky)
call RHS4(ukx,uky,bkx,bky,nonlinbkx)
call RHS5(ukx,uky,bkx,bky,nonlinbky)

! Dealiasing
! nonlinrhok = kill * nonlinrhok
! nonlinukx  = kill * nonlinukx
! nonlinuky  = kill * nonlinuky
! nonlinbkx  = kill * nonlinbkx
! nonlinbky  = kill * nonlinbky

deallocate(tmpk1, tmpk2, tmp1, tmp2)

END SUBROUTINE RHS

subroutine RHS1(rhok,ukx,uky,nonlinrhok)
double complex, intent(in) :: rhok(Nh,N), ukx(Nh,N), uky(Nh,N)
double complex, intent(out) :: nonlinrhok(Nh,N)
! double complex :: divu(Nh,N)
! double complex, allocatable :: tmpk1(:,:), tmpk2(:,:)
! rhot = -divu - rho*(divu) - ux*rhox - uy*rhoy

! Add linear terms
call divergence(ukx,uky,tmpk2)

nonlinrhok = -tmpk2

tmpk1 = rhok
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! rho*(uxdx + uydy)
tmp1 = tmp1*tmp2
call FFT_PS(tmp1,tmpk1)
nonlinrhok = nonlinrhok - tmpk1

tmpk1 = ukx
call derivex(rhok,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! ux*rhodx
tmp1 = tmp1*tmp2 
call FFT_PS(tmp1,tmpk1)
nonlinrhok = nonlinrhok - tmpk1

tmpk1 = uky
call derivey(rhok,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! uy*rhody
tmp1 = tmp1*tmp2 
call FFT_PS(tmp1,tmpk1)
nonlinrhok = nonlinrhok - tmpk1

end subroutine RHS1

SUBROUTINE RHS2(ukx,uky,bkx,bky,nonlinukx)
double complex :: ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double complex :: nonlinukx(Nh,N)
! uxt = -ux*uxdx - uy*uxdy - by*(curl(b))

! Add linear terms
nonlinukx = 0.

tmpk1 = ukx
call derivex(ukx,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! ux*uxdx
tmp1 = tmp1*tmp2 
call FFT_PS(tmp1,tmpk1)
nonlinukx = nonlinukx - tmpk1

tmpk1 = uky 
call derivey(ukx,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! uy*uxdy
tmp1 = tmp1*tmp2 
call FFT_PS(tmp1,tmpk1)
nonlinukx = nonlinukx - tmpk1

tmpk1 = bky
call curl(bkx,bky,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! by*(bydx - bxdy)
tmp1 = tmp1*tmp2 
call FFT_PS(tmp1,tmpk1)
nonlinukx = nonlinukx - tmpk1

END SUBROUTINE RHS2

SUBROUTINE RHS3(rhok,ukx,uky,bkx,bky,nonlinuky)
double complex :: rhok(Nh,N), ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double complex :: nonlinuky(Nh,N)
! uyt = cur(b) - ux*uydx - uy*uydy + (bx-rho)*(curl(b))

! Add linear terms
call curl(bkx,bky,tmpk2)
nonlinuky = tmpk2

! Nonlinear terms
tmpk1 = bkx - rhok
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! (bx-rho)*(bxdy - bydx)
tmp1 = tmp1*tmp2
call FFT_PS(tmp1,tmpk1)
nonlinuky = nonlinuky + tmpk1

tmpk1 = ukx
call derivex(uky,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! ux*uydx
tmp1 = tmp1*tmp2 
call FFT_PS(tmp1,tmpk1)
nonlinuky = nonlinuky - tmpk1

tmpk1 = uky
call derivey(uky,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! uy*uydy
tmp1 = tmp1*tmp2 
call FFT_PS(tmp1,tmpk1)
nonlinuky = nonlinuky - tmpk1

END SUBROUTINE RHS3

SUBROUTINE RHS4(ukx,uky,bkx,bky,nonlinbkx)
double complex :: ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double complex :: nonlinbkx(Nh,N)
! bxt = -dyuy + by*uxdy - ux*bxdx - uy*bxdy - bx*uydy

! Add linear terms
call derivey(uky,tmpk1)
nonlinbkx = -tmpk1

tmpk1 = bky
call derivey(ukx,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! by*uxdy
tmp1 = tmp1*tmp2 
call FFT_PS(tmp1,tmpk1)
nonlinbkx = nonlinbkx + tmpk1

tmpk1 = ukx
call derivex(bkx,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! ux*bxdx
tmp1 = tmp1*tmp2 
call FFT_PS(tmp1,tmpk1)
nonlinbkx = nonlinbkx - tmpk1

tmpk1 = uky
call derivey(bkx,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! uy*bxdy
tmp1 = tmp1*tmp2 
call FFT_PS(tmp1,tmpk1)
nonlinbkx = nonlinbkx - tmpk1

tmpk1 = bkx
call derivey(uky,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! bx*uydy
tmp1 = tmp1*tmp2 
call FFT_PS(tmp1,tmpk1)
nonlinbkx = nonlinbkx - tmpk1

END SUBROUTINE RHS4

SUBROUTINE RHS5(ukx,uky,bkx,bky,nonlinbky)
double complex :: ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double complex :: nonlinbky(Nh,N)
! byt = uydx + bx*uydx - ux*bydx - uy*bydy - by*uxdx

! Add linear terms 
call derivex(uky,tmpk1)
nonlinbky = tmpk1

tmpk1 = bkx
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! bx*uydx
tmp1 = tmp1*tmp2 
call FFT_PS(tmp1,tmpk1)
nonlinbky = nonlinbky + tmpk1

tmpk1 = ukx
call derivex(bky,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! ux*bydx
tmp1 = tmp1*tmp2
call FFT_PS(tmp1,tmpk1)
nonlinbky = nonlinbky - tmpk1

tmpk1 = uky
call derivey(bky,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! uy*bydy
tmp1 = tmp1*tmp2 
call FFT_PS(tmp1,tmpk1)
nonlinbky = nonlinbky - tmpk1

tmpk1 = bky
call derivex(ukx,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! by*uxdx
tmp1 = tmp1*tmp2
call FFT_PS(tmp1,tmpk1)
nonlinbky = nonlinbky - tmpk1

END SUBROUTINE RHS5

SUBROUTINE dissipation(rhok,ukx,uky,bkx,bky)
double complex, intent(inout) :: rhok(Nh,N), ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
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
        rhok(j,i) = rhok(j,i)*dissip_nu
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
if (aa .ne. aa) then
    print *, "ERROR: NaN detected"
    stop 1
endif

END SUBROUTINE check_nan

END MODULE cMHD_mod
