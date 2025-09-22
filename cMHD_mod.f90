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

contains

SUBROUTINE RHS(rhok,ukx,uky,bkx,bky,nonlinrhok,nonlinukx,nonlinuky,nonlinbkx,nonlinbky)
double complex :: rhok(Nh,N), ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double complex :: nonlinrhok(Nh,N), nonlinukx(Nh,N), nonlinuky(Nh,N) 
double complex :: nonlinbkx(Nh,N), nonlinbky(Nh,N)
double complex :: dissip_nu, dissip_eta
integer :: i,j
double precision :: kx3, ky3, k4

! call init_rhs
allocate (tmp1(N,N), tmp2(N,N), tmp3(N,N))
allocate (tmpk1(Nh,N), tmpk2(Nh,N), tmpk3(Nh,N))

call RHS1(rhok,ukx,uky,nonlinrhok)
call RHS2(ukx,uky,bkx,bky,nonlinukx)
call RHS3(rhok,ukx,uky,bkx,bky,nonlinuky)
call RHS4(ukx,uky,bkx,bky,nonlinbkx)
call RHS5(ukx,uky,bkx,bky,nonlinbky)

! dissipation + dealiasing
! Hyperviscosity (bilaplacian) + hypoviscosity + dispersion
! TODO: Dissipation in rho? Do we need dissipation at large scale?
!$omp parallel do private(i,j,kx3,ky3,k4,dissip_nu,dissip_eta) 
do i = 1, N
    kx3 = kx(i)**3
    do j = 1, Nh
        k4 = kd(j,i)*kd(j,i)
        ky3 = ky(j)**3
        dissip_nu = (k4*nu  + alpha + imag*disp*(kx3+ky3))
        dissip_eta= (k4*eta + alpha + imag*disp*(kx3+ky3))
        nonlinrhok(j,i) = kill(j,i)*(nonlinrhok(j,i) - dissip_nu*rhok(j,i))
        nonlinukx (j,i) = kill(j,i)*(nonlinukx (j,i) - dissip_nu*ukx(j,i))
        nonlinuky (j,i) = kill(j,i)*(nonlinuky (j,i) - dissip_nu*uky(j,i))
        nonlinbkx (j,i) = kill(j,i)*(nonlinbkx (j,i) - dissip_eta*bkx(j,i))
        nonlinbky (j,i) = kill(j,i)*(nonlinbky (j,i) - dissip_eta*bky(j,i))
    end do
end do

! nonlinrhok = nonlinrhok - (kd*kd*nu + alpha)*rhok
! nonlinukx = nonlinukx - (kd*kd*nu + alpha)*ukx
! nonlinuky = nonlinuky - (kd*kd*nu + alpha)*uky
! nonlinbkx = nonlinbkx - (kd*kd*eta + alpha)*bkx
! nonlinbky = nonlinbky - (kd*kd*eta + alpha)*bky

! Dealiasing
! nonlinrhok = kill * nonlinrhok
! nonlinukx  = kill * nonlinukx
! nonlinuky  = kill * nonlinuky
! nonlinbkx  = kill * nonlinbkx
! nonlinbky  = kill * nonlinbky

deallocate(tmpk1, tmpk2, tmpk3, tmp1, tmp2, tmp3)
! call end_rhs

END SUBROUTINE RHS


SUBROUTINE RHS1(rhok,ukx,uky,nonlinrhok)
double complex, intent(in) :: rhok(Nh,N), ukx(Nh,N), uky(Nh,N)
double complex, intent(out) :: nonlinrhok(Nh,N)
! rhot = -divu - rho*(divu) - ux*rhox - uy*rhoy

! call init_rhs

! Add linear terms
call divergence(ukx,uky,tmpk2)
nonlinrhok = -tmpk2

tmpk1 = rhok
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! -rho*(uxdx + uydy)
tmp3 = -tmp1*tmp2

tmpk1 = ukx
call derivex(rhok,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! -ux*rhodx
tmp3 = tmp3 - tmp1*tmp2 

tmpk1 = uky
call derivey(rhok,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! -uy*rhody
tmp3 = tmp3 - tmp1*tmp2 

call FFT_PS(tmp3,tmpk3)
nonlinrhok = nonlinrhok + tmpk3

! call end_rhs

END SUBROUTINE RHS1

SUBROUTINE RHS2(ukx,uky,bkx,bky,nonlinukx)
double complex :: ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double complex :: nonlinukx(Nh,N)
! uxt = -ux*uxdx - uy*uxdy - by*(curl(b))

! call init_rhs

tmpk1 = ukx
call derivex(ukx,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! -ux*uxdx
tmp3 = -tmp1*tmp2 

tmpk1 = uky 
call derivey(ukx,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! -uy*uxdy
tmp3 = tmp3 - tmp1*tmp2 

tmpk1 = bky
call curl(bkx,bky,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! -by*(bydx - bxdy)
tmp3 = tmp3 - tmp1*tmp2 

call FFT_PS(tmp3,tmpk3)
nonlinukx = tmpk3

! No linear terms

! call end_rhs

END SUBROUTINE RHS2

SUBROUTINE RHS3(rhok,ukx,uky,bkx,bky,nonlinuky)
double complex :: rhok(Nh,N), ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double complex :: nonlinuky(Nh,N)
! uyt = curl(b) - ux*uydx - uy*uydy + (bx-rho)*(curl(b))

! call init_rhs

! Add linear terms
call curl(bkx,bky,tmpk2)
nonlinuky = tmpk2

! Nonlinear terms
tmpk1 = bkx - rhok
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! (bx-rho)*(bxdy - bydx)
tmp3 = tmp1*tmp2

tmpk1 = ukx
call derivex(uky,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! -ux*uydx
tmp3 = tmp3 - tmp1*tmp2 

tmpk1 = uky
call derivey(uky,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! -uy*uydy
tmp3 = tmp3 - tmp1*tmp2 

call FFT_PS(tmp3,tmpk3)

nonlinuky = nonlinuky + tmpk3

! call end_rhs

END SUBROUTINE RHS3

SUBROUTINE RHS4(ukx,uky,bkx,bky,nonlinbkx)
double complex :: ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double complex :: nonlinbkx(Nh,N)
! bxt = -dyuy + by*uxdy - ux*bxdx - uy*bxdy - bx*uydy

! call init_rhs

! Add linear terms
call derivey(uky,tmpk1)
nonlinbkx = -tmpk1

tmpk1 = bky
call derivey(ukx,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! by*uxdy
tmp3 = tmp1*tmp2 

tmpk1 = ukx
call derivex(bkx,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! -ux*bxdx
tmp3 = tmp3 - tmp1*tmp2 

tmpk1 = uky
call derivey(bkx,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! -uy*bxdy
tmp3 = tmp3 - tmp1*tmp2 

tmpk1 = bkx
call derivey(uky,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! -bx*uydy
tmp3 = tmp3 - tmp1*tmp2 

call FFT_PS(tmp3,tmpk3)
nonlinbkx = nonlinbkx + tmpk3

! call end_rhs

END SUBROUTINE RHS4

SUBROUTINE RHS5(ukx,uky,bkx,bky,nonlinbky)
double complex :: ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double complex :: nonlinbky(Nh,N)
! byt = uydx + bx*uydx - ux*bydx - uy*bydy - by*uxdx

! call init_rhs

! Add linear terms 
call derivex(uky,tmpk2)
nonlinbky = tmpk2

tmpk1 = bkx
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! bx*uydx
tmp3 = tmp1*tmp2 

tmpk1 = ukx
call derivex(bky,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! -ux*bydx
tmp3 = tmp3 - tmp1*tmp2

tmpk1 = uky
call derivey(bky,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! -uy*bydy
tmp3 = tmp3 - tmp1*tmp2 

tmpk1 = bky
call derivex(ukx,tmpk2)
call FFT_SP(tmpk1,tmp1)
call FFT_SP(tmpk2,tmp2)
! -by*uxdx
tmp3 = tmp3 - tmp1*tmp2

call FFT_PS(tmp3,tmpk3)
nonlinbky = nonlinbky + tmpk3

! call end_rhs

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

! SUBROUTINE init_rhs
! allocate(buffer(2*Nh, N)) ! factor 2 since complex = 2 reals

! ! Associate views
! call c_f_pointer(c_loc(buffer), tmp, [N,N])
! call c_f_pointer(c_loc(buffer), tmpk, [Nh,N])

! END SUBROUTINE init_rhs

! SUBROUTINE end_rhs
! deallocate(buffer)
! END SUBROUTINE end_rhs

END MODULE cMHD_mod