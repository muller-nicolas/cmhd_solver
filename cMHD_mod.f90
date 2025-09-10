MODULE cMHD_mod

use, intrinsic :: iso_c_binding
use parameters
use fftw_mod
use spectral_mod
implicit none

double complex :: rhokdx(Nh,N), rhokdy(Nh,N), ukxdx(Nh,N), ukxdy(Nh,N), ukydx(Nh,N), ukydy(Nh,N)
double complex :: bkxdx(Nh,N), bkxdy(Nh,N), bkydx(Nh,N), bkydy(Nh,N)
double complex :: tmpk1(Nh,N), tmpk2(Nh,N), tmpk3(Nh,N), tmpk4(Nh,N)
double complex :: field1(Nh,N), field2(Nh,N), field3(Nh,N), field4(Nh,N), field5(Nh,N)
double precision :: tmp1(N,N), tmp2(N,N), tmp3(N,N), tmp4(N,N)

contains

SUBROUTINE RHS(rhok,ukx,uky,bkx,bky,nonlinrhok,nonlinukx,nonlinuky,nonlinbkx,nonlinbky)
double complex :: rhok(Nh,N), ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double complex :: nonlinrhok(Nh,N), nonlinukx(Nh,N), nonlinuky(Nh,N) 
double complex :: nonlinbkx(Nh,N), nonlinbky(Nh,N)

call RHS1(rhok,ukx,uky,nonlinrhok)
call RHS2(ukx,uky,bkx,bky,nonlinukx)
call RHS3(rhok,ukx,uky,bkx,bky,nonlinuky)
call RHS4(ukx,uky,bkx,bky,nonlinbkx)
call RHS5(ukx,uky,bkx,bky,nonlinbky)

! print *, "rho: ", sum(abs(nonlinrhok))

! Dealiasing
nonlinrhok = kill * nonlinrhok
nonlinukx  = kill * nonlinukx
nonlinuky  = kill * nonlinuky
nonlinbkx  = kill * nonlinbkx
nonlinbky  = kill * nonlinbky

END SUBROUTINE RHS

subroutine RHS1(rhok,ukx,uky,nonlinrhok)
double complex :: rhok(Nh,N), ukx(Nh,N), uky(Nh,N)
double complex :: nonlinrhok(Nh,N), divu(Nh,N)

field1 = rhok
field2 = ukx
field3 = uky

! Add linear terms
call derivex(ukx,ukxdx)
call derivey(uky,ukydy)

divu = ukxdx + ukydy

nonlinrhok = -divu

call FFT_SP(field1,tmp1)
call FFT_SP(divu,tmp2)
! rho*(uxdx + uydy)
tmp4 = tmp1*(tmp2) 
call FFT_PS(tmp4,tmpk1)

call derivex(rhok,rhokdx)
call FFT_SP(field2,tmp1)
call FFT_SP(rhokdx,tmp2)
! ux*rhodx
tmp3 = tmp1*tmp2 
call FFT_PS(tmp3,tmpk2)

call derivey(rhok,rhokdy)
call FFT_SP(field3,tmp1)
call FFT_SP(rhokdy,tmp2)
! uy*rhody
tmp3 = tmp1*tmp2 
call FFT_PS(tmp3,tmpk3)

nonlinrhok = nonlinrhok - (tmpk1 + tmpk2 + tmpk3)

end subroutine RHS1

SUBROUTINE RHS2(ukx,uky,bkx,bky,nonlinukx)
double complex :: ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double complex :: nonlinukx(Nh,N)

field2 = ukx
field3 = uky 
field5 = bky

! Add linear terms
nonlinukx = 0.

call derivex(ukx,ukxdx)
call FFT_SP(field2,tmp1)
call FFT_SP(ukxdx,tmp2)
! ux*uxdx
tmp3 = tmp1*tmp2 
call FFT_PS(tmp3,tmpk1)

call derivey(ukx,ukxdy)
call FFT_SP(field3,tmp1)
call FFT_SP(ukxdy,tmp2)
! uy*uxdy
tmp3 = tmp1*tmp2 
call FFT_PS(tmp3,tmpk2)

call derivex(bky,bkydx)
call derivey(bkx,bkxdy)
call FFT_SP(field5,tmp1)
call FFT_SP(bkydx,tmp2)
call FFT_SP(bkxdy,tmp3)
! by*(bydx - bxdy)
tmp4 = tmp1*(tmp2 - tmp3) 
call FFT_PS(tmp4,tmpk3)

nonlinukx = nonlinukx - (tmpk1 + tmpk2 + tmpk3)

END SUBROUTINE RHS2

SUBROUTINE RHS3(rhok,ukx,uky,bkx,bky,nonlinuky)
double complex :: rhok(Nh,N), ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double complex :: nonlinuky(Nh,N)

field1 = rhok
field2 = ukx
field3 = uky
field4 = bkx

! Add linear terms
call derivex(bky,bkydx)
call derivey(bkx,bkxdy)
nonlinuky = bkydx - bkxdy

! Nonlinear terms
call derivex(uky,ukydx)
call FFT_SP(field2,tmp1)
call FFT_SP(ukydx,tmp2)
! ux*uydx
tmp3 = tmp1*tmp2 
call FFT_PS(tmp3,tmpk1)

call derivey(uky,ukydy)
call FFT_SP(field3,tmp1)
call FFT_SP(ukydy,tmp2)
! uy*uydy
tmp3 = tmp1*tmp2 
call FFT_PS(tmp3,tmpk2)

call derivex(bky,bkydx)
call derivey(bkx,bkxdy)
call FFT_SP(field4,tmp1)
call FFT_SP(bkxdy,tmp2)
call FFT_SP(bkydx,tmp3)
! bx*(bxdy - bydx)
tmp4 = tmp1*(tmp2 - tmp3) 
call FFT_PS(tmp4,tmpk3)

call FFT_SP(field1,tmp1)
! rho*(bxdy - bydx)
tmp4 = -tmp1*(tmp2-tmp3) 
call FFT_PS(tmp4,tmpk4)

nonlinuky = nonlinuky - (tmpk1 + tmpk2 + tmpk3 + tmpk4)

END SUBROUTINE RHS3

SUBROUTINE RHS4(ukx,uky,bkx,bky,nonlinbkx)
double complex :: ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double complex :: nonlinbkx(Nh,N)

field2 = ukx
field3 = uky
field4 = bkx
field5 = bky

! Add linear terms
call derivey(uky,ukydy)
nonlinbkx = -ukydy

call derivey(ukx,ukxdy)
call FFT_SP(field5,tmp1)
call FFT_SP(ukxdy,tmp2)
! by*uxdy
tmp3 = tmp1*tmp2 
call FFT_PS(tmp3,tmpk1)

call derivex(bkx,bkxdx)
call FFT_SP(field2,tmp1)
call FFT_SP(bkxdx,tmp2)
! ux*bxdx
tmp3 = tmp1*tmp2 
call FFT_PS(tmp3,tmpk2)

call derivey(bkx,bkxdy)
call FFT_SP(field3,tmp1)
call FFT_SP(bkxdy,tmp2)
! uy*bxdy
tmp3 = tmp1*tmp2 
call FFT_PS(tmp3,tmpk3)

call derivey(uky,ukydy)
call FFT_SP(field4,tmp1)
call FFT_SP(ukydy,tmp2)
! bx*uydy
tmp3 = tmp1*tmp2 
call FFT_PS(tmp3,tmpk4)

nonlinbkx = nonlinbkx + (tmpk1 - tmpk2 - tmpk3 - tmpk4)

END SUBROUTINE RHS4

SUBROUTINE RHS5(ukx,uky,bkx,bky,nonlinbky)

double complex :: ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
double complex :: nonlinbky(Nh,N)

field2 = ukx
field3 = uky
field4 = bkx
field5 = bky

! Add linear terms 
call derivex(uky,ukydx)
nonlinbky = ukydx

! call derivex(uky,ukydx)
call FFT_SP(field4,tmp1)
call FFT_SP(ukydx,tmp2)
! bx*uydx
tmp3 = tmp1*tmp2 
call FFT_PS(tmp3,tmpk1)

call derivex(bky,bkydx)
call FFT_SP(field2,tmp1)
call FFT_SP(bkydx,tmp2)
! ux*bydx
tmp3 = tmp1*tmp2
call FFT_PS(tmp3,tmpk2)

call derivey(bky,bkydy)
call FFT_SP(field3,tmp1)
call FFT_SP(bkydy,tmp2)
! uy*bydy
tmp3 = tmp1*tmp2 
call FFT_PS(tmp3,tmpk3)

call derivex(ukx,ukxdx)
call FFT_SP(field5,tmp1)
call FFT_SP(ukxdx,tmp2)
! by*uxdx
tmp3 = tmp1*tmp2
call FFT_PS(tmp3,tmpk4)

nonlinbky = nonlinbky + (tmpk1 - tmpk2 - tmpk3 - tmpk4)

END SUBROUTINE RHS5

SUBROUTINE dissipation(rhok,ukx,uky,bkx,bky)
double complex :: rhok(Nh,N), ukx(Nh,N), uky(Nh,N), bkx(Nh,N), bky(Nh,N)
integer :: i,j

! Implicit method for dissipation term
! Hypoviscosity (bilaplacian) + dispersion
! TODO: Dissipation in rho?
do i = 1, N
    do j = 1, Nh
        rhok(j,i) = rhok(j,i)*exp(-(kd(j,i)*nu*kd(j,i)+imag*disp*(kx(i)**3+ky(j)**3))*deltaT)
        ukx(j,i) = ukx(j,i)*exp(-(kd(j,i)*nu*kd(j,i)+imag*disp*(kx(i)**3+ky(j)**3))*deltaT)
        uky(j,i) = uky(j,i)*exp(-(kd(j,i)*nu*kd(j,i)+imag*disp*(kx(i)**3+ky(j)**3))*deltaT)
        bkx(j,i) = bkx(j,i)*exp(-(kd(j,i)*eta*kd(j,i)+imag*disp*(kx(i)**3+ky(j)**3))*deltaT)
        bky(j,i) = bky(j,i)*exp(-(kd(j,i)*eta*kd(j,i)+imag*disp*(kx(i)**3+ky(j)**3))*deltaT)
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
