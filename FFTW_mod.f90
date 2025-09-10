module FFTW_mod
use, intrinsic :: iso_c_binding
use parameters
implicit none
include "fftw3.f"

public

integer (kind=8) :: plan_for, plan_back
! integer :: N, Nh

save
contains

SUBROUTINE init_fftw
! integer, intent(in) :: N, Nh
double precision :: aa(N,N)
double complex :: bb(Nh,N)

call dfftw_plan_dft_r2c_2d_(plan_for, N, N, aa, bb, FFTW_MEASURE)
call dfftw_plan_dft_c2r_2d_(plan_back, N, N, bb, aa, FFTW_MEASURE)

RETURN
END SUBROUTINE init_fftw

SUBROUTINE FFT_PS(Ain, Aout)
double precision :: Ain(N,N)
double complex :: Aout(Nh,N)

call dfftw_execute_dft_r2c(plan_for, Ain, Aout)

RETURN
END SUBROUTINE FFT_PS

SUBROUTINE FFT_SP(Ain, Aout)
double complex :: Ain(Nh,N)
double precision :: Aout(N,N), norm

norm = 1. / real(N*N)

call dfftw_execute_dft_c2r(plan_back, Ain, Aout)
Aout = Aout*norm

RETURN
END SUBROUTINE FFT_SP


SUBROUTINE end_fftw
call dfftw_destroy_plan(plan_back)
call dfftw_destroy_plan(plan_for)
END SUBROUTINE end_fftw

END MODULE FFTW_mod
