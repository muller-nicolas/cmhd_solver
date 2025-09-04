module FFTW_mod
use, intrinsic :: iso_c_binding
implicit none
include "fftw3.f"

public

integer (kind=8) :: plan_for, plan_back

save
contains

SUBROUTINE init_fftw(aa_in, bb_in, N, Nh)
integer, intent(in) :: N, Nh
double precision, intent(in) :: aa_in(N,N)
complex*16, intent(in) :: bb_in(Nh,N)

call dfftw_plan_dft_r2c_2d_(plan_for, N, N, aa_in, bb_in, FFTW_ESTIMATE)
call dfftw_plan_dft_c2r_2d_(plan_back, N, N, bb_in, aa_in, FFTW_ESTIMATE)

RETURN
END SUBROUTINE init_fftw

SUBROUTINE end_fftw
call dfftw_destroy_plan(plan_back)
call dfftw_destroy_plan(plan_for)
END SUBROUTINE end_fftw

END MODULE FFTW_mod
