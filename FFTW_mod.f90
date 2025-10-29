module FFTW_mod
use, intrinsic :: iso_c_binding
use omp_lib
use parameters
implicit none
include "fftw3.f"

public

integer (kind=8) :: plan_for, plan_back
integer, parameter :: Nh = N/2+1  ! Number of wavenumbers
integer, parameter :: Na = int(N/3)   ! Dealiasing cutoff

save
contains

SUBROUTINE init_fftw
double precision, allocatable :: aa(:,:)
double complex, allocatable :: bb(:,:)
integer :: void

allocate (aa(N,N))
allocate (bb(Nh,N))

! Initialize FFTW threading
call dfftw_init_threads(void)
call dfftw_plan_with_nthreads(nthreads)

! Create plans
call dfftw_plan_dft_r2c_2d_(plan_for, N, N, aa, bb, FFTW_MEASURE)
call dfftw_plan_dft_c2r_2d_(plan_back, N, N, bb, aa, FFTW_MEASURE)

deallocate (aa, bb)

RETURN
END SUBROUTINE init_fftw

SUBROUTINE FFT_PS(Ain, Aout)
double precision, intent(in) :: Ain(N,N)
double complex, intent(out) :: Aout(Nh,N)

call dfftw_execute_dft_r2c(plan_for, Ain, Aout)

RETURN
END SUBROUTINE FFT_PS

SUBROUTINE FFT_SP(Ain, Aout)
double complex, intent(in) :: Ain(Nh,N)
double precision, intent(out) :: Aout(N,N)
double precision :: norm

norm = 1. / real(N*N)

call dfftw_execute_dft_c2r(plan_back, Ain, Aout)
Aout = Aout*norm

RETURN
END SUBROUTINE FFT_SP


SUBROUTINE end_fftw
call dfftw_destroy_plan(plan_back)
call dfftw_destroy_plan(plan_for)
call dfftw_cleanup_threads()
END SUBROUTINE end_fftw

END MODULE FFTW_mod
