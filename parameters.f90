MODULE parameters
implicit none

public
integer, parameter :: N = 256
integer, parameter :: Nh = N/2+1
integer, parameter :: Na = N/3  ! partie entière
double precision, parameter :: deltaT = 1.d-3
double precision, parameter :: deltaTi = deltaT/10.d0
double precision, parameter :: nu = 1.d-7

integer :: istore_sp = 100
integer :: istore_fields = 100
integer, parameter :: ndeltaT = 1000
integer, parameter :: inrj = 1
integer, parameter :: ispec = 100  !*********must be a multiple of inrj
integer, parameter :: ifields = 100  !*********must be a multiple of inrj
! integer :: irestart = 1000
double precision, parameter :: kinj = 3.
double precision, parameter :: eta = nu
double precision, parameter :: disp = 0.d-5  ! without dispersion => 0.d0
double precision, parameter :: a = 1.d0  !*********a=0. => linear equations; a=1. non-linear equations
double precision, parameter :: amp = 1.d-3
double precision, parameter :: off = 0.d0   !*********forcing => off = 1.d0
integer, parameter :: sts = 0 ! 1 if spatio-temporal spectra are saved, else 0
integer, parameter :: ists = 10 ! time step for spatio-temporal spectra
integer, parameter :: nrestart = 1    !*********for a restart => nrestart = 0
integer, parameter :: seed = 800

! Constants
double precision, parameter :: pi = 3.141592653589793238d0
double complex :: imag = (0.0d0,1.0d0)
double precision :: dk = 2.*pi

save

end module parameters
