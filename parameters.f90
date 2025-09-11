MODULE parameters
implicit none

public
integer, parameter :: N = 256
integer, parameter :: Nh = N/2+1
integer, parameter :: Na = N/3  ! partie entière
double precision, parameter :: deltaT0 = 5.d-4
double precision, parameter :: deltaTi = deltaT0/10.d0
double precision, parameter :: nu0  = 5.d-8 ! viscosity
double precision :: eta = nu0 ! magnetic viscosity

integer :: istore_sp = 100
integer :: istore_fields = 100
integer, parameter :: ndeltaT = 1000
integer, parameter :: inrj = 1
integer, parameter :: ispec = 100  !*********must be a multiple of inrj
integer, parameter :: ifields = 100  !*********must be a multiple of inrj
! integer :: irestart = 1000
double precision, parameter :: kinj = 3.
double precision, parameter :: disp = 0.d-5  ! without dispersion => 0.d0
! double precision, parameter :: a = 1.d0  !*********a=0. => linear equations; a=1. non-linear equations
double precision, parameter :: amp = 5.d-2
double precision, parameter :: off = 0.d-2   !*********forcing => off = 1.d0
integer, parameter :: sts = 0 ! 1 if spatio-temporal spectra are saved, else 0
integer, parameter :: ists = 10 ! time step for spatio-temporal spectra
integer, parameter :: nrestart = 1    !*********for a restart => nrestart = 0
integer, parameter :: seed = 800  ! seed for random number generator

! Constants
double precision, parameter :: pi = 3.141592653589793238d0
double complex :: imag = (0.0d0,1.0d0)
double precision :: dk = 2.*pi

save

end module parameters
