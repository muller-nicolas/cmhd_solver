MODULE parameters
implicit none

public
integer, parameter :: N = 128
integer, parameter :: Nh = N/2+1
integer, parameter :: Na = N/3  ! partie entière
double precision, parameter :: deltaT0 = 2.5d-3
double precision, parameter :: deltaTi = deltaT0/10.d0
double precision, parameter :: nu0  = 8.d-9 ! viscosity
double precision, parameter :: alpha  = 1.d0 ! damping (large-scale)
double precision :: eta = nu0 ! magnetic viscosity

integer :: istore_sp = 100
integer :: istore_fields = 100
integer, parameter :: ndeltaT = 10000
integer, parameter :: inrj = 100
integer, parameter :: ispec = 1000  !*********must be a multiple of inrj
integer, parameter :: ifields = 1000  !*********must be a multiple of inrj
! integer :: irestart = 1000
double precision, parameter :: kinj = 3.
double precision, parameter :: disp = 0.d-2  ! without dispersion => 0.d0
double precision, parameter :: amp = 0.d-3
double precision, parameter :: off = 1.d0   !*********forcing => off = 1.d0
double precision, parameter :: width = 5.   !*********forcing width (gaussian forcing)
integer, parameter :: sts = 1 ! 1 if spatio-temporal spectra are saved, else 0
integer, parameter :: ists = 10 ! time step for spatio-temporal spectra
integer, parameter :: nrestart = 1    !*********for a restart => nrestart = 0
integer, parameter :: seed = 100  ! seed for random number generator

! Constants
double precision, parameter :: pi = 3.141592653589793238d0
double complex :: imag = (0.0d0,1.0d0)
double precision :: dk = 1.

save

end module parameters
