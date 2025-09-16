MODULE parameters
implicit none

public
integer, parameter :: N = 512 ! Resolution
integer, parameter :: Nh = N/2+1  ! Number of wavenumbers
integer, parameter :: Na = int(N/3)   
integer :: nthreads = 4  ! Number of omp threads
double precision, parameter :: deltaT0 = 1.d-3
! double precision, parameter :: deltaTi = deltaT0/10.d0
double precision, parameter :: nu0  = 1.d-7 ! viscosity
double precision, parameter :: alpha  = 1.d-1 ! damping (large-scale)
double precision :: eta = nu0 ! magnetic viscosity

integer :: istore_sp = 100
integer :: istore_fields = 100
integer, parameter :: ndeltaT = 10
integer, parameter :: inrj = 10
integer, parameter :: ispec = 10  !*********must be a multiple of inrj
integer, parameter :: ifields = 10  !*********must be a multiple of inrj
! integer :: irestart = 1000
double precision, parameter :: kinj = 3.
double precision, parameter :: disp = 0.d-2  ! without dispersion => 0.d0
double precision, parameter :: amp = 0.d-2
double precision, parameter :: famp = 1.d-5   !*********forcing => famp > 0.d0
double precision, parameter :: width = 3.   !*********forcing width (gaussian forcing)
double precision, parameter :: corr0 = 1.d-3   !*********forcing correlation 
integer, parameter :: sts = 0 ! 1 if spatio-temporal spectra are saved, else 0
integer, parameter :: ists = 10 ! time step for spatio-temporal spectra
integer, parameter :: nrestart = 1    !*********for a restart => nrestart = 0
integer, parameter :: seed = 100  ! seed for random number generator

! Constants
double precision, parameter :: pi = 3.141592653589793238d0
double complex :: imag = (0.0d0,1.0d0)
! double precision :: dk = 2*pi
double precision :: dk = 1.

save

end module parameters
