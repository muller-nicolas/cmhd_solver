MODULE parameters
implicit none

public
! Simulation parameters
integer, parameter :: N = 128 ! Resolution
integer, parameter :: Nh = N/2+1  ! Number of wavenumbers
integer, parameter :: Na = int(N/3)   
integer :: nthreads = 4  ! Number of omp threads
double precision, parameter :: deltaT0 = 5.d-3
integer, parameter :: nrestart = 0    !*********for a restart => nrestart = 1
integer, parameter :: seed = 100  ! seed for random number generator NOT WORKING YET
integer :: istore_sp = 100
integer :: istore_fields = 100

! Physical parameters
double precision, parameter :: nu0  = 1.d-5 ! hyperviscosity (n=2)
! double precision, parameter :: nu0  = 5.d-3 ! viscosity (n=1)
double precision, parameter :: alpha  = 0.d-2 ! damping (large-scale)
double precision, parameter :: cspeed  = 0.3   ! speed of sound
double precision, parameter :: gamma  = 1.4   ! polytropic index
double precision :: eta = nu0 ! magnetic viscosity
integer, parameter :: ndeltaT = 10000
integer, parameter :: inrj = 100
integer, parameter :: ispec = 100    !*********must be a multiple of inrj
integer, parameter :: ifields = 100  !*********must be a multiple of inrj

! IC and forcing parameters
double precision, parameter :: kinj = 2.
double precision, parameter :: disp = 0.d-4  ! without dispersion => 0.d0
double precision, parameter :: amp = 1.d-1
double precision, parameter :: famp = 1.d-1   !*********forcing => famp > 0.d0
double precision, parameter :: width = 5.     !*********forcing width (gaussian forcing)
double precision, parameter :: corr0 = 1.d0  !*********forcing correlation 
integer, parameter :: sts = 1   ! 1 if spatio-temporal spectra are saved, else 0
integer, parameter :: ists = 10 ! time step for spatio-temporal spectra

! Constants
double precision, parameter :: pi = 3.141592653589793238d0
double complex :: imag = (0.0d0,1.0d0)
double precision :: dk = 1.

save

end module parameters
