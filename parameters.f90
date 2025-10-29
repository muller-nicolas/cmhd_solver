MODULE parameters
implicit none

public

! Simulation parameters
integer, parameter :: N = 128 ! Resolution
double precision, parameter :: deltaT0 = 2.d-3
integer, parameter :: nrestart = 0    !*********for a restart => nrestart = 1
integer :: istore_sp = 100
integer :: istore_fields = 100
integer, parameter :: seed = 100  ! seed for random number generator NOT WORKING YET
integer :: nthreads = 4  ! Number of omp threads

! Physical parameters
logical, parameter :: hyperviscosity = .true. ! hyperviscosity = .true. for \nabla^4 dissipation
double precision, parameter :: nu0  = 1.d-5 ! hyperviscosity (n=2)
! double precision, parameter :: nu0  = 4.d-3 ! viscosity (n=1)
double precision, parameter :: alpha  = 0.d-2 ! damping (large-scale)
double precision, parameter :: cspeed  = 0.3   ! speed of sound
double precision, parameter :: gamma  = 1.4   ! polytropic index
double precision :: eta = nu0 ! magnetic viscosity
integer, parameter :: ndeltaT = 100000
integer, parameter :: inrj = 500
integer, parameter :: ispec = 2000    !*********must be a multiple of inrj
integer, parameter :: ifields = 2000  !*********must be a multiple of inrj

! IC and forcing parameters
double precision, parameter :: kinj = 2.
double precision, parameter :: disp = 1.d-3  ! without dispersion => 0.d0
double precision, parameter :: amp = 0.d-2
double precision, parameter :: famp = 1.d-1   !*********forcing => famp > 0.d0
double precision, parameter :: width = 1.     !*********forcing width (gaussian forcing)
double precision, parameter :: corr0 = 1.d0   !*********forcing correlation 
double precision, parameter :: CompAmp = 0.d0  !********* Compressible amplitude: 0.d0 => incompressible forcing, 1.d0 => compressible forcing
integer, parameter :: sts = 1   ! 1 if spatio-temporal spectra are saved, else 0
integer, parameter :: ists = 10 ! time step for spatio-temporal spectra

! Constants
double precision, parameter :: pi = 3.141592653589793238d0
double complex :: imag = (0.0d0,1.0d0)
double precision :: dk = 1.

save

end module parameters
