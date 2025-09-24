MODULE parameters
implicit none

public
integer, parameter :: N = 128 ! Resolution
integer, parameter :: Nh = N/2+1  ! Number of wavenumbers
integer, parameter :: Na = int(N/3)   
integer :: nthreads = 1  ! Number of omp threads
double precision, parameter :: deltaT0 = 2.5d-3
double precision, parameter :: nu0  = 2.d-5 ! viscosity
double precision, parameter :: alpha  = 0.d-2 ! damping (large-scale)
double precision :: eta = nu0 ! magnetic viscosity

integer :: istore_sp = 100
integer :: istore_fields = 100
integer, parameter :: ndeltaT = 20000
integer, parameter :: inrj = 20
integer, parameter :: ispec = 200    !*********must be a multiple of inrj
integer, parameter :: ifields = 200  !*********must be a multiple of inrj
double precision, parameter :: kinj = 2.
double precision, parameter :: disp = 0.d-3  ! without dispersion => 0.d0
double precision, parameter :: amp = 1.d-2
double precision, parameter :: famp = 1.d-1   !*********forcing => famp > 0.d0
double precision, parameter :: width = 5.   !*********forcing width (gaussian forcing)
double precision, parameter :: corr0 = 1.d-1   !*********forcing correlation 
integer, parameter :: sts = 1 ! 1 if spatio-temporal spectra are saved, else 0
integer, parameter :: ists = 10 ! time step for spatio-temporal spectra
integer, parameter :: nrestart = 0    !*********for a restart => nrestart = 1
integer, parameter :: seed = 100  ! seed for random number generator

! Constants
double precision, parameter :: pi = 3.141592653589793238d0
double complex :: imag = (0.0d0,1.0d0)
double precision :: dk = 1.

save

end module parameters
