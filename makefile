# Makefile
# Compiler and flags
FC = gfortran
FFLAGS = -O3 -fopenmp -march=native -funroll-loops -ftree-vectorize
DEBUG = #-g -fbacktrace -fcheck=all -Wall
LDFLAGS = -L/usr/local/Cellar/fftw/3.3.10_1/lib -I/usr/local/Cellar/fftw/3.3.10_1/include
LIBS = -lfftw3 -lfftw3_threads -lm
PYTHON = /usr/local/opt/python@3.9/bin/python3.9

# Adastra
# FC = ftn
# FFLAGS = -O2 -fopenmp
# DEBUG = 
# LDFLAGS = -L/opt/cray/pe/fftw/3.3.10.6/x86_genoa/lib -I/opt/cray/pe/fftw/3.3.10.6/x86_genoa/include
# LIBS = -lfftw3 -lfftw3_threads -lm

# Cholesky
# FC = ifx
# FFLAGS = -O3 -qopenmp
# DEBUG = #-O0 -g -check all -fpe0 -warn all -traceback
# LDFLAGS = -I$(FFTW_INC) -L$(FFTW_LIB)
# LIBS = -lfftw3 -lfftw3_threads

# Targets
# all: run_cmhd run_spectrum1 run_spectrum2 run_spectrum3
all: CMHD2D

clean:
	rm -f out_* *.mod STS* restart-* field_* # *.exe 

compile_all: CMHD2D 

SRC = parameters.f90 FFTW_mod.f90 adaptive_mod.f90 spectral_mod.f90 cMHD_mod.f90 outputs_mod.f90 CMHD2D.f90

# Build and run CMHD2D
CMHD2D: $(SRC)
	$(FC) $(FFLAGS) $(LDFLAGS) $(DEBUG) $(LIBS) $(SRC) -o CMHD2D.exe

run_cmhd: CMHD2D
	./CMHD2D.exe

