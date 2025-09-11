# Makefile
# Compiler and flags
FC = gfortran
FFLAGS = -O3
LDFLAGS = -L/usr/local/Cellar/fftw/3.3.10_1/lib -I/usr/local/Cellar/fftw/3.3.10_1/include
LIBS = -lfftw3 # -lfftw3_threads
PYTHON = /usr/local/opt/python@3.9/bin/python3.9

# Adastra
# LDFLAGS = -L/opt/cray/pe/fftw/3.3.10.6/x86_genoa/lib -I/opt/cray/pe/fftw/3.3.10.6/x86_genoa/include

# Targets
# all: run_cmhd run_spectrum1 run_spectrum2 run_spectrum3
all: CMHD2D

clean:
	rm -f out_* *.mod STS* restart-* *.exe 

plot_spectra: run_spectrum1 run_spectrum2 run_spectrum3

compile_all: CMHD2D 

SRC = parameters.f90 FFTW_mod.f90 adaptive_mod.f90 spectral_mod.f90 cMHD_mod.f90 outputs_mod.f90 CMHD2D.f90

# Build and run CMHD2D
CMHD2D: $(SRC)
	$(FC) $(FFLAGS) $(LDFLAGS) $(LIBS) $(SRC) -o CMHD2D.exe

run_cmhd: CMHD2D
	./CMHD2D

run_spectrum1: 
	$(PYTHON) plot-spectrum.py

run_spectrum2: 
	$(PYTHON) plot-spectrum2.py

run_spectrum3: 
	$(PYTHON) plot-spectrum3.py

