# Makefile
# Compiler and flags
FC = gfortran
FFLAGS = -O3
LDFLAGS = -L/usr/local/Cellar/fftw/3.3.10_1/lib -I/usr/local/Cellar/fftw/3.3.10_1/include
LIBS = -lfftw3 # -lfftw3_threads
PYTHON = /usr/local/opt/python@3.9/bin/python3.9

# Targets
all: run_cmhd run_spectrum1 run_spectrum2 run_spectrum3

clean:
	rm -f out* *.mod restart-* CMHD2D-v8 spectrum-anim spectrum-anim2 spectrum-anim3

spectra: run_spectrum1 run_spectrum2 run_spectrum3

SRC = FFTW_mod.f90 CMHD2D-v8.f90

# Build and run CMHD2D-v8
CMHD2D-v8: FFTW_mod.f90 CMHD2D-v8.f90
	$(FC) $(FFLAGS) $(LDFLAGS) $(LIBS) $(SRC) -o CMHD2D-v8 

run_cmhd: CMHD2D-v8
	./CMHD2D-v8

# Spectrum 1
spectrum-anim: spectrum-anim.f90
	$(FC) $(FFLAGS) -o $@ $<

run_spectrum1: spectrum-anim
	./spectrum-anim &
	$(PYTHON) plot-spectrum.py

# Spectrum 2
spectrum-anim2: spectrum-anim2.f90
	$(FC) $(FFLAGS) -o $@ $<

run_spectrum2: spectrum-anim2
	./spectrum-anim2
	$(PYTHON) plot-spectrum2.py

# Spectrum 3
spectrum-anim3: spectrum-anim3.f90
	$(FC) $(FFLAGS) -o $@ $<

run_spectrum3: spectrum-anim3
	./spectrum-anim3
	$(PYTHON) plot-spectrum3.py

