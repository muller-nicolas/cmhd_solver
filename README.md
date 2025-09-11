# 2D Compressible MHD code

This is a pseudo-spectral code for the 2D compressible MHD equations. 
The main approximations are:
- small $\beta$ (hydrodynamic pressure gradients are neglected)
- small density perturbations (a Taylor expansion is used)

## Compilation

### Dependencies

- FFTW
- Python

### Compilation

```bash
make
```

## Running

```bash
./CMHD2D.exe
```

## Output

The code produces several files in the `out_` directory:

- `out_parameter`: contains the parameters of the simulation
- `out_energy`: contains the energy of the fields
- `out_spectrumEU-1D-101`, `out_spectrumEU-1D-102`, ... : contains the spectrum of the energy of the velocity field
- `out_spectrumEB-1D-101`, `out_spectrumEB-1D-102`, ... : contains the spectrum of the energy of the magnetic field
- `out_spectrumrho-1D-101`, `out_spectrumrho-1D-102`, ... : contains the spectrum of the density
- `out_wz-2D-101`, `out_wz-2D-102`, ... : contains the vorticity 
- `out_jz-2D-101`, `out_jz-2D-102`, ... : contains the current density
- `out_divb-2D-101`, `out_divb-2D-102`, ... : contains the divergence of the magnetic field (Check that is zero)
- `out_divu-2D-101`, `out_divu-2D-102`, ... : contains the divergence of the velocity field
- `out_rho-2D-101`, `out_rho-2D-102`, ... : contains the density
- `out_deltaT`: contains the time step
- `out_nu`: contains the viscosity
- `out_time`: contains the time

## Plotting

Plots scripts are generated using the Python library `matplotlib`.
