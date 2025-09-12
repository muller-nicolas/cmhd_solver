import numpy as np
import matplotlib.pyplot as plt

# --- PARAMETERS TO SET ---
filename = "STS_Ebxx_x"  # The binary file written by your Fortran code
N = 128                  # Number of points in ky direction (len(Ebyy(1,:)))
# N = N//2+1
nt = None                # Set to None to infer from file size
dtype = np.float64       # Fortran double precision

# --- LOAD DATA ---
data = np.fromfile(filename, dtype=dtype)

# Infer number of timesteps if not specified
if nt is None:
    if data.size % N != 0:
        raise ValueError("File size not divisible by N, check N or Fortran write format.")
    nt = data.size // N

# Reshape into (nt, N)
spectrum = data.reshape((nt, N))
sps = np.fft.fft(spectrum, axis=0)
spectrum = np.abs(np.fft.fftshift(sps))**2

# --- PLOT ---
plt.figure(figsize=(10, 6))
extent = [0, N-1, 0, nt-1]  # ky index (x), time index (y)
plt.imshow(np.log(spectrum), aspect='auto', origin='lower', extent=extent, cmap='inferno')
plt.colorbar(label='Magnetic Energy Spectrum')
plt.xlabel("ky index")
plt.ylabel("Time step")
plt.title("Spatiotemporal Spectrum (kx=0, ky range)")
plt.tight_layout()
plt.show()
