import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import sys

path = "./"
name = 'divb' # 'rho' or 'jz' or 'wz' or 'divu' or 'divb'
ista = 100
iend = 149
iskip = 1
gamma = 1.4
cs = 1
P = np.loadtxt('out_parameter')
N = int(P[6]) # N
L = 2*np.pi
dx = L/N
k = np.fft.fftfreq(N)*N
ky,kx = np.meshgrid(k,k)

cmap = plt.cm.coolwarm

if len(sys.argv)>1:
    iend = int(sys.argv[1])
if len(sys.argv)>2:
    iskip = int(sys.argv[2])

def load_fields(path, name, time, N):
    if name=="wz":
        fieldx = load_fields(path, 'ux', time, N)
        fieldy = load_fields(path, 'uy', time, N)
        field = derivex(fieldy) - derivey(fieldx)
    elif name=="jz":
        fieldx = load_fields(path, 'bx', time, N)
        fieldy = load_fields(path, 'by', time, N)
        field = derivex(fieldy) - derivey(fieldx)
    elif name=="divu":
        fieldx = load_fields(path, 'ux', time, N)
        fieldy = load_fields(path, 'uy', time, N)
        field = derivex(fieldx) + derivey(fieldy)
    elif name=="divb":
        fieldx = load_fields(path, 'bx', time, N)
        fieldy = load_fields(path, 'by', time, N)
        field = derivex(fieldx) + derivey(fieldy)
    else:
        filename = path + 'field_' + name + f'-2D-{time}'
        field = np.fromfile(filename, dtype=np.float64).reshape((N,N))
    return field

def derivex(field):
    fk = np.fft.fft(field, axis=0)
    dfk = 1j*kx*fk
    out = np.real(np.fft.ifft(dfk, axis=0))
    return out

def derivey(field):
    fk = np.fft.fft(field, axis=1)
    dfk = 1j*ky*fk
    out = np.real(np.fft.ifft(dfk, axis=1))
    return out

import numpy as np

def energy_decomposition(ux, uy):
    """
    Compute incompressible and compressible kinetic energies from a 2D velocity field.
    
    Parameters
    ----------
    ux, uy : 2D numpy arrays
        Velocity components (must have the same shape).
        
    Returns
    -------
    Ei : float
        Incompressible kinetic energy.
    Ec : float
        Compressible kinetic energy.
    """
    assert ux.shape == uy.shape, "ux and uy must have the same shape"
    
    nx, ny = ux.shape

    # Fourier transforms
    ux_hat = np.fft.fft2(ux)
    uy_hat = np.fft.fft2(uy)
    
    # Wavenumber grids
    kx = np.fft.fftfreq(nx).reshape(-1, 1)
    ky = np.fft.fftfreq(ny).reshape(1, -1)
    
    kx, ky = np.meshgrid(ky.flatten(), kx.flatten())  # reorder to match array layout
    ksq = kx**2 + ky**2
    ksq[ksq == 0] = 1.0  # avoid division by zero
    
    # Dot product k·u_hat
    k_dot_u = kx * ux_hat + ky * uy_hat
    
    # Compressible (curl-free) projection
    ux_c_hat = (k_dot_u * kx) / ksq
    uy_c_hat = (k_dot_u * ky) / ksq
    
    # Incompressible (div-free) projection
    ux_i_hat = ux_hat - ux_c_hat
    uy_i_hat = uy_hat - uy_c_hat
    
    # Back to real space
    ux_c = np.fft.ifft2(ux_c_hat).real
    uy_c = np.fft.ifft2(uy_c_hat).real
    ux_i = np.fft.ifft2(ux_i_hat).real
    uy_i = np.fft.ifft2(uy_i_hat).real
    
    # Energies
    Ec = 0.5 * np.mean(ux_c**2 + uy_c**2)
    Ei = 0.5 * np.mean(ux_i**2 + uy_i**2)
    
    return Ei, Ec


nfiles = (iend-ista)//iskip + 1

ekin = np.zeros(nfiles) 
emag = np.zeros(nfiles) 
eint = np.zeros(nfiles)
div = np.zeros(nfiles)
Ei = np.zeros(nfiles)
Ec = np.zeros(nfiles)

for i in range(ista,iend+1,iskip):
    ux = load_fields(path, 'ux', i, N)
    uy = load_fields(path, 'uy', i, N)
    divu = derivex(ux) + derivey(uy)
    Ei[i-ista], Ec[i-ista] = energy_decomposition(ux, uy)
    div[i-ista] = np.mean(divu**2)
    ekin[i-ista] = 0.5*np.mean(ux**2 + uy**2)
    bx = load_fields(path, 'bx', i, N)
    by = load_fields(path, 'by', i, N)
    emag[i-ista] = 0.5*np.mean(bx**2 + by**2)
    rho = load_fields(path, 'rho', i, N)
    eint[i-ista] = cs**2/(gamma*(gamma-1))*np.mean((1+rho)**gamma - 1)

fig, ax = plt.subplots()

t = np.arange(0, nfiles, 1)
etot = ekin + emag + eint
ax.plot(t, ekin, label='E_kin')
ax.plot(t, Ei, label='E_inc')
ax.plot(t, Ec, label='E_com')
ax.plot(t, emag, label='E_mag')
ax.plot(t, eint, label='E_int')
ax.plot(t, etot/3, label='E_tot')

ax.legend()

plt.figure()
plt.plot(div, label='divu^2')

plt.show()
