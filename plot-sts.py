import matplotlib.pyplot as plt
import numpy as np

ti = 100 # Bounds of data to process
tf = -1 # Bounds of data to process
vmin = -30
vmax = -15
P = np.loadtxt('out_parameter')
N = int(P[6]) # N
kmax = N//3
cs = 0.3
b0 = 1.0
beta = cs**2 / b0**2
disp = (1.e-3 ) / 5

cmap = plt.cm.Blues

filenames_u = ['STS_uxkx', 'STS_uxky'] #, 'STS_uykx', 'STS_uyky']
filenames_b = ['STS_bxkx', 'STS_bxky'] #, 'STS_bykx', 'STS_byky']
filenames_Ak= ['STS_Ak1kx', 'STS_Ak1ky'] #, 'STS_Ak2kx', 'STS_Ak2ky']
filenames = filenames_Ak
kpar = 5
kper = 5 

def spatiotemporal(filename, ti, tf):
    print(filename)
    time = np.loadtxt('STS_time')
    time = time - time[0]
    Nt = len(time)
    sts = np.fromfile(filename)
    sts = sts[::2] + 1j*sts[1::2]
    # filename2 = filename.replace('u','b') # For z+ and z-
    # sts2 = np.fromfile(filename2)
    # sts2 = sts2[::2] + 1j*sts2[1::2]
    # sts = sts - sts2
    if len(sts)==Nt*N:
        sts = sts.reshape(Nt,N)
        FFF = sts[ti:tf,:]
        Omega = np.fft.fftshift(np.fft.fft(FFF, axis=0))
        k = np.fft.fftshift(np.fft.fftfreq(N)*N)
    else:
        sts = sts.reshape(Nt,N//2+1)
        FFF = sts[ti:tf,:]
        Omega = np.fft.fftshift(np.fft.fft(FFF, axis=0), axes=0)
        # k = np.fft.fftshift(np.fft.fftfreq(N//2)*N//2)
        k = np.linspace(0, N//2, N//2+1)
        # print(k)
    sts = sts[ti:tf,:]
    time = time[ti:tf]
    T = time[-1] - time[0]
    # FFF[:,0] /= 100 # atenuate kx=0

    # Omega = (np.fft.fft(FFF, axis=0))
    Omega = np.abs(Omega / np.sqrt(np.sum(np.abs(Omega))**2))**2 # Normalization
    time_omega = time[ti:tf]
    T = time_omega[-1] - time_omega[0]
    w = np.arange(0,(len(time)-1)//2) * 2*np.pi/T # Frequencies
    w = np.concatenate((-np.flip(w[1:]), w))
    return time, Omega, k, w

# Main

fig, axes = plt.subplots(1,2,figsize=(12,7), sharex=True, sharey=True, tight_layout=True)
# ax = [ axes[i,j] for i in range(2) for j in range(2) ]
ax = axes
for i,filename in enumerate(filenames):
    title = filename[4:]
    time, Omega, k, w = spatiotemporal(filename, ti, tf)
    OM = Omega.shape[0] // 2
    OM2 = len(w) // 2
    dT = Omega.shape[0]
    sz = dT // 5

    plot = ax[i].imshow(np.log(Omega[OM-sz:OM+sz+1,:]), extent=(k[0],k[-1],w[OM2-sz],w[OM2+sz]), aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)

    ax[i].set_ylabel(r'$\omega$')

    ax[i].set_title(title + f" k={kpar}")
ax[0].set_xlabel(r'$k_\parallel$')
ax[1].set_xlabel(r'$k_\perp$')

ax[0].set_xlim(left=0, right=kmax+0)

dispersion = lambda disp,kparallel,kperp: disp*(kparallel**3 + kperp**3)
# omega_fast = lambda kperp,kparallel,beta: b0*np.sqrt(kparallel**2 + kperp**2)/np.sqrt(2)*np.sqrt((1+beta) + np.sqrt( (1+beta)**2 - 4*beta*kparallel**2/(kparallel**2 + kperp**2)) + dispersion(disp,kparallel,kperp)**2)
# omega_slow = lambda kperp,kparallel,beta: b0*np.sqrt(kparallel**2 + kperp**2)/np.sqrt(2)*np.sqrt((1+beta) - np.sqrt( (1+beta)**2 - 4*beta*kparallel**2/(kparallel**2 + kperp**2)) + dispersion(disp,kparallel,kperp)**2)
# omega_alfven = lambda kperp,kparallel: np.sqrt(b0**2 * kparallel**2 + dispersion(disp,kparallel,kperp)**2)

omega_fast = lambda kperp,kparallel,beta: b0*np.sqrt(kparallel**2 + kperp**2)/np.sqrt(2)*np.sqrt((1+beta) + np.sqrt( (1+beta)**2 - 4*beta*kparallel**2/(kparallel**2 + kperp**2))) + dispersion(disp,kparallel,kperp)
omega_slow = lambda kperp,kparallel,beta: b0*np.sqrt(kparallel**2 + kperp**2)/np.sqrt(2)*np.sqrt((1+beta) - np.sqrt( (1+beta)**2 - 4*beta*kparallel**2/(kparallel**2 + kperp**2))) + dispersion(disp,kparallel,kperp)
omega_alfven = lambda kperp,kparallel: np.sqrt(b0**2 * kparallel**2) + dispersion(disp,kparallel,kperp)

# omega_fast2 = lambda kperp,kparallel,beta: b0**2*(kparallel**2 + kperp**2)/2*((1+beta) + np.sqrt( (1+beta)**2 - 4*beta*kparallel**2/(kparallel**2 + kperp**2)))
# omega_slow2 = lambda kperp,kparallel,beta: b0**2*(kparallel**2 + kperp**2)/2*((1+beta) - np.sqrt( (1+beta)**2 - 4*beta*kparallel**2/(kparallel**2 + kperp**2)))
# omega_fast = lambda kperp,kparallel,beta: np.sqrt(np.abs(omega_fast2(kperp,kparallel,beta)))
# omega_slow = lambda kperp,kparallel,beta: np.sqrt(np.abs(omega_slow2(kperp,kparallel,beta)))

k = np.linspace(0, kmax, kmax+1)
# x = np.sqrt(k**2 + kper**2)
ax[0].plot(k, omega_fast(kper,k,beta), 'k--', label=r'Fast')
ax[0].plot(k, -omega_fast(kper,k,beta), 'k:')
ax[0].plot(k, omega_slow(kper,k,beta), 'r--', label=r'Slow')
ax[0].plot(k, -omega_slow(kper,k,beta), 'r:')
ax[0].plot(k, omega_alfven(kper,k), 'C1--', label=r'Alfven')
ax[0].plot(k, -omega_alfven(kper,k), 'C1:')

# x = np.sqrt(k**2 + kpar**2)
ax[1].plot(k, omega_fast(k,kpar,beta), 'k--', label=r'Fast')
ax[1].plot(k, -omega_fast(k,kpar,beta), 'k:')
ax[1].plot(k, omega_slow(k,kpar,beta), 'r--', label=r'Slow')
ax[1].plot(k, -omega_slow(k,kpar,beta), 'r:')
# ax[1].plot(k, x, 'C1--', label=r'Alfven')
# ax[1].plot(k, -x, 'C1:')

# ax[0].plot(x, x*np.sqrt(1+beta), 'k--', label=r'$\sqrt{1+\beta}k$ (F)')
# ax[0].plot(x, -x*np.sqrt(1+beta), 'k:')
# ax[1].plot(x, x*np.sqrt(1+beta), 'k--')
# ax[1].plot(x, -x*np.sqrt(1+beta), 'k:')

# ax[0].plot(x, x*np.sqrt(beta), 'r--', label=r'$\sqrt{\beta}k$ (S)')
# ax[0].plot(x, -x*np.sqrt(beta), 'r:')
# ax[1].plot(x, x*np.sqrt(beta), 'r--')
# ax[1].plot(x, -x*np.sqrt(beta), 'r:')

ax[0].legend()
ax[1].legend()

#fig.savefig("figs/spatiotemporal_spectrum.pdf")
plt.show()


