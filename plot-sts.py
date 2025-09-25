import matplotlib.pyplot as plt
import numpy as np

ti = 100 # Bounds of data to process
tf = -1 # Bounds of data to process
vmin = -20
vmax = -10
P = np.loadtxt('out_parameter')
N = int(P[6]) # N
kmax = N//3

cmap = plt.cm.Blues

filenames_x = ['STS_Euxx_x', 'STS_Euyy_x', 'STS_Ebxx_x', 'STS_Ebyy_x']
filenames_y = ['STS_Euxx_y', 'STS_Euyy_y', 'STS_Ebxx_y', 'STS_Ebyy_y']
filenames_u = ['STS_uxkx', 'STS_uxky', 'STS_uykx', 'STS_uyky']
filenames_b = ['STS_bxkx', 'STS_bxky', 'STS_bykx', 'STS_byky']
filenames_r = ['STS_rhokx', 'STS_rhoky']
filenames = filenames_b

def spatiotemporal(filename, ti, tf):
    print(filename)
    time = np.loadtxt('STS_time')
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
    Omega = np.abs(Omega / np.sqrt(np.sum(np.abs(Omega))**2))**2
    time_omega = time[ti:tf]
    T = time_omega[-1] - time_omega[0]
    w = np.arange(0,(len(time)-1)//2) * 2*np.pi/T
    w = np.concatenate((-np.flip(w[1:]), w))
    return time, Omega, k, w

# Main

fig, axes = plt.subplots(2,2,figsize=(12,7), sharex=True, sharey=True, tight_layout=True)
ax = [ axes[i,j] for i in range(2) for j in range(2) ]
for i,filename in enumerate(filenames):
    title = filename[4:]
    time, Omega, k, w = spatiotemporal(filename, ti, tf)
    OM = Omega.shape[0] // 2
    OM2 = len(w) // 2
    dT = Omega.shape[0]
    sz = dT // 3

    plot = ax[i].imshow(np.log(Omega[OM-sz:OM+sz,:]), extent=(k[0],k[-1],w[OM2-sz],w[OM2+sz]), aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
    # plot = ax[i].imshow(np.log(Omega[OM-sz:OM+sz,:]), extent=(k[0],k[-1],w[OM2-sz],w[OM2+sz]), aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)

    ax[i].set_title(title)
    ax[i].set_xlabel(r'$k$')
    ax[i].set_ylabel(r'$\omega$')

ax[0].set_xlim(left=0, right=kmax+0)
# ax.set_ylim(bottom=-2, top=N//2)
# fig.colorbar(plot)

x = np.linspace(0, kmax, kmax+1)
ax[0].plot(x, x, 'k--', label='k')
ax[1].plot(x, x, 'k--')
ax[2].plot(x, x, 'k--')
ax[3].plot(x, x, 'k--')
ax[0].plot(x, 2*x, 'r--', label='2k')
ax[1].plot(x, 2*x, 'r--')
ax[2].plot(x, 2*x, 'r--')
ax[3].plot(x, 2*x, 'r--')
# ax[0].plot(x, np.sqrt(2)*x, 'C1--', label=r'$\sqrt{2}k$')
# ax[1].plot(x, np.sqrt(2)*x, 'C1--')
# ax[2].plot(x, np.sqrt(2)*x, 'C1--')
# ax[3].plot(x, np.sqrt(2)*x, 'C1--')

ax[0].legend()

#fig.savefig("figs/spatiotemporal_spectrum.pdf")
plt.show()


