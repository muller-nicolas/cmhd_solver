import matplotlib.pyplot as plt
import numpy as np

ti = 10 # Bounds of data to process
tf = -1 # Bounds of data to process
vmin = -80
vmax = -10
P = np.loadtxt('out_parameter')
N = int(P[6]) # N

cmap = plt.cm.Blues

filenames_x = ['STS_Euxx_x', 'STS_Euyy_x', 'STS_Ebxx_x', 'STS_Ebyy_x']
filenames_y = ['STS_Euxx_y', 'STS_Euyy_y', 'STS_Ebxx_y', 'STS_Ebyy_y']
filenames = filenames_y

def spatiotemporal(filename, ti, tf):
    time = np.loadtxt('STS_time')
    Nt = len(time)
    sts = np.fromfile(filename)
    if len(sts)==Nt*N:
        sts = sts.reshape(Nt,N)
    else:
        sts = sts.reshape(Nt,N//2+1)
    sts = sts[ti:tf,:]
    time = time[ti:tf]
    T = time[-1] - time[0]
    FFF = sts[ti:tf,:]
    # FFF[:,0] /= 100 # atenuate kx=0

    Omega = np.fft.fftshift(np.fft.fft(FFF, axis=0))
    Omega = np.abs(Omega / np.sqrt(np.sum(np.abs(Omega))**2))**2
    time_omega = time[ti:tf]
    T = time_omega[-1] - time_omega[0]
    k = np.fft.fftshift(np.fft.fftfreq(N)*N)
    w = np.arange(0,(len(time)-1)//2) * 2*np.pi/T
    w = np.concatenate((-np.flip(w[1:]), w))
    return time, Omega, k, w

# Main

fig, axes = plt.subplots(2,2,figsize=(12,7), sharex=True, sharey=True, tight_layout=True)
ax = [ axes[i,j] for i in range(2) for j in range(2) ]
for i,filename in enumerate(filenames):
    time, Omega, k, w = spatiotemporal(filename, ti, tf)
    OM = Omega.shape[0] // 2
    OM2 = len(w) // 2
    dT = Omega.shape[0]
    sz = dT // 10

    plot = ax[i].imshow(np.log(Omega[OM-sz:OM+sz,:]), extent=(k[0],k[-1],w[OM2-sz],w[OM2+sz]), aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)

    ax[i].set_title(filename)
    ax[i].set_xlabel(r'$k$')
    ax[i].set_ylabel(r'$\omega$')

ax[0].set_xlim(left=0, right=N//3+0)
# ax.set_ylim(bottom=-2, top=N//2)
# fig.colorbar(plot)

#fig.savefig("figs/spatiotemporal_spectrum.pdf")
plt.show()


