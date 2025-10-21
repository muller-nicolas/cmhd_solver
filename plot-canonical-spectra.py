import matplotlib.pyplot as plt
import numpy as np
import sys

ista = 100
iend = 109
iskip = 1
if len(sys.argv)>1:
    iend = int(sys.argv[1])
if len(sys.argv)>2:
    iskip = int(sys.argv[2])

filenames_Ek = [f'out_spectrumEAk-1D-{i:03d}' for i in range(ista,iend+1,iskip)]
filenames_Hk = [f'out_spectrumHAk-1D-{i:03d}' for i in range(ista,iend+1,iskip)]

P = np.loadtxt('out_parameter')
kinj = P[3]
N = int(P[6]) # kmax
Nh = N//2+1
kmax = N//3
k = np.arange(0, kmax)

nfiles = len(filenames_Ek)
colors = [ plt.cm.viridis(i/nfiles) for i in range(nfiles) ]

fig, axes = plt.subplots(1,2,figsize=(15,7), sharex=True, sharey=True, tight_layout=True)
# ax = [ axes[i,j] for i in range(2) for j in range(2) ]
ax = axes

for i in range(nfiles):
    filename_EAk = filenames_Ek[i]
    filename_HAk = filenames_Hk[i]

    Se = np.loadtxt(filename_EAk, dtype=np.float64)  # load data from file
    Sh = np.loadtxt(filename_HAk, dtype=np.float64)  # load data from file
    Sh = np.abs(Sh)

    ax[0].loglog(k[1:kmax],Se[1:kmax],label=f't{i+1}', color=colors[i])
    ax[1].loglog(k[1:kmax],Sh[1:kmax],label=f't{i+1}', color=colors[i])

ax[0].set_xlabel(r'$k$')
ax[1].set_xlabel(r'$k$')
ax[0].set_ylabel(r'$E(k)$')
ax[1].set_ylabel(r'$|H(k)|$')

ax[0].set_ylim(bottom=1.e-11)

ax[0].axvline(kinj, color='k', linestyle='--')
ax[1].axvline(kinj, color='k', linestyle='--')

ksta = int(kinj) - 0
kend = kmax
x = k[ksta:kmax]
ax[0].plot(x, (x/x[0])**(-2.)*Se[ksta], color='k', ls='--',linewidth=1., label=r'$k^{-2}$')
ax[1].plot(x, (x/x[0])**(-2.)*Sh[ksta], color='k', ls='--',linewidth=1., label=r'$k^{-2}$')
ax[0].plot(x, (x/x[0])**(-3/2)*Se[ksta], color='r', ls='--',linewidth=1., label=r'$k^{-3/2}$')
ax[1].plot(x, (x/x[0])**(-3/2)*Sh[ksta], color='r', ls='--',linewidth=1., label=r'$k^{-3/2}$')

ax[0].legend(ncol=2)

# plt.savefig("Figure-spectrum-Eu.png")
plt.show()
