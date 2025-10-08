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

filenames_Ak1 = [f'out_spectrumAk1-1D-{i:03d}' for i in range(ista,iend+1,iskip)]
filenames_Ak2 = [f'out_spectrumAk2-1D-{i:03d}' for i in range(ista,iend+1,iskip)]
filenames_Ek = [f'out_spectrumEAk-1D-{i:03d}' for i in range(ista,iend+1,iskip)]
filenames_Hk = [f'out_spectrumHAk-1D-{i:03d}' for i in range(ista,iend+1,iskip)]

P = np.loadtxt('out_parameter')
kinj = P[3]
N = int(P[6]) # kmax
Nh = N//2+1
kmax = N//3
k = np.arange(0, kmax)

nfiles = len(filenames_Ak1)
colors = [ plt.cm.viridis(i/nfiles) for i in range(nfiles) ]

fig, axes = plt.subplots(2,2,figsize=(15,7), sharex=True, sharey=True, tight_layout=True)
ax = [ axes[i,j] for i in range(2) for j in range(2) ]

for i in range(nfiles):
    filename_Ak1 = filenames_Ak1[i]
    filename_Ak2 = filenames_Ak2[i]
    filename_EAk = filenames_Ak1[i]
    filename_HAk = filenames_Ak2[i]

    Su = np.loadtxt(filename_Ak1, dtype=np.float64)  # load data from file
    Sb = np.loadtxt(filename_Ak2, dtype=np.float64)  # load data from file
    Se = np.loadtxt(filename_EAk, dtype=np.float64)  # load data from file
    Sh = np.loadtxt(filename_HAk, dtype=np.float64)  # load data from file

    ax[0].loglog(k[1:kmax],Su[1:kmax],label=f't{i+1}', color=colors[i])
    ax[1].loglog(k[1:kmax],Sb[1:kmax],label=f't{i+1}', color=colors[i])
    ax[2].loglog(k[1:kmax],Se[1:kmax],label=f't{i+1}', color=colors[i])
    ax[3].loglog(k[1:kmax],Sh[1:kmax],label=f't{i+1}', color=colors[i])

ax[0].set_xlabel(r'$k$')
ax[1].set_xlabel(r'$k$')
ax[0].set_ylabel(r'$E_{A^+}(k)$')
ax[1].set_ylabel(r'$E_{A^-}(k)$')

ax[0].set_ylim(1.e-11,1.e1)

ax[0].axvline(kinj, color='k', linestyle='--')
ax[1].axvline(kinj, color='k', linestyle='--')

ksta = int(kinj) - 0
kend = kmax
x = k[ksta:kmax]
ax[0].plot(x, (x/x[0])**(-3.)*Su[ksta], color='k', ls='--',linewidth=1.)
ax[1].plot(x, (x/x[0])**(-3.)*Sb[ksta], color='k', ls='--',linewidth=1.)
ax[0].plot(x, (x/x[0])**(-3/2)*Su[ksta], color='k', ls='--',linewidth=1., label=r'$k^{-3/2}$')
ax[1].plot(x, (x/x[0])**(-3/2)*Sb[ksta], color='k', ls='--',linewidth=1., label=r'$k^{-3/2}$')

ax[0].legend(ncol=2)

# plt.savefig("Figure-spectrum-Eu.png")
plt.show()
