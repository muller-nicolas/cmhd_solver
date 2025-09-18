import matplotlib.pyplot as plt
import numpy as np
import sys

fig, ax = plt.subplots(1,3,figsize=(15,5), sharex=True, sharey=True, tight_layout=True)

ista = 100
iend = 109
iskip = 1
if len(sys.argv)>1:
    iend = int(sys.argv[1])
if len(sys.argv)>2:
    iskip = int(sys.argv[2])

filenames_EU = [f'out_spectrumEU-1D-{i:03d}' for i in range(ista,iend+1,iskip)]
filenames_EB = [f'out_spectrumEB-1D-{i:03d}' for i in range(ista,iend+1,iskip)]
filenames_rho = [f'out_spectrumrho-1D-{i:03d}' for i in range(ista,iend+1,iskip)]

P = np.loadtxt('out_parameter')
kinj = P[3]
N = int(P[6]) # kmax
Nh = N//2+1
kmax = N//3
k = np.arange(0, kmax)

nfiles = len(filenames_EU)
colors = [ plt.cm.viridis(i/nfiles) for i in range(nfiles) ]

for i in range(nfiles):
    filename_u = filenames_EU[i]
    filename_b = filenames_EB[i]
    filename_r = filenames_rho[i]

    Su = np.loadtxt(filename_u)  # load data from file
    Sb = np.loadtxt(filename_b)  # load data from file
    Sr = np.loadtxt(filename_r)  # load data from file

    ax[0].loglog(k[1:kmax],Su[1:kmax],label=f't{i+1}', color=colors[i])
    ax[1].loglog(k[1:kmax],Sb[1:kmax],label=f't{i+1}', color=colors[i])
    ax[2].loglog(k[1:kmax],Sr[1:kmax],label=f't{i+1}', color=colors[i])

ax[0].set_xlabel(r'$k$')
ax[1].set_xlabel(r'$k$')
ax[2].set_xlabel(r'$k$')
ax[0].set_ylabel(r'$E_u(k)$')
ax[1].set_ylabel(r'$E_b(k)$')
ax[2].set_ylabel(r'$E_{\rho}(k)$')

ax[0].set_ylim(1.e-20,1.e-1)

ax[0].axvline(kinj, color='k', linestyle='--')
ax[1].axvline(kinj, color='k', linestyle='--')
ax[2].axvline(kinj, color='k', linestyle='--')

ksta = int(kinj) - 0
kend = kmax
x = k[ksta:kmax]
ax[0].plot(x, (x/x[0])**(-3.)*Su[ksta], color='k', ls='--',linewidth=1.)
ax[1].plot(x, (x/x[0])**(-3.)*Sb[ksta], color='k', ls='--',linewidth=1.)
ax[2].plot(x, (x/x[0])**(-3.)*Sr[ksta], color='k', ls='--',linewidth=1.)

ax[0].legend(ncol=2)

# plt.savefig("Figure-spectrum-Eu.png")
plt.show()
