import matplotlib.pyplot as plt
import numpy as np
import sys

fig, ax = plt.subplots(1,1,figsize=(12,6), sharex=True, sharey=True, tight_layout=True)

ista = 100
iend = 199
iskip = 10
if len(sys.argv)>1:
    iend = int(sys.argv[1])
if len(sys.argv)>2:
    iskip = int(sys.argv[2])

filenames_EU = [f'out_spectrumEU-1D-{i:03d}' for i in range(ista,iend+1,iskip)]

P = np.loadtxt('out_parameter')
kinj = P[3]
N = int(P[6]) # kmax
Nh = N//2+1
kmax = N//3
k = np.arange(0, Nh)

nfiles = len(filenames_EU)
colors = [ plt.cm.viridis(i/nfiles) for i in range(nfiles) ]

for i in range(nfiles):
    filename_u = filenames_EU[i]

    Su = np.loadtxt(filename_u, dtype=np.float64)  # load data from file
    enstrophy = k**2 * Su

    ax.loglog(k[1:kmax],enstrophy[1:kmax],label=f't{i+1}', color=colors[i])

ax.set_xlabel(r'$k$')
ax.set_ylabel(r'$\Omega(k)$')

ax.set_ylim(bottom=1.e-8)

ax.axvline(kinj, color='k', linestyle='--')

ksta = int(kinj) - 1
kend = kmax
x = k[ksta:kmax]
# ax.plot(x, (x/x[0])**(-2.)*Su[ksta], color='k', ls='--',linewidth=1.)
# ax[0].plot(x, (x/x[0])**(-2.)*Su[ksta], color='k', ls='--',linewidth=1.)

# ax[0].legend(ncol=2)

# plt.savefig("Figure-spectrum-Eu.png")
plt.show()
