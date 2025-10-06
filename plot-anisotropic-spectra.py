import matplotlib.pyplot as plt
import numpy as np
import sys

fig, ax = plt.subplots(1,2,figsize=(12,6), sharex=True, sharey=True, tight_layout=True)

ista = 150
iend = 199
if len(sys.argv)>1:
    ista = int(sys.argv[1])
if len(sys.argv)>2:
    iend = int(sys.argv[2])

filenames_EUpara = [f'out_spectrumEUpara-1D-{i:03d}' for i in range(ista,iend+1)]
filenames_EUperp = [f'out_spectrumEUperp-1D-{i:03d}' for i in range(ista,iend+1)]
filenames_EBpara = [f'out_spectrumEBpara-1D-{i:03d}' for i in range(ista,iend+1)]
filenames_EBperp = [f'out_spectrumEBperp-1D-{i:03d}' for i in range(ista,iend+1)]

nfiles = len(filenames_EUpara)

P = np.loadtxt('out_parameter')
kinj = P[3]
N = int(P[6]) # kmax
Nh = N//2+1
kmax = N//3
k = np.arange(0, kmax)

def mean_spectrum(filenames):
    nfiles = len(filenames)
    spec = np.zeros(Nh)
    for i in range(nfiles):
        filename = filenames[i]
        spec = spec + np.loadtxt(filename, dtype=np.float128)  # load data from file
    return spec/nfiles

colors = [ plt.cm.viridis(i/nfiles) for i in range(nfiles) ]

Supara = mean_spectrum(filenames_EUpara)
Superp = mean_spectrum(filenames_EUperp)
Sbpara = mean_spectrum(filenames_EBpara)
Sbperp = mean_spectrum(filenames_EBperp)

ax[0].loglog(k[1:kmax],Supara[1:kmax]/1,label= "Eu")
ax[0].loglog(k[1:kmax],Sbpara[1:kmax],label= "Eb")

ax[1].loglog(k[1:kmax],Superp[1:kmax]/1,label= "Eu")
ax[1].loglog(k[1:kmax],Sbperp[1:kmax],label= "Eb")

ax[0].set_xlabel(r'$k_\parallel$')
ax[1].set_xlabel(r'$k_\perp$')
ax[0].set_ylabel(r'$E(k)$')
ax[1].set_ylabel(r'$E(k)$')

# ax[0].set_ylim(1.e-8,1.e2)
ax[0].set_ylim(bottom=1.e-8)

ax[0].axvline(kinj, color='k', linestyle='--')
ax[1].axvline(kinj, color='k', linestyle='--')

ksta = int(kinj) - 1
kend = kmax
x = k[ksta:kmax]
ax[0].plot(x, (x/x[0])**(-3)*Sbpara[ksta], color='k', ls='--',linewidth=1., label=r'$k^{-3}$')
ax[0].plot(x, (x/x[0])**(-2.)*Sbpara[ksta], color='r', ls='--',linewidth=1., label=r'$k^{-2}$')
# ax[0].plot(x, (x/x[0])**(-2)*Supara[ksta], color='k', ls='--',linewidth=1.)
# ax[0].plot(x, (x/x[0])**(-5/2.)*Supara[ksta], color='r', ls='--',linewidth=1., label=r'$k^{-5/2}$')

ax[1].plot(x, (x/x[2])**(-5/3)*Sbperp[ksta+2], color='k', ls='--',linewidth=1., label=r'$k^{-5/3}$')
ax[1].plot(x, (x/x[2])**(-3/2)*Sbperp[ksta+2], color='r', ls='--',linewidth=1., label=r'$k^{-3/2}$')
ax[1].plot(x, (x/x[2])**(-7/3)*Sbperp[ksta+2], color='b', ls='--',linewidth=1., label=r'$k^{-7/3}$')
# ax[1].plot(x, (x/x[2])**(-2)*Sbperp[ksta+2], color='violet', ls='--',linewidth=1., label=r'$k^{-2}$')

ax[0].legend(ncol=1)
ax[1].legend(ncol=1)

# plt.savefig("Figure-averaged-spectrum.png")
plt.show()
