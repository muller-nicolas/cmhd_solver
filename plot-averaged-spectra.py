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

filenames_EU = [f'out_spectrumEU-1D-{i:03d}' for i in range(ista,iend+1)]
filenames_Ek = [f'out_spectrumEk-1D-{i:03d}' for i in range(ista,iend+1)]
filenames_EB = [f'out_spectrumEB-1D-{i:03d}' for i in range(ista,iend+1)]
filenames_rho = [f'out_spectrumrho-1D-{i:03d}' for i in range(ista,iend+1)]
filenames_Ak1 = [f'out_spectrumAk1-1D-{i:03d}' for i in range(ista,iend+1)]
filenames_Ak2 = [f'out_spectrumAk2-1D-{i:03d}' for i in range(ista,iend+1)]
filenames_EAk = [f'out_spectrumEAk-1D-{i:03d}' for i in range(ista,iend+1)]
filenames_HAk = [f'out_spectrumHAk-1D-{i:03d}' for i in range(ista,iend+1)]

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

nfiles = len(filenames_EU)
colors = [ plt.cm.viridis(i/nfiles) for i in range(nfiles) ]

Su = mean_spectrum(filenames_EU)
Sk = mean_spectrum(filenames_Ek)
Sb = mean_spectrum(filenames_EB)
Sr = mean_spectrum(filenames_rho)

SA1 = mean_spectrum(filenames_Ak1)
SA2 = mean_spectrum(filenames_Ak2)
SAe = mean_spectrum(filenames_EAk)
SAh = mean_spectrum(filenames_HAk)

ax[0].loglog(k[1:kmax],Su[1:kmax]/1,label= "Eu")
ax[0].loglog(k[1:kmax],Sb[1:kmax],label= "Eb")
ax[0].loglog(k[1:kmax],Sr[1:kmax],label= "Erho")
ax[0].loglog(k[1:kmax],Sk[1:kmax],label= "Ekin")

ax[1].loglog(k[1:kmax],SA1[1:kmax],label="EA1")
ax[1].loglog(k[1:kmax],SA2[1:kmax],label="EA2")
ax[1].loglog(k[1:kmax],SAe[1:kmax],label="EAe")
ax[1].loglog(k[1:kmax],np.abs(SAh[1:kmax]),label="EAh")

ax[0].set_xlabel(r'$k$')
ax[1].set_xlabel(r'$k$')
ax[0].set_ylabel(r'$E(k)$')
ax[1].set_ylabel(r'$E(k)$')

# ax[0].set_ylim(1.e-8,1.e2)
ax[0].set_ylim(bottom=1.e-8)

ax[0].axvline(kinj, color='k', linestyle='--')
ax[1].axvline(kinj, color='k', linestyle='--')

ksta = int(kinj) - 1
kend = kmax
x = k[ksta:kmax]
ax[0].plot(x, (x/x[0])**(-2.)*Su[ksta], color='k', ls='--',linewidth=1., label=r'$k^{-2}$')
ax[0].plot(x, (x/x[0])**(-2)*Sb[ksta], color='k', ls='--',linewidth=1.)
ax[0].plot(x, (x/x[0])**(-5/2.)*Su[ksta], color='r', ls='--',linewidth=1., label=r'$k^{-5/2}$')

ax[1].plot(x, (x/x[2])**(-2.)*SAe[ksta+2], color='k', ls='--',linewidth=1., label=r'$k^{-2}$')

ax[0].legend(ncol=1)
ax[1].legend(ncol=1)

# plt.savefig("Figure-averaged-spectrum.png")
plt.show()
