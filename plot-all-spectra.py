import matplotlib.pyplot as plt
import numpy as np
import sys

fig, ax = plt.subplots(1,3,figsize=(15,5), sharex=True, sharey=True, tight_layout=True)

ista = 100
iend = 199
iskip = 10
if len(sys.argv)>1:
    iend = int(sys.argv[1])
if len(sys.argv)>2:
    iskip = int(sys.argv[2])

filenames_EU = [f'out_spectrumEU-1D-{i:03d}' for i in range(ista,iend+1,iskip)]
filenames_EB = [f'out_spectrumEB-1D-{i:03d}' for i in range(ista,iend+1,iskip)]
filenames_rho = [f'out_spectrumrho-1D-{i:03d}' for i in range(ista,iend+1,iskip)]
filenames_EUinc = [f'out_spectrumEUinc-1D-{i:03d}' for i in range(ista,iend+1,iskip)]
filenames_EUcom = [f'out_spectrumEUcom-1D-{i:03d}' for i in range(ista,iend+1,iskip)]

P = np.loadtxt('out_parameter')
kinj = P[3]
N = int(P[6]) # kmax
Nh = N//2+1
kmax = N//3
k = np.arange(0, kmax)

nfiles = len(filenames_EU)
colors = [ plt.cm.viridis(i/nfiles) for i in range(nfiles) ]

ekin = np.zeros(nfiles) 
emag = np.zeros(nfiles) 
eint = np.zeros(nfiles)

def load_and_plot_spectrum(ax, filenames, ls='-'):
    nfiles = len(filenames)
    colors = [ plt.cm.viridis(i/nfiles) for i in range(nfiles) ]
    energy = np.zeros(nfiles)
    for i in range(nfiles):
        filename = filenames[i]
        spec = np.loadtxt(filename, dtype=np.float64)  # load data from file
        energy[i] = np.sum(spec[1:])
        ax.loglog(k[1:kmax],spec[1:kmax], label=f't{i+1}', color=colors[i], ls=ls)

    return energy, spec

ekin, Su = load_and_plot_spectrum(ax[0], filenames_EU)
emag, Sb = load_and_plot_spectrum(ax[1], filenames_EB)
eint, Sr = load_and_plot_spectrum(ax[2], filenames_rho)
einc, Suinc = load_and_plot_spectrum(ax[0], filenames_EUinc, ls='--')
ecom, Sucom = load_and_plot_spectrum(ax[0], filenames_EUcom, ls=':')

ax[0].set_xlabel(r'$k$')
ax[1].set_xlabel(r'$k$')
ax[2].set_xlabel(r'$k$')
ax[0].set_ylabel(r'$E_u(k)$')
ax[1].set_ylabel(r'$E_b(k)$')
ax[2].set_ylabel(r'$E_{\rho}(k)$')

# ax[0].set_ylim(1.e-11,1.e1)
ax[0].set_ylim(bottom=1.e-10)

ax[0].axvline(kinj, color='k', linestyle='--')
ax[1].axvline(kinj, color='k', linestyle='--')
ax[2].axvline(kinj, color='k', linestyle='--')

ksta = int(kinj) - 1
kend = kmax
x = k[ksta:kmax]
ax[0].plot(x, (x/x[0])**(-2.)*Su[ksta], color='k', ls='--',linewidth=1., label=r'$k^{-2}$')
ax[1].plot(x, (x/x[0])**(-2.)*Sb[ksta], color='k', ls='--',linewidth=1.)
ax[2].plot(x[2:], (x[2:]/x[2])**(-3.)*Sr[ksta+2], color='k', ls='--',linewidth=1.)
# ax[0].plot(x, (x/x[0])**(-2.)*Su[ksta], color='k', ls='--',linewidth=1.)
ax[1].plot(x, (x/x[0])**(-2.)*Sb[ksta], color='k', ls='--',linewidth=1.)
ax[2].plot(x, (x/x[0])**(-2.)*Sr[ksta], color='k', ls='--',linewidth=1.)

ax[0].legend(ncol=2)

print(Su[0])
print(Sb[0])
print(Sr[0])

plt.figure()
plt.plot(ekin, label='E_kin')
plt.plot(einc, label='E_inc')
plt.plot(ecom, label='E_com')
plt.plot(emag, label='E_mag')
# plt.plot(eint, label='E_int')
etot = ekin + emag #+ eint
plt.plot(etot/2, label='E_tot')
plt.legend()

# plt.savefig("Figure-spectrum-Eu.png")
plt.show()
