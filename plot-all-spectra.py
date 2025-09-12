import matplotlib.pyplot as plt
import numpy as np
import sys

fig, ax = plt.subplots(1,3,figsize=(15,5), sharex=True, sharey=True, tight_layout=True)

ista = 100
iend = 109
iskip = 2
if len(sys.argv)>1:
    iend = int(sys.argv[1])

filenames_EU = [f'out_spectrumEU-1D-{i:03d}' for i in range(ista,iend+1,iskip)]
filenames_EB = [f'out_spectrumEB-1D-{i:03d}' for i in range(ista,iend+1,iskip)]
filenames_rho = [f'out_spectrumrho-1D-{i:03d}' for i in range(ista,iend+1,iskip)]

P = np.loadtxt('out_parameter')
kmax = int(P[6]/2+1) # kmax
k = np.linspace(1, kmax, kmax)

def lissage(S,L):
    res = np.copy(S) # duplication des valeurs
    for i in range (2,len(S)-2): # toutes les valeurs sauf la première et la dernière
        L_g = min(i,L) # nombre de valeurs disponibles à gauche
        L_d = min(len(S)-i-1,L) # nombre de valeurs disponibles à droite
        Li=min(L_g,L_d)
        res[i]=np.sum(S[i-Li:i+Li+1])/(2*Li+1)
    return res

a=2
xp = 0

nfiles = len(filenames_EU)
colors = [ plt.cm.viridis(i/nfiles) for i in range(nfiles) ]

# TODO: Check this normalization, as well as the lissage function
norm = np.pi*k*k**(xp)/(kmax*kmax*kmax*kmax)
for i in range(nfiles):
    filename_u = filenames_EU[i]
    filename_b = filenames_EB[i]
    filename_r = filenames_rho[i]

    Su = np.loadtxt(filename_u)  # load data from file
    Sb = np.loadtxt(filename_b)  # load data from file
    Sr = np.loadtxt(filename_r)  # load data from file

    # Ru=lissage(Su,a)*norm
    # Rb=lissage(Sb,a)*norm
    # Rr=lissage(Sr,a)*norm
    Ru=Su
    Rb=Sb
    Rr=Sr

    ax[0].loglog(k,Ru,label=f't{i+1}', color=colors[i])
    ax[1].loglog(k,Rb,label=f't{i+1}', color=colors[i])
    ax[2].loglog(k,Rr,label=f't{i+1}', color=colors[i])

ax[0].set_xlabel(r'$k$')
ax[1].set_xlabel(r'$k$')
ax[2].set_xlabel(r'$k$')
ax[0].set_ylabel(r'$E_u(k)$')
ax[1].set_ylabel(r'$E_b(k)$')
ax[2].set_ylabel(r'$E_{\rho}(k)$')

# plt.grid(True,linestyle=':',which="both")
ax[0].set_xlim(right=kmax)
ax[0].set_ylim(1.e-20,1.e-1)

ax[0].axvline(kmax*2/3, color='k', linestyle='--')
ax[1].axvline(kmax*2/3, color='k', linestyle='--')
ax[2].axvline(kmax*2/3, color='k', linestyle='--')

ksta = 1
kend = kmax
x = k[ksta:kmax]
ax[0].plot(x, (x/x[0])**(-6.)*Ru[ksta], color='k', ls='--',linewidth=1.)
ax[1].plot(x, (x/x[0])**(-6.)*Rb[ksta], color='k', ls='--',linewidth=1.)
ax[2].plot(x, (x/x[0])**(-6.)*Rr[ksta], color='k', ls='--',linewidth=1.)

ax[0].legend(ncol=2)

# plt.savefig("Figure-spectrum-Eu.png")
plt.show()
