import matplotlib.pyplot as plt
import numpy as np
from pylab import *

fig, ax = subplots(1,3,figsize=(15,6), sharex=True, sharey=True, tight_layout=True)

ista = 100
iend = 119
iskip = 2
filenames_EU = [f'out_spectrumEU-1D-{i:03d}' for i in range(ista,iend+1,iskip)]
filenames_EB = [f'out_spectrumEB-1D-{i:03d}' for i in range(ista,iend+1,iskip)]
filenames_rho = [f'out_spectrumrho-1D-{i:03d}' for i in range(ista,iend+1,iskip)]

P = np.loadtxt('out_parameter')
N = int(P[5]/2+1)
x = np.linspace(1, N, N)

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

norm = np.pi*x*x**(xp)/(N*N*N*N)
for i in range(len(filenames_EU)):
    filename_u = filenames_EU[i]
    filename_b = filenames_EB[i]
    filename_r = filenames_rho[i]

    Su = np.loadtxt(filename_u)  # load data from file
    Sb = np.loadtxt(filename_b)  # load data from file
    Sr = np.loadtxt(filename_r)  # load data from file

    Ru=lissage(Su,a)*norm
    Rb=lissage(Sb,a)*norm
    Rr=lissage(Sr,a)*norm

    ax[0].loglog(x,Ru,label=f't{i+1}')
    ax[1].loglog(x,Rb,label=f't{i+1}')
    ax[2].loglog(x,Rr,label=f't{i+1}')

ax[0].set_xlabel(r'$k$')
ax[1].set_xlabel(r'$k$')
ax[2].set_xlabel(r'$k$')
ax[0].set_ylabel(r'$E_u(k)$')
ax[1].set_ylabel(r'$E_b(k)$')
ax[2].set_ylabel(r'$E_{\rho}(k)$')

# plt.grid(True,linestyle=':',which="both")
ax[0].set_xlim(right=N)
ax[0].set_ylim(1.e-30,1.e-2)

#ks = np.array([2.5e0,1.5e1])
#plt.loglog(ks, ks ** (-6.) * 1e5, 'k',linewidth=3.)
#kl = np.array([1.e0,1.e1])
#plt.ylim(1e-12,1e9)

# plt.legend()
# plt.savefig("Figure-spectrum-Eu.png")
plt.show()
