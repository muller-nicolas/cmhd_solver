import matplotlib.pyplot as plt
import numpy as np
from pylab import *

figure(figsize=(15,8))

filenames = ['out_spectrumEB-1D-101',
            #  'out_spectrumEB-1D-106', 
            #  'out_spectrumEB-1D-107',
            #  'out_spectrumEB-1D-108',
            #  'out_spectrumEB-1D-109'
            ]

ista = 100
iend = 101
iskip = 1
filenames = [f'out_spectrumEB-1D-{i:03d}' for i in range(ista,iend+1,iskip)]

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

for i,filename in enumerate(filenames):
    S = np.loadtxt(filename)
    R=lissage(S,a)
    plt.loglog(x,R*np.pi*x*x**(xp)/(N*N*N*N),label=f't{i+1}')

plt.xlabel('$k$')
plt.ylabel('$E_b(k)$')
plt.grid(True,linestyle=':',which="both")
plt.xlim(1.,N)
plt.ylim(1.e-30,1.e-2)
#ks = np.array([2.5e0,1.5e1])
#plt.loglog(ks, ks ** (-6.) * 1e5, 'k',linewidth=3.)
#kl = np.array([1.e0,1.e1])
#plt.ylim(1e-12,1e9)

plt.legend()
plt.savefig("Figure-spectrum-Eb.png")
plt.show()
