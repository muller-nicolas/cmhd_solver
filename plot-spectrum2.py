import matplotlib.pyplot as plt
import numpy as np
from pylab import *

figure(figsize=(15,8))

S1 = np.loadtxt('out_spectrumEB-2D-100')
S2 = np.loadtxt('out_spectrumEB-2D-106')
S3 = np.loadtxt('out_spectrumEB-2D-107')
S4 = np.loadtxt('out_spectrumEB-2D-108')
S5 = np.loadtxt('out_spectrumEB-2D-109')
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

R1=lissage(S1,a)
R2=lissage(S2,a)
R3=lissage(S3,a)
R4=lissage(S4,a)
R5=lissage(S5,a)

xp = 0

plt.xlabel('$k$')
plt.ylabel('$E_b(k)$')
plt.grid(True,linestyle=':',which="both")
plt.xlim(1.,N)
plt.ylim(1.e-30,1.e-2)
plt.loglog(x,R1*np.pi*x*x**(xp)/(N*N*N*N),label='t1')
plt.loglog(x,R2*np.pi*x*x**(xp)/(N*N*N*N),label='t2')
plt.loglog(x,R3*np.pi*x*x**(xp)/(N*N*N*N),label='t3')
plt.loglog(x,R4*np.pi*x*x**(xp)/(N*N*N*N),label='t4')
plt.loglog(x,R5*np.pi*x*x**(xp)/(N*N*N*N),label='t5')
#ks = np.array([2.5e0,1.5e1])
#plt.loglog(ks, ks ** (-6.) * 1e5, 'k',linewidth=3.)
#kl = np.array([1.e0,1.e1])
#plt.ylim(1e-12,1e9)

plt.legend()
plt.savefig("Figure-spectrum-Eb.png")
plt.show()
