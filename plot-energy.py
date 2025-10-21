import matplotlib.pyplot as plt
import numpy as np
# from pylab import *

plt.figure(figsize=(10,8), tight_layout=True)

P = np.loadtxt('out_parameter')
EU = np.loadtxt('out_EU')
EB = np.loadtxt('out_EB')
Erho = np.loadtxt('out_Erho')
Eint = np.loadtxt('out_Eint')
divu = np.loadtxt('out_divu')
divb = np.loadtxt('out_divb')
time = np.loadtxt('out_time')

Etot = Erho + EB + Eint
# Etot = Erho + EB #+ Eint
cs = 1.0
b0 = 1

Nt = len(EU)
urms = np.mean(np.sqrt(EU[Nt//2:]))
brms = np.mean(np.sqrt(EB[Nt//2:]))
krms = np.mean(np.sqrt(Erho[Nt//2:]))
print('urms',urms)
print('brms',brms)
print('krms',krms)
print('Ms', urms/cs)
print('Ma', urms/b0)

DeltaT = P[0]
Nt = P[1]
N = int(P[1]/P[2])
ki = P[3]
x = np.linspace(0, Nt*DeltaT*ki, N)
x = time
                          
plt.subplot(211)
plt.title('Energies')
# plt.plot(x,EU,label='<$u^2/2$>')
# plt.plot(x,Urms,label='<$u^2/2$>')
plt.plot(x,Erho,label=r'<$\rho u^2/2$>')
plt.plot(x,EB,label='<$b^2/2$>')
plt.plot(x,Eint,label=r'<$\rho e$>')
plt.plot(x,Etot/3,label='$E^{tot}/3$')
plt.xlabel('$t$')
plt.grid(True,linestyle=':', linewidth=1)
plt.legend()

plt.subplot(212)
plt.title('Compressibilities')
plt.plot(x,divu,label='<div u>')
plt.plot(x,divb,label='<div b>')
plt.xlabel('$t$')
plt.grid(True,linestyle=':', linewidth=1)
plt.legend()

# Etot = Erho + EB
print(Etot[-1]/Etot[0])

# plt.savefig("Figure-Time.png",dpi=600)
plt.show()
