import matplotlib.pyplot as plt
import numpy as np
from pylab import *

figure(figsize=(10,8))

P = np.loadtxt('out_parameter')
EU = np.loadtxt('out_EU')
EB = np.loadtxt('out_EB')
Erho = np.loadtxt('out_Erho')
divu = np.loadtxt('out_divu')
divb = np.loadtxt('out_divb')
time = np.loadtxt('out_time')

DeltaT = P[0]
Nt = P[1]
N = int(P[1]/P[2])
ki = P[3]
x = np.linspace(0, Nt*DeltaT*ki, N)
x = time

plt.gcf().subplots_adjust(left=0.1,bottom=0.1,right=0.9,
                          top=0.9,wspace=0.2,hspace=0.3)
                          
plt.subplot(211)
plt.title('Energies')
#plt.ylim(0,6.e-7)
#plt.xlim(0,3.9)
plt.plot(x,Erho,label='<\u03C1$u^2/2$>')
plt.plot(x,EB,label='<$b^2/2$>')
plt.plot(x,(Erho+EB)/2.,label='$E^{tot}/2$')
plt.xlabel('$t$')
plt.grid(True,linestyle=':', linewidth=1)
plt.legend()

plt.subplot(212)
plt.title('Compressibilities')
#plt.ylim(0.,1.1)
#plt.xlim(0,3.9)
plt.plot(x,divu,label='<div u>')
plt.plot(x,divb,label='<div b>')
plt.xlabel('$t$')
plt.grid(True,linestyle=':', linewidth=1)
plt.legend()

plt.savefig("Figure-Time.png",dpi=600)
plt.show()
