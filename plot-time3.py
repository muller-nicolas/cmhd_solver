import matplotlib.pyplot as plt
import numpy as np
from pylab import *

figure(figsize=(10,8))

P = np.loadtxt('out_parameter')
Erho2 = np.loadtxt('out_Erho2')
#uxm = np.loadtxt('out_uxm')
#bxm = np.loadtxt('out_bxm')
#uxdm = np.loadtxt('out_uxdm')

DeltaT = P[0]
Nt = P[1]
N = int(P[1]/P[2])
ki = P[3]
x = np.linspace(0, Nt*DeltaT*ki, N)

plt.gcf().subplots_adjust(left=0.1,bottom=0.1,right=0.9,
                          top=0.9,wspace=0.2,hspace=0.3)
                          
plt.subplot(211)
#plt.ylim(0,6.e-7)
#plt.xlim(0,3.9)
plt.plot(x,Erho2,label='mass')
plt.xlabel('$t$')
plt.grid(True,linestyle=':', linewidth=1)
plt.legend()

plt.subplot(212)
#plt.ylim(0.,1.1)
#plt.xlim(0,3.9)
#plt.plot(x,uxm,label='uxm')
#plt.plot(x,bxm,label='bxm')
#plt.plot(x,uxdm,label='uxdm')
plt.xlabel('$t$')
plt.grid(True,linestyle=':', linewidth=1)
plt.legend()

# plt.savefig("Figure-Time3.png",dpi=600)
plt.show()
