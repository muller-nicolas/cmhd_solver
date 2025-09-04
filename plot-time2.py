import matplotlib.pyplot as plt
import numpy as np
from pylab import *

figure(figsize=(10,7))

P = np.loadtxt('out_parameter')
EU = np.loadtxt('out_deltaT')
EB = np.loadtxt('out_nu')
time = np.loadtxt('out_time')

DeltaT = P[0]
Nt = P[1]
N = int(P[1]/P[2])
ki = P[3]
x = np.linspace(0, Nt*DeltaT*ki, N)

plt.subplot(311)
#plt.ylim(1,5.5)
plt.plot(time,EU,label='deltaT')
plt.grid(True,linestyle=':', linewidth=1)
plt.legend()

plt.subplot(312)
plt.plot(time,EB,label='nu')
plt.xlabel('$t$')
plt.grid(True,linestyle=':', linewidth=1)
#plt.ylim(0,0.6)
plt.legend()

plt.subplot(313)
plt.plot(x,time,label='time')
plt.xlabel('$t$')
plt.grid(True,linestyle=':', linewidth=1)
#plt.ylim(0,0.6)
plt.legend()

# plt.savefig("Figure-Time2.png",dpi=600)
plt.show()
