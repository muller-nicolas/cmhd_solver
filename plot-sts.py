import matplotlib.pyplot as plt
import numpy as np

N = 128
ti = 10 # Bounds of data to process
tf = -1 # Bounds of data to process
vmin = -100
vmax = -15

cmap = plt.cm.Blues

filename = 'STS_Euxx_x'
time = np.loadtxt('STS_time')
Nt = len(time)
sts = np.fromfile(filename)
# st = sts[::2] + 1j*sts[1::2]
st = sts[:] 
sts = st.reshape(Nt,N)
print(Nt)

time = time[ti:]
sts = sts[ti:,:]
T = time[-1] - time[0]

FFF = sts[ti:tf,:]
# FFF[:,0] /= 100 # atenuate kx=0
#print(np.max(FFF))
#print(np.min(FFF))
#FFF[:,0] = 0 # set kz=0 to 0

Omega = np.fft.fftshift(np.fft.fft(FFF, axis=0))
Omega = np.abs(Omega / np.sqrt(np.sum(np.abs(Omega))**2))**2
time_omega = time[ti:tf]
T = time_omega[-1] - time_omega[0]
k = np.fft.fftshift(np.fft.fftfreq(N)*N)
w = np.arange(0,(len(time)-1)//2) * 2*np.pi/T
#print(w.shape)
w = np.concatenate((-np.flip(w[1:]), w))
#print(w.shape)
#print(k.shape)
#print(np.max(Omega))
#print(np.min(Omega))
#print(Omega.shape)
#print(w.shape)

OM = Omega.shape[0] // 2
OM2 = len(w) // 2
dT = Omega.shape[0]
sz = dT // 10
#plt.imshow(np.log(Omega), extent=(k[0], k[-1], w[0], w[-1]), aspect=.02)

fig, ax = plt.subplots(figsize=(8,5))

plot = ax.imshow(np.log(Omega[OM-sz:OM+sz,:]), extent=(k[0],k[-1],w[OM2-sz],w[OM2+sz]), aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
ax.set_xlabel(r'$k$')
# ax.set_title(labels[i])

ax.set_ylabel(r'$\omega$')
# ax.set_xlim(left=0, right=N//3-1)
# ax.set_ylim(bottom=-2, top=N//2)
#fig.colorbar(plot)

fig.tight_layout()
#fig.savefig("figs/spatiotemporal_spectrum.pdf")
plt.show()


