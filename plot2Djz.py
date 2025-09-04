# from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
# import math
from matplotlib import cm
# from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

time = 119
P = np.loadtxt('out_parameter')
N = int(P[5])
size=1.#2.*np.pi
L = np.loadtxt(f'out_jz-2D-{time}')
Z0 = np.reshape(L,(N,N))
Z1 = pow(abs(Z0)/abs(Z0).max(),1)*Z0/abs(Z0)
#Z1 = pow(abs(Z0)/abs(Z0).max(),3)*Z0/abs(Z0)

X = np.arange(0, size, size/N)
Y = np.arange(0, size, size/N)
X, Y = np.meshgrid(X, Y)

A=100
vmin, vmax = -1/A, 1/A
vmin, vmax = -0.1, 0.1
fig = plt.figure(figsize=(10,8))
#im = plt.imshow(Z1, cmap=cm.coolwarm, extent=[0, size, 0, size],
#    vmin=-1/A, vmax=1/A, interpolation='bicubic')
im = plt.imshow(Z1, cmap=cm.coolwarm, extent=[0, size, 0, size],
    vmin=vmin, vmax=vmax, interpolation='bicubic')

plt.title('$jz(x,y)$')
plt.xlabel('$x$')
plt.ylabel('$y$')
cbar = fig.colorbar(im)

# plt.savefig(f"Figure-Plot2D-jz-{time}.png",dpi=1200)
plt.show()
