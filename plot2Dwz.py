import matplotlib.pyplot as plt
import numpy as np

time = 136
P = np.loadtxt('out_parameter')
N = int(P[6])
size=1.#2.*np.pi

def load_fields(path, name, time, N):
    filename = path + 'out_' + name + f'-2D-{time}'
    field = np.fromfile(filename, dtype=np.float64).reshape((N,N))
    return field

Z0 = load_fields("./", "wz", time, N)
Z1 = pow(abs(Z0)/abs(Z0).max(),1)*Z0/abs(Z0)
#Z1 = pow(abs(Z0)/abs(Z0).max(),3)*Z0/abs(Z0)

X = np.arange(0, size, size/N)
Y = np.arange(0, size, size/N)
X, Y = np.meshgrid(X, Y)

A=100
vmin, vmax = -1/A, 1/A
vmin, vmax = -0.0001, 0.0001
fig = plt.figure(figsize=(10,8))
#im = plt.imshow(Z1, cmap=cm.coolwarm, extent=[0, size, 0, size],
#    vmin=-1/A, vmax=1/A, interpolation='bicubic')
im = plt.imshow(Z1, cmap=plt.cm.coolwarm, extent=[0, size, 0, size],
    vmin=vmin, vmax=vmax, interpolation='bicubic')

plt.title('$wz(x,y)$')
plt.xlabel('$x$')
plt.ylabel('$y$')
cbar = fig.colorbar(im)

# plt.savefig(f"Figure-Plot2D-wz-{time}.png",dpi=1200)
plt.show()
