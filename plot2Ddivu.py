import matplotlib.pyplot as plt
import numpy as np

time = 119
P = np.loadtxt('out_parameter')
N = int(P[6])
size=2.*np.pi

def load_fields(path, name, time, N):
    filename = path + 'out_' + name + f'-2D-{time}'
    field = np.fromfile(filename, dtype=np.float64).reshape((N,N))
    return field

Z0 = load_fields("./", "divu", time, N)
#Z1 = pow(abs(Z0)/abs(Z0).max(),1)*Z0/abs(Z0)
#Z1 = pow(abs(Z0)/abs(Z0).max(),3)*Z0/abs(Z0)

X = np.arange(0, size, size/N)
Y = np.arange(0, size, size/N)
X, Y = np.meshgrid(X, Y)

A=1000
fig = plt.figure(figsize=(10,8))
#im = plt.imshow(Z1, cmap=cm.coolwarm, extent=[0, size, 0, size],
#    vmin=-1/A, vmax=1/A, interpolation='bicubic')
im = plt.imshow(Z0, cmap=cm.coolwarm, extent=[0, size, 0, size])
#    ,vmin=-1/A, vmax=1/A, interpolation='bicubic')

plt.title('$div_u(x,y)$')
plt.xlabel('$x$')
plt.ylabel('$y$')
cbar = fig.colorbar(im)

# plt.savefig(f"Figure-Plot2D-divu-{time}.png",dpi=1200)
plt.show()
