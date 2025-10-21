import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import sys

path = "./"
name = 'divb' # 'rho' or 'jz' or 'wz' or 'divu' or 'divb'
ista = 100
iend = 199
iskip = 1
gamma = 1.4
cs = 1
P = np.loadtxt('out_parameter')
N = int(P[6]) # N
L = 2*np.pi
dx = L/N
k = np.fft.fftfreq(N)*N
ky,kx = np.meshgrid(k,k)

cmap = plt.cm.coolwarm

if len(sys.argv)>1:
    iend = int(sys.argv[1])
if len(sys.argv)>2:
    iskip = int(sys.argv[2])

def load_fields(path, name, time, N):
    if name=="wz":
        fieldx = load_fields(path, 'ux', time, N)
        fieldy = load_fields(path, 'uy', time, N)
        field = derivex(fieldy) - derivey(fieldx)
    elif name=="jz":
        fieldx = load_fields(path, 'bx', time, N)
        fieldy = load_fields(path, 'by', time, N)
        field = derivex(fieldy) - derivey(fieldx)
    elif name=="divu":
        fieldx = load_fields(path, 'ux', time, N)
        fieldy = load_fields(path, 'uy', time, N)
        field = derivex(fieldx) + derivey(fieldy)
    elif name=="divb":
        fieldx = load_fields(path, 'bx', time, N)
        fieldy = load_fields(path, 'by', time, N)
        field = derivex(fieldx) + derivey(fieldy)
    else:
        filename = path + 'field_' + name + f'-2D-{time}'
        field = np.fromfile(filename, dtype=np.float64).reshape((N,N))
    return field

def derivex(field):
    fk = np.fft.fft(field, axis=0)
    dfk = 1j*kx*fk
    out = np.real(np.fft.ifft(dfk, axis=0))
    return out

def derivey(field):
    fk = np.fft.fft(field, axis=1)
    dfk = 1j*ky*fk
    out = np.real(np.fft.ifft(dfk, axis=1))
    return out

nfiles = (iend-ista)//iskip + 1

ekin = np.zeros(nfiles) 
emag = np.zeros(nfiles) 
eint = np.zeros(nfiles)
div = np.zeros(nfiles)

for i in range(ista,iend+1,iskip):
    ux = load_fields(path, 'ux', i, N)
    uy = load_fields(path, 'uy', i, N)
    divu = derivex(ux) + derivey(uy)
    div[i-ista] = np.mean(divu**2)
    ekin[i-ista] = 0.5*np.mean(ux**2 + uy**2)
    bx = load_fields(path, 'bx', i, N)
    print(bx)
    by = load_fields(path, 'by', i, N)
    emag[i-ista] = 0.5*np.mean(bx**2 + by**2)
    rho = load_fields(path, 'rho', i, N)
    eint[i-ista] = cs**2/(gamma*(gamma-1))*np.mean((1+rho)**gamma - 1)

fig, ax = plt.subplots()

t = np.arange(0, nfiles, 1)
etot = ekin + emag + eint
print(len(t))
print(len(ekin))
ax.plot(t, ekin, label='E_kin')
ax.plot(t, emag, label='E_mag')
ax.plot(t, eint, label='E_int')
ax.plot(t, etot/3, label='E_tot')

ax.legend()

plt.figure()
plt.plot(div, label='divu^2')

plt.show()
