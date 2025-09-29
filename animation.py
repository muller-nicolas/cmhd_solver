import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import sys

path = "./"
name = 'divb' # 'rho' or 'jz' or 'wz' or 'divu' or 'divb'
ista = 100
iend = 109
iskip = 1
P = np.loadtxt('out_parameter')
N = int(P[6]) # N
L = 2*np.pi
dx = L/N
k = np.fft.fftfreq(N)*N
ky,kx = np.meshgrid(k,k)

cmap = plt.cm.coolwarm

if len(sys.argv)>1:
    name = sys.argv[1]
if len(sys.argv)>2:
    iend = int(sys.argv[2])

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
        filename = path + 'out_' + name + f'-2D-{time}'
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

def run_animation():
    anim_running = True

    # Function to allow pause/play with click
    def onClick(event):
        nonlocal anim_running
        if anim_running:
            anim.event_source.stop()
            anim_running = False
        else:
            anim.event_source.start()
            anim_running = True

    def animFunc(frame):
        field = load_fields(path, name, frame, N)
        im.set_array(field)
        ax.set_title(frame)
        return im,

    fig.canvas.mpl_connect('button_press_event', onClick)

    anim = animation.FuncAnimation(fig, animFunc, interval=100, frames=range(ista,iend+1),
            blit=False, repeat_delay=1000, repeat=True)

fig = plt.figure(1)
ax = plt.axes(xlim=(0, N), ylim=(0, N))
field = load_fields(path, name, (ista+iend+2)//2, N)
vmin = -np.max(field) * 1.0
vmax = -vmin
im = plt.imshow(field, animated=False, vmin=vmin, vmax=vmax, cmap=cmap)
plt.title(name)
plt.colorbar()

run_animation()

#ani = animation.FuncAnimation(fig, run_animation, interval=50, blit=True,
#                                repeat_delay=1000,repeat=True)

#ani.save('dipolar_supersolid.gif', writer='imagemagick')

plt.draw()
plt.show()
