import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import sys

path = "./"
name = 'b' # 'rho' or 'jz' or 'wz' or 'divu' or 'divb'
ista = 100
iend = 109
iskip = 1
cs = 0.3
b0 = 1

P = np.loadtxt('out_parameter')
N = int(P[6]) # N
L = 2*np.pi
dx = L/N
k = np.fft.fftfreq(N)*N
ky,kx = np.meshgrid(k,k)

cmap = plt.cm.plasma

if len(sys.argv)>1:
    name = sys.argv[1]
if len(sys.argv)>2:
    iend = int(sys.argv[2])
if len(sys.argv)>3:
    iskip = int(sys.argv[3])

def load_fields(path, name, time, N):
    if name=="u":
        filename = path + f'field_ux-2D-{time}'
        fieldx = np.fromfile(filename, dtype=np.float64).reshape((N,N))
        filename = path + f'field_uy-2D-{time}'
        fieldy = np.fromfile(filename, dtype=np.float64).reshape((N,N))
    elif name=="b":
        filename = path + f'field_bx-2D-{time}'
        fieldx = np.fromfile(filename, dtype=np.float64).reshape((N,N))
        filename = path + f'field_by-2D-{time}'
        fieldy = np.fromfile(filename, dtype=np.float64).reshape((N,N))
    else:
        print("ERROR: no proper fields were selected")
    return fieldx, fieldy

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
        fieldx, fieldy = load_fields(path, name, frame, N)
        field = fieldx**2 + fieldy**2
        # fluct = np.mean(field)
        if name=='u':
            print(np.max(np.sqrt(field))/cs)
        if name=='b':
            print(np.max(np.sqrt(field))/b0)
        im.set_array(field)
        ax.set_title(frame)
        return im,

    fig.canvas.mpl_connect('button_press_event', onClick)

    anim = animation.FuncAnimation(fig, animFunc, interval=100, frames=range(ista,iend+1,iskip),
            blit=False, repeat_delay=1000, repeat=True)

nfiles = (iend-ista)//iskip + 1

fig = plt.figure(1)
ax = plt.axes(xlim=(0, N), ylim=(0, N))
fieldx, fieldy = load_fields(path, name, ista+(iend-ista)//2, N)
field = fieldx**2 + fieldy**2
vmin = 0
vmax = np.max(field) * 1.0
im = plt.imshow(field, animated=False, vmin=vmin, vmax=vmax, cmap=cmap)
plt.title(name)
plt.colorbar()

run_animation()

#ani = animation.FuncAnimation(fig, run_animation, interval=50, blit=True,
#                                repeat_delay=1000,repeat=True)

#ani.save('dipolar_supersolid.gif', writer='imagemagick')

plt.draw()
plt.show()
