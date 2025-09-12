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

cmap = plt.cm.coolwarm

if len(sys.argv)>1:
    name = sys.argv[1]
if len(sys.argv)>2:
    iend = int(sys.argv[2])

def load_fields(path, name, time, N):
    filename = path + 'out_' + name + f'-2D-{time}'
    field = np.fromfile(filename, dtype=np.float64).reshape((N,N))
    return field

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
field = load_fields(path, name, (ista+iend)//2, N)
vmin = np.min(field) * 0.8
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
