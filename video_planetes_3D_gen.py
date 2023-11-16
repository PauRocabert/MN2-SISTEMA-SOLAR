import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import datetime 

plt.style.use('dark_background')

file = open(r'trayectorias.txt','r')
Lines = file.readlines()


list_headers = Lines[0].split(' ')

N_cossos = 7
N_dim = 3
N_dies = len(Lines)-1

DADES = np.ndarray(shape=(N_cossos,N_dim, N_dies))
for cos in range(N_cossos):
    DADES[cos] = np.transpose(np.array([[float(line.split(' ')[4*cos+1]), float(line.split(' ')[4*cos+2]),float(line.split(' ')[4*cos+3])]for line in Lines[1:]]))
        


fig = plt.figure()
ax = fig.add_subplot(projection="3d")
colors_planetes =['yellow', 'darkorange','lime', 'darkturquoise', 'crimson', 'violet', 'cyan']
colors = ['white', 'orange','limegreen', 'mediumturquoise', 'red','darkmagenta', 'navy']

PLANETES = [ax.plot([],[], 'o', markersize= 5, color = colors_planetes[i])[0] for i in range(N_cossos-2)]
traces = [ax.plot([],[], '-', lw = 0.7, alpha = 0.9, color = colors[i])[0] for i in range(N_cossos-2) ]
patches = PLANETES + traces


def animate(i):
    i = int(i)
    for cos, object in enumerate(zip(traces, PLANETES)):
        object[0].set_data(DADES[cos][0:2,:i])
        object[0].set_3d_properties(DADES[cos][2,:i])
        object[1].set_data(DADES[cos][0][i], DADES[cos][1][i])
        object[1].set_3d_properties(DADES[cos][2][i])
    return   patches



ax.set(xlim3d=(-1.5, 1.5), xlabel='X')
ax.set(ylim3d=(-1.5, 1.5), ylabel='Y')
ax.set(zlim3d=(-0.3, 0.3), zlabel='Z')

fig.patch.set_facecolor('black') 
ani = animation.FuncAnimation(fig, animate, frames = N_dies,interval = 1, blit = True)
dpi = 200
writer = animation.writers['ffmpeg'](fps=30)
ani.save('video_planetes_trayectorias_estandard.mp4',writer=writer,dpi=dpi)
