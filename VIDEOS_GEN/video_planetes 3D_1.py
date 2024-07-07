import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D


file = open(r'sistema_solar_posicions_postopt.txt','r')
Lines = file.readlines()

list_headers = Lines[0].split(' ')
print(Lines[0][1])

sol         =np.array([[float(line.split(' ')[list_headers.index('x_Sol')]), float(line.split(' ')[list_headers.index('y_Sol')]), float(line.split(' ')[list_headers.index('z_Sol')])] for line in Lines[1:]])
terra       =np.array([[float(line.split(' ')[list_headers.index('x_Terra')]), float(line.split(' ')[list_headers.index('y_Terra')]), float(line.split(' ')[list_headers.index('z_Terra')])] for line in Lines[1:]])
mart        =np.array([[float(line.split(' ')[list_headers.index('x_Mart')]), float(line.split(' ')[list_headers.index('y_Mart')]), float(line.split(' ')[list_headers.index('z_Mart')])] for line in Lines[1:]])
jupiter     =np.array([[float(line.split(' ')[list_headers.index('x_Jupiter')]), float(line.split(' ')[list_headers.index('y_Jupiter')]), float(line.split(' ')[list_headers.index('z_Jupiter')])] for line in Lines[1:]])
sol_T     = np.transpose(sol)
terra_T   = np.transpose(terra)
mart_T    = np.transpose(mart)
jupiter_T = np.transpose(jupiter)


fig = plt.figure()
ax = fig.add_subplot(projection="3d")
colors = ['#fdaa48', '#b0ff9d', '#d0fefe', '#ffffe4' ]
tr2, = ax.plot([],[],[], '-', lw = 1.5, alpha = 0.9, color = colors[1]) 
tr3, = ax.plot([],[],[], '-', lw = 1.5, alpha = 0.9, color = colors[2])
tr4, = ax.plot([],[],[], '-', lw = 1.5, alpha = 0.9, color = colors[3])

SOL, = ax.plot([],[],[], 'o', color = 'orange')
TERRA, = ax.plot([],[],[], 'o', color = 'blue') 
MART, = ax.plot([],[],[], 'o', color = 'red') 
JUPITER, = ax.plot([],[],[], 'o', color = 'orange') 

def animate(i):
    tr2.set_data(terra_T[0:2,:i])
    tr2.set_3d_properties(terra_T[2,:i])
    tr3.set_data(mart_T[0:2,:i])
    tr3.set_3d_properties(mart_T[2,:i])
    tr4.set_data(jupiter_T[0:2,:i])
    tr4.set_3d_properties(jupiter_T[2,:i])
    SOL.set_data(sol_T[0][i], sol_T[1][i])
    SOL.set_3d_properties(sol_T[2][i])
    TERRA.set_data(terra_T[0][i], terra_T[1][i])
    TERRA.set_3d_properties(terra_T[2][i])
    MART.set_data(mart_T[0][i], mart_T[1][i])
    MART.set_3d_properties(mart_T[2][i])
    JUPITER.set_data(jupiter_T[0][i],jupiter_T[1][i])
    JUPITER.set_3d_properties(jupiter_T[2][i]) 
    return tr2, tr3, tr4, SOL, TERRA, MART, JUPITER,

    

# Setting the axes properties
ax.set(xlim3d=(-5.5, 5.5), xlabel='X')
ax.set(ylim3d=(-5.5, 5.5), ylabel='Y')
ax.set(zlim3d=(-0.1, 0.1), zlabel='Z')

ax.patch.set_facecolor('black') 
ani = animation.FuncAnimation(fig, animate,frames = len(sol_T[0]) , blit=False)
dpi = 200
writer = animation.writers['ffmpeg'](fps=30)
ani.save('video_planetes3D_noaxis.mp4',writer=writer,dpi=dpi)
