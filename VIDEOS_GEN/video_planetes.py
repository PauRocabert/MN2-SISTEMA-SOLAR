import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import datetime 


file = open(r'sistema_solar_3D_video.txt','r')
Lines = file.readlines()

dia1= '2021-09-13'
dia_inici = datetime.datetime.strptime(dia1, '%Y-%m-%d')

list_headers = Lines[0].split(' ')
print(Lines[0])

time        =np.array([float(line.split(' ')[list_headers.index('t')]) for line in Lines[1:]])
sol         =np.array([[float(line.split(' ')[list_headers.index('x_Sol')]), float(line.split(' ')[list_headers.index('y_Sol')])]for line in Lines[1:]])
terra       =np.array([[float(line.split(' ')[list_headers.index('x_Terra')]), float(line.split(' ')[list_headers.index('y_Terra')])] for line in Lines[1:]])
mart        =np.array([[float(line.split(' ')[list_headers.index('x_Mart')]), float(line.split(' ')[list_headers.index('y_Mart')]), ] for line in Lines[1:]])
jupiter     =np.array([[float(line.split(' ')[list_headers.index('x_Jupiter')]), float(line.split(' ')[list_headers.index('y_Jupiter')])] for line in Lines[1:]])

sol_T       = np.transpose(sol)
terra_T     = np.transpose(terra)
mart_T      = np.transpose(mart)
jupiter_T   = np.transpose(jupiter)

fig = plt.figure()
ax = fig.add_subplot(111, xlim=(-6,6), ylim=(-6, 6))
ax.set_aspect('equal', adjustable='box')
colors = ['#fdaa48', '#b0ff9d', '#d0fefe', '#ffffe4' ]


SOL, = ax.plot([],[], 'o', color = 'white', markevery = [False, True])
TERRA, = ax.plot([],[], 'o', color = 'blue', markevery = [False, True]) 
MART, = ax.plot([],[], 'o', color = 'red', markevery = [False, True]) 
JUPITER, = ax.plot([],[], 'o', color = 'orange', markevery = [False, True]) 
tr2, = plt.plot([],[], '-', lw = 0.7, alpha = 0.9, color = colors[0]) 
tr3, = plt.plot([],[], '-', lw = 0.7, alpha = 0.9, color = colors[1])
tr4, = plt.plot([],[], '-', lw = 0.7, alpha = 0.9, color = colors[2])



def init():
    SOL.set_data([],[])
    TERRA.set_data([],[])
    MART.set_data([],[])
    JUPITER.set_data([],[])
    tr2.set_data([],[])
    tr3.set_data([],[])
    tr4.set_data([],[])
    return  tr2, tr3, tr4, SOL, TERRA, MART, JUPITER
def animate(i):
    i = int(i)
    tr2.set_data(terra_T[0][:i], [terra_T[1][:i]])
    tr3.set_data(mart_T[0][:i], [mart_T[1][:i]])
    tr4.set_data(jupiter_T[0][:i], jupiter_T[1][:i])
    SOL.set_data([0,sol_T[0][i]],[0, sol_T[1][i]])
    TERRA.set_data([0,terra_T[0][i]],[0, terra_T[1][i]])
    MART.set_data([0,mart_T[0][i]],[0, mart_T[1][i]])
    JUPITER.set_data([0,jupiter_T[0][i]],[0, jupiter_T[1][i]])
    return   tr2, tr3, tr4,SOL, TERRA, MART, JUPITER

plt.axis('off')
fig.patch.set_facecolor('black') 
index_list = np.linspace(0,782)
ani = animation.FuncAnimation(fig, animate, frames = len(sol_T[0]),interval = 10, blit = True, init_func = init)
dpi = 200
writer = animation.writers['ffmpeg'](fps=15)
ani.save('video_planetes.mp4',writer=writer,dpi=dpi)
