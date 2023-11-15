import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d



file = open(r'sistema_solar_posicions_postopt.txt','r')
Lines = file.readlines()

list_headers = Lines[0].split(' ')
print(Lines[0])

sol         =np.array([[float(line.split(' ')[list_headers.index('x_Sol')]), float(line.split(' ')[list_headers.index('y_Sol')]), float(line.split(' ')[list_headers.index('z_Sol')])] for line in Lines[1:]])
terra       =np.array([[float(line.split(' ')[list_headers.index('x_Terra')]), float(line.split(' ')[list_headers.index('y_Terra')]), float(line.split(' ')[list_headers.index('z_Terra')])] for line in Lines[1:]])
mart        =np.array([[float(line.split(' ')[list_headers.index('x_Mart')]), float(line.split(' ')[list_headers.index('y_Mart')]), float(line.split(' ')[list_headers.index('z_Mart')])] for line in Lines[1:]])
jupiter     =np.array([[float(line.split(' ')[list_headers.index('x_Jupiter')]), float(line.split(' ')[list_headers.index('y_Jupiter')]), float(line.split(' ')[list_headers.index('z_Jupiter')])] for line in Lines[1:]])
sol_T     = np.transpose(sol)
terra_T   = np.transpose(terra)
mart_T    = np.transpose(mart)
jupiter_T = np.transpose(jupiter)

ax = plt.figure().add_subplot(projection='3d')

ax.plot(terra_T[0], terra_T[1], terra_T[2], label = 'Terra')
ax.plot(sol_T[0], sol_T[1],sol_T[2], label = 'Sol')
ax.plot(mart_T[0], mart_T[1],mart_T[2], label ='Mart')
ax.plot(jupiter_T[0], jupiter_T[1],jupiter_T[2], label = 'Jupiter')
ax.set_title(f'Trajectòries Sistema solar')
ax.legend()
plt.show()