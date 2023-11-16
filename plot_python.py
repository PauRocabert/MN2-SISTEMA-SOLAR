import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import datetime 


file = open(r'trayectorias.txt','r')
Lines = file.readlines()


list_headers = Lines[0].split(' ')

N_cossos =  int((len(list_headers)-1)/4)
N_dim = 3
N_dies = len(Lines)-1

DADES = np.ndarray(shape=(N_cossos,N_dim, N_dies))
fig = plt.figure()
ax = fig.add_subplot(projection="3d")

for cos in range(N_cossos):
    DADES[cos] = np.transpose(np.array([[float(line.split(' ')[4*cos+1]), float(line.split(' ')[4*cos+2]),float(line.split(' ')[4*cos+3])]for line in Lines[1:]]))
    ax.plot(DADES[cos][0],DADES[cos][1], DADES[cos][2])
    
plt.title(f'Trajectòries Sistema solor 13/9/2021 - 3/11/2023 RG-4 dt ={1} MIN')
ax.legend()
plt.show() 
    



    
