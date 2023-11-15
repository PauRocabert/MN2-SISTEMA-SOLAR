import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

file = open(r'trayectorias.txt', 'r')
Lines = file.readlines()

k=len(Lines[0].split(" "))-2

dt = float(Lines[2].split(' ')[0]) - float(Lines[1].split(' ')[0])
t = np.array([float(line.split(' ')[0])]for line in Lines[1:])
labels=["Sol", "Terra", "Mart", "Jupiter", "Mercuri","Venus","Saturn"]

ax = plt.figure().add_subplot(projection='3d')

for i in range(int(k/4)):
    print(i)
    p =np.array([[float(line.split(' ')[4*i+1]), float(line.split(' ')[4*i+2]), float(line.split(' ')[4*i+3])] for line in Lines[1:]])
    p_T = (np.transpose(p))
    ax.plot((p_T[0]), (p_T[1]), (p_T[2]), label = labels[i])


plt.title(f'Trajectòries Sistema solor 13/9/2021 - 3/11/2023 RG-4 dt ={round(dt*365*24*60/(2*np.pi),2)} MIN')
ax.legend()


plt.show() 