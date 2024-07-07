import numpy as np
import matplotlib.pyplot as plt


file = open(r'sistema_solar_posicions_postopt.txt','r')
Lines = file.readlines()

list_headers = Lines[0].split(' ')
print(Lines[0])

sol         =np.array([[float(line.split(' ')[list_headers.index('x_Sol')]), float(line.split(' ')[list_headers.index('y_Sol')]), float(line.split(' ')[list_headers.index('E_Sol')])] for line in Lines[1:]])
terra       =np.array([[float(line.split(' ')[list_headers.index('x_Terra')]), float(line.split(' ')[list_headers.index('y_Terra')]), float(line.split(' ')[list_headers.index('E_Terra')])] for line in Lines[1:]])
mart        =np.array([[float(line.split(' ')[list_headers.index('x_Mart')]), float(line.split(' ')[list_headers.index('y_Mart')]), float(line.split(' ')[list_headers.index('E_Mart')])] for line in Lines[1:]])
jupiter     =np.array([[float(line.split(' ')[list_headers.index('x_Jupiter')]), float(line.split(' ')[list_headers.index('y_Jupiter')]), float(line.split(' ')[list_headers.index('E_Jupiter')])] for line in Lines[1:]])
sol_T     = np.transpose(sol)
terra_T   = np.transpose(terra)
mart_T    = np.transpose(mart)
jupiter_T = np.transpose(jupiter)


fig, axes = plt.subplots(2)
axes[0].plot(terra_T[0], terra_T[1], label = 'Terra')
axes[0].plot(sol_T[0], sol_T[1], label = 'Sol')
axes[0].plot(mart_T[0], mart_T[1], label ='Mart')
axes[0].plot(jupiter_T[0], jupiter_T[1], label = 'Jupiter')
axes[0].set_title(f'Trajectòries Sistema solor 13/9/2021 - 3/11/2023 RK-4 dt ={1} min')
axes[0].legend()
axes[1].plot(sol_T[2], label ='E_sol')
axes[1].plot(terra_T[2], label ='E_Terra')
axes[1].plot(mart_T[2], label='E_mart')
axes[1].plot(jupiter_T[2], label='E_jupiter')
axes[1].legend()
plt.show() 

