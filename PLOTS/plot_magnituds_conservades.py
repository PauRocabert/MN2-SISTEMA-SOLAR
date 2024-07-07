import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import datetime 

plt.style.use('seaborn-darkgrid')

file = open(r'trajectories_MAGNITUDS_4C_terramig.txt','r')
Lines = file.readlines()
list_headers = Lines[0].split(' ')

N_cossos =  int((len(list_headers)-8)/5)
N_dim = 3
N_dies = len(Lines)-1

labels = [word.split('_')[1] for i,word in enumerate(list_headers[1:-7]) if i%5 == 0]
print(labels)

fEnergies = open('Error_hora_terra2.txt','w')
for label in labels: fEnergies.write(' '+label)

DADES = np.ndarray(shape=(N_cossos,N_dim, N_dies))
ERROR_PLANETES = np.ndarray(shape=(N_cossos, N_dies))


fig = plt.figure()
ax = fig.add_subplot(projection="3d")
time = np.array([float(line.split(' ')[0])for line in Lines[1:]])
fEnergies.write('\nmaxErr')
for cos in range(N_cossos):
    ERROR_PLANETES[cos] = np.array([float(line.split(' ')[5*(cos+1)]) for line in Lines[1:]])
    fEnergies.write(' '+str(max(ERROR_PLANETES[cos])))
    DADES[cos] = np.transpose(np.array([[float(line.split(' ')[5*cos+1]), float(line.split(' ')[5*cos+2]),float(line.split(' ')[5*cos+3])]for line in Lines[1:]]))
    ax.plot(DADES[cos][0],DADES[cos][1], DADES[cos][2], label = labels[cos])

plt.title(f'Trajectòries sistema solar 13/9/2021 - 3/11/2023  dt ={1} MIN')
ax.set_xlabel('x (UA)')
ax.set_ylabel('y (UA)')
ax.set_zlabel(' z (UA)')
ax.legend()
plt.show() 


fig = plt.figure()
ax = fig.add_subplot()
E0 = np.array([float(Lines[1].split(' ')[5*cos+4]) for cos in range(N_cossos)])
fEnergies.write('\nmaxErrRelatiu')
for cos in range(N_cossos):
    if E0[cos] !=0:
        ax.plot(time, ERROR_PLANETES[cos]/E0[cos], label = labels[cos])
        fEnergies.write(' '+str(max(ERROR_PLANETES[cos]/E0[cos])))
    else:
        fEnergies.write(' '+str(0.0))
    

plt.legend()
ax.set_title(r"Diferència de l'energia respecte l'energia inicial (dt = 1 min)", fontweight= 'book')
ax.title.set_size(15)
ax.set_xlabel('temps (dies)')
ax.set_ylabel(r'$\frac{\Delta E}{E_0}$')
plt.show()



E0_total = float(Lines[1].split(' ')[-1][:-1])
ERRORS = np.array([float(line.split(' ')[-1]) - E0_total for line in Lines[1:]])/E0_total
fEnergies.write('\nErrorTotalRelatiuEnergia '+str(max(abs(ERRORS))))
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(time, ERRORS)
ax.set_title(r"Diferència de l'energia total respecte l'energia total inicial (dt = 1 min)", fontweight= 'book')
ax.title.set_size(15)
ax.set_xlabel('temps (dies)')
ax.set_ylabel(r'$\frac{\Delta E}{E_0}$')
plt.show()


fig = plt.figure()
ax = fig.add_subplot()
Lz0 =float(Lines[1].split(' ')[-2])  
Lz = np.array([float(line.split(' ')[-2])  for line in Lines[1:]])
Ly0 =float(Lines[1].split(' ')[-4])  
Ly = np.array([float(line.split(' ')[-4])  for line in Lines[1:]])
Lx0 =float(Lines[1].split(' ')[-6])  
Lx = np.array([float(line.split(' ')[-6])  for line in Lines[1:]])
L0 =np.sqrt(Lx0**2+Ly0**2+Lz0**2)
L = np.sqrt((Lx**2+Ly**2+Lz**2))
Error_L = abs(L-L0)/L0
fEnergies.write('\nErrorTotalRelatiuMoment angular '+str(max(abs(Error_L))))
ax.plot(time,Error_L)
ax.set_xlabel('temps (Dies)')
ax.set_ylabel(r'$\frac{\Delta L}{L}$')
ax.set_title(r"Diferència del mòdul del moment angular total respecte de l'inicial (dt = 1 dia)", fontweight= 'book')
ax.title.set_size(15)
plt.show()



fig = plt.figure()
ax = fig.add_subplot()
Pz0 =float(Lines[1].split(' ')[-3])  
Pz = np.array([float(line.split(' ')[-3])  for line in Lines[1:]])
Py0 =float(Lines[1].split(' ')[-5])  
Py = np.array([float(line.split(' ')[-5])  for line in Lines[1:]])
Px0 =float(Lines[1].split(' ')[-7])  
Px = np.array([float(line.split(' ')[-7])  for line in Lines[1:]])
P0 =np.sqrt(Px0**2+Py0**2+Pz0**2)
P = np.sqrt((Px**2+Py**2+Pz**2))
Error_P = abs(P-P0)/P0
fEnergies.write('\nErrorTotalRelatiuMoment lineal '+str(max(abs(Error_P))))
ax.plot(time,Error_P)
ax.set_xlabel('temps (Dies)')
ax.set_ylabel(r'$\frac{\Delta P}{P}$')
ax.set_title(r"Diferència del mòdul del moment lineal total respecte de l'inicial (dt = 1 dia)", fontweight= 'book')
ax.title.set_size(15)
plt.show()



