import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import math



f = open("data/delta_E.txt")
n=0
data=[[],[],[]]
for line in f:
	data_type=1
	for num in line.split(", "):
		data[data_type].append(float(num))
		data_type+=1
	data[0].append(n)
	n+=1

totalEn = [sum(x) for x in zip(data[1], data[2])]
index = 0
colors=['r','k']
labels=['Kinetic Energy','Potential Energy']

plt.plot(data[0],data[1],color=colors[index],label=labels[index])
index += 1
plt.plot(data[0],data[2],color=colors[index],label=labels[index])
plt.xlabel('Simulation Time (10^4 Seconds)')
plt.ylabel('Energy (Joules)')
plt.title('Energy Vs Simulation Time')
#plt.ylim([-0.5*10**36.0,0.5*10**36.0])
plt.plot(data[0],totalEn,color='k',label='Total Energy')

plt.legend(loc='right')
plt.savefig('images/Plots/EnergyvsTimesteps.png', bbox_inches='tight')
plt.show()
plt.close()
