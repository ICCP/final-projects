import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import math



f = open("data/delta_E.txt")

data=[[[],[],[],[],[]]]
int_method=0
time_step=0
for line in f:
	if line == "-\n":
		data.append([[],[],[],[],[]])
		int_method+=1
	else:
		data_type=0
		for num in line.split(", "):
			data[int_method][data_type].append(float(num))
			data_type+=1
		time_step+=1
data = data[:-1]
int_method = 0

index = 0
colors=['r','k']
labels=['Runga-Kutta-4','Euler']
for int_method in data:
	plt.plot(int_method[1],int_method[0],color=colors[index],label=labels[index])
	plt.xlabel('10,000s or timesteps')
	plt.ylabel('Fractional Energy Change')
	index+=1
plt.legend(loc='lower left')
plt.savefig('images/Plots/EnergyVSTimesteps.png', bbox_inches='tight')
plt.close()

index=0
for int_method in data:
	slope = (data[index][1][-1]-data[index][1][0])/(data[index][2][-1]-data[index][2][0])
	slope = str(round(1/slope,3))
	plt.plot(int_method[2],int_method[1],color=colors[index],label=(labels[index]+": ms per timestep="+slope))
	plt.xlabel('Realtime (ms)')
	plt.ylabel('Timesteps')
	index+=1
plt.legend(loc='lower left')
plt.savefig('images/Plots/TimestepsVSRealtime.png', bbox_inches='tight')
plt.close()

index = 0
for int_method in data:
	plt.plot(int_method[1],int_method[3],color=colors[index],label=labels[index])
	plt.xlabel('10,000s or timesteps')
	plt.ylabel('Kinetic Energy')
	index+=1
plt.legend(loc='lower left')
plt.savefig('images/Plots/KineticVSTimesteps.png', bbox_inches='tight')
plt.close()

index = 0
for int_method in data:
	plt.plot(int_method[1],int_method[4],color=colors[index],label=labels[index])
	plt.xlabel('10,000s or timesteps')
	plt.ylabel('Potential Energy')
	index+=1
plt.legend(loc='lower left')
plt.savefig('images/Plots/PotentialVSTimesteps.png', bbox_inches='tight')
plt.close()

#print(data[0][1])
#print(data[0][4])