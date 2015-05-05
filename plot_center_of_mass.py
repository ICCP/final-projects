import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import math



f = open("data/center_of_mass.txt")

data=[[],[],[],[]]
int_method=0
time_step=0
for line in f:
	data_type=1
	data[0].append(time_step)
	for num in line.split(","):
		data[data_type].append(float(num))
		data_type+=1
	time_step+=10

f.close()
index = 1


for index in range(1,4):
	if(index == 2):
		data[0]=data[0][:-1]
		#data[index]=data[index][:-1]
	print("len data[0]: " + str(len(data[0])))
	print("len data[index]: " + str(len(data[index])))
	print("data[0][0]: " + str(data[0][100]))
	print("data[index][0]: " + str(data[index][100]))
	plt.plot(data[0],data[index])
	plt.xlabel('10,000s or timesteps')
	plt.ylabel('Distance from center (m)')
	plt.savefig('images/Plots/CoMVSTime'+str(index)+'.png', bbox_inches='tight')
	plt.close()

