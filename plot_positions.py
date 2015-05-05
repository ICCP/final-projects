import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import math



f = open("data/raw_positions.txt")
numSats = 5

n = 0
data = [[]]

n2 = 0
for line in f:
	if line != '-\n':
		data_aux = []
		for x in line.split(","):
			data_aux.append(float(x)/(10.0**12.0))
		data[-1].append(data_aux)

	else:
		n+=1
		data_aux = []
		data.append([])


colors = ['#FFFF00','#663300','#99CC00','#0000FF','#FF0000','#FFCC33','#CC99FF','#006633','#000066','#000000']
sizes = [10**7,1.66012,24.4781,30.03467,3.22713,9545.8,2858.12,436.576,515.028,.066082]

colors += (['#FFFFFF']*numSats)
sizes += ([10**7]*numSats)

for i in range(0,len(sizes)):
	sizes[i]=sizes[i]**(2/3)
	#artificial scaling for viewability
	sizes[i]=sizes[i]**.3

for i in range(0,n):
	if i%1 == 0:
		print ('saving figure ',i)

		#fig = plt.figure()
		#ax = fig.add_subplot(111, projection='3d')
		# for j in range(0,len(data[i])):
		# 		ax.scatter(data[i][j][0], data[i][j][1], data[i][j][2],c=colors[j],s=sizes[j])
		# ax.set_xlabel('X (billion kilometers)')
		# ax.set_ylabel('Y (billion kilometers)')
		# ax.set_zlabel('Z (billion kilometers)')
		# plt.xlim([-8,8])
		# plt.ylim([-8,8])
		# ax.set_zlim([-8,8])
		# plt.title('Position at timestep '+str(i))
		# plt.savefig('images/3D/PositionData'+str(i)+'.png', bbox_inches='tight')
		# plt.close()

		fig = plt.figure()
		ax = fig.add_subplot(111)
		for j in range(0,len(data[i])):
				ax.scatter(data[i][j][0], data[i][j][1],c=colors[j],s=sizes[j])
		ax.set_xlabel('X (billion kilometers)')
		ax.set_ylabel('Y (billion kilometers)')
		plt.xlim([-.3,.3])
		plt.ylim([-.3,.3])
		plt.title('Position at timestep '+str(i))
		plt.savefig('images/TopDown/TopDown'+str(i)+'.png', bbox_inches='tight')
		plt.close()

		# fig = plt.figure()
		# ax = fig.add_subplot(111)
		# for j in range(0,len(data[i])):
		# 		ax.scatter(data[i][j][0], data[i][j][2],c=colors[j],s=sizes[j])
		# ax.set_xlabel('X (billion kilometers)')
		# ax.set_ylabel('Z (billion kilometers)')
		# plt.xlim([-8,8])
		# plt.ylim([-8,8])
		# plt.title('Position at timestep '+str(i))
		# plt.savefig('images/Front/Front'+str(i)+'.png', bbox_inches='tight')
		# plt.close()


