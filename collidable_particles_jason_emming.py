#Collidable Particles
#Programmed by Jason Emming

from __future__ import division
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import random
import time
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### INITIAL  CONDITIONS ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

global G,N,T,L,h,skip_factor

#Setup determines the initialized conditions for the particles
options = ["random", "sun", "test1", "test2", "test3", "test4", "test5"]
setup = "definitly not a setup"

try:
	while options.count(setup.lower()) < 1: 
		print "Options: random, sun, test1, test2, test3, test4, test5"
		setup = raw_input("Initialize the setup: ")
		print '\n'

except:
	print "You done goofed, the setup of 'test3' has been chosen for you"
	setup = 'test3' 

G = 1 						#Gravitational Constant

T = 10000					#Temperature of particles
L = 6						#Radius of initialized distribution

f = 800						#Number of times smaller each particle's radius is than the init radius
p_size = L/f 				#Sets the radius of each particle as some fraction of init radius

h_min = 0.00001				#Lowest timestep needed (so far)
h=h_min						#Timestep

top_v = 10000		#Highest velocity ever recorded
 	
standard_time = 1 	#A useful unit to vary simulation time

skip_factor = h_min**(-1.0)//600			#Determines how many frames are played back per second on animation
max_frames = standard_time*h_min**(-1.0)	#Runs animation for 


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ###   PARTICLE  CLASS   ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


class Particle(object):

	def __init__(self, pos = [0,0,0], mom=np.array([0,0,0]), size=0, m=150.0, collis=False ,color='ro'):
		self.position = pos
		self.momenta = mom
		self.radius = size
		self.mass = m
		self.collis = collis
		self.color = color

	"""Accepts particle's position; returns net acceleration"""
	def acceleration(self, position1):
		ax, ay, az = 0, 0, 0
		for p2 in particles:
			if p2 != self:
				r = distance(self,p2)
				a = G*p2.mass / r**2.0
				a_x = a * (p2.position[0] - position1[0])/r
				a_y = a * (p2.position[1] - position1[1])/r
				a_z = a * (p2.position[2] - position1[2])/r
				ax += a_x
				ay += a_y
				az += a_z
		return (ax, ay, az)		

	""" Calculates the Kinetic energy of the particle """
	def kinetic(self):
		ke = 0.5*self.mass*np.sum(self.momenta**2.0)
		return ke

	""" Calculates the Potential energy between two particles """
	def potential(self,particle2):
		r = distance(self,particle2)
		if r == 0.0: 
			r = 0.00000000000000000001
		pe = -G*self.mass*particle2.mass/r
		return pe

	""" Adds the particle to the other particle (argument) """
	def add(self, particle2):
		#pos = (self.position + particle2.position)/2
		pos = (self.mass*self.position+particle2.mass*particle2.position)/(self.mass+particle2.mass)
		momenta = (self.momenta*self.mass + particle2.momenta*particle2.mass)/(self.mass + particle2.mass)
		size = (self.radius**3 + particle2.radius**3)**(1/3.0)
		m = self.mass + particle2.mass
		return Particle(pos,momenta,size,m)

	""" Prints the properties of the particle """
	def print_particle(self):
		print "{:<10} {:>40}".format("Position", round(self.position[0],3), round(self.position[1],3), round(self.position[2],3))
		print "{:<10} {:>40}".format("Momentum", round(self.momenta[0],3), round(self.momenta[1],3), round(self.momenta[2],3))
		print "{:<10} {:>40}".format("Radius",self.radius)
		print "{:<10} {:>40}".format("Mass", self.mass)
		print "{:<10} {:>40}".format("Kinetic", self.kinetic())

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ###  FUNCTIONS  ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


#Returns x randomly generated numbers from 0 up to b as a list
def randomish(x,b):
	rand_list = []
	for i in range(0,x):
		rand_list.append(random.uniform(0,b))
	return rand_list


#Returns a gaussian distribution of momenta as a list for N particles at T temperature
def momentum(N,T):
	a=0						#Value the distribution centers around
	b=(2*T)**(1.0/2)		#Standard deviation of momentums. Same spread used from MD Project.

	momenta = np.random.normal(a,b,(N,3))
	avg = np.zeros((3), dtype=float)

	for i in range(0, 3):
		avg[i] = sum(momenta[:,i]/N)
		momenta[:,i] = [j - avg[i] for j in momenta[:,i]]

	return momenta


#Returns a randomized list of positions inside a sphere of radius L centered at the orgin
def position(N):
	position_list = []

	while len(position_list) != N:
		rand = randomish(3,L)
		if ((rand[0]**2 + rand[1]**2 + rand[2]**2) <= L**2.0):
			r_int = [random.randint(1,2),random.randint(1,2),random.randint(1,2)]
			for i in range(len(rand)):
				rand[i] = rand[i]*(-1.0)**r_int[i]
			position_list.append(rand)

	return position_list


#Returns an array of N particles with relevant properties initialized
def initialize(N,T,L):
	momenta_list = momentum(N,T)
	position_list = position(N)

	particles = np.zeros(N,dtype=Particle)

	for i in range(len(particles)):
		p = Particle()

		p.position = position_list[i]
		p.momenta = momenta_list[i]/p.mass
		p.radius = p_size

		particles[i] = p

	return particles


#Calculates and returns the distance between two particles
def distance(particle1, particle2):
	pos1, pos2 = particle1.position, particle2.position
	distance = ((pos1[0]-pos2[0])**2 + (pos1[1]-pos2[1])**2 + (pos1[2]-pos2[2])**2)**(1/2.0)
	return distance


#If two particles are close enough they will collide and combine their mass and velocity vectors
def collision(particle1,particle2):
	d = distance(particle1, particle2)
	if d <= (particle1.radius + particle2.radius):
		particle1.collis, particle2.collis = True, True
		new_particle = particle1.add(particle2)
		print "COLLISION!"
		return new_particle


#Adds 3D vectors together
def vec_add(v1,v2):
	v = np.zeros(3)
	for i in range(3):
		v[i] = v1[i] + v2[i]
	return v


#Multiplies a 3D vector, v, by a scalar, s
def vec_scalar(v,s):
	v = np.array(v)
	for i in range(3):
		v[i] = v[i]*s
	return v


#Initializes the animation plot
def init():
	path = ax.scatter3D([],[],[])
	return path


#Updates the position of the particles for the animation
def update_position(i):
	ax.clear()

	ax.set_xlim3d([-L, L])
	ax.set_xlabel('X')

	ax.set_ylim3d([-L, L])
	ax.set_ylabel('Y')

	ax.set_zlim3d([-L, L])
	ax.set_zlabel('Z')

	x = data[0][i*skip_factor]
	y = data[1][i*skip_factor]
	z = data[2][i*skip_factor]
	path = ax.scatter3D(x,y,z)

	return path


#Runs Runge-Kutta Method on a single particle
def rk4(p):
	position = p.position
	
	#Calc kv1
	a = p.acceleration(position)
	kv1 = a

	#Calc kr1 
	kr1 = p.momenta

	#Calc kv2
	new_position = vec_add(position, vec_scalar(kr1,h/2.0))
	a = p.acceleration(new_position)
	kv2 = a

	#Calc kr2
	kr2 = vec_add(p.momenta, vec_scalar(kv2, h/2.0))

	#Calc kv3
	new_position = vec_add(position, vec_scalar(kr2,h/2.0))
	a = p.acceleration(new_position)
	kv3 = a

	#Calc kr3
	kr3 = vec_add(p.momenta, vec_scalar(kv3, h/2.0))

	#Calc kv4
	new_position = vec_add(position, vec_scalar(kr3, h))
	a = p.acceleration(new_position)
	kv4 = a

	#Calc kr4
	kr4 = vec_add(p.momenta, vec_scalar(kv4, h))

	#Update final position
	p.position = vec_add(position, phi(kr1,kr2,kr3,kr4))

	#Update final velocity
	p.momenta = vec_add(p.momenta, phi(kv1,kv2,kv3,kv4))	

	return p.position, p.momenta


#Part of the Runge-Kutta Method
def phi(k1,k2,k3,k4):
	k_sum = vec_add( vec_add( k1, vec_scalar(k2, 2.0) ), vec_add(vec_scalar(k3,2.0),k4 ) )
	phi = vec_scalar(k_sum, h/6.0)	

	return phi


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ###   MAIN  ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


#Runs the program and generates the animation
def main():

	#Initalizes global variables, which are used throughout the code
	global init_particles, particles, ax, data, new_particles

	#Initializes a system of particles based on 'setup' condition
	if setup == "random":
		N = 25
		init_particles = initialize(N,T,L)

	#Initializes a solar-like system with 5 planets (Mercury to Jupiter)
	elif setup == "sun":
		N = 6
		init_particles = initialize(N,T,L)
		init_particles[0].mass = 333000
		init_particles[0].momenta = np.array([0,0,0])
		init_particles[0].position = np.array([0,0,0])
		init_particles[1].mass = 0.0553
		init_particles[1].momenta = np.array([0,1000,0])
		init_particles[1].position = np.array([0.3529,0,0])
		init_particles[2].mass = 0.8153
		init_particles[2].momenta = np.array([-750,0,0])
		init_particles[2].position = np.array([0,0.719,0])		
		init_particles[3].mass = 1
		init_particles[3].momenta = np.array([0,666,0])
		init_particles[3].position = np.array([1,0,0])
		init_particles[4].mass = 0.10748
		init_particles[4].momenta = np.array([-500,0,0])
		init_particles[4].position = np.array([0,1.5,0])
		init_particles[5].mass = 318
		init_particles[5].momenta = np.array([0,276,0])
		init_particles[5].position = np.array([5.3,0,0])

	#Initializes a central mass with two orbiting bodies: one massive, other small
	elif setup == 'test1':
		N = 3
		init_particles = initialize(N,T,L)
		init_particles[0].mass = 1000
		init_particles[0].momenta = np.array([0,0,0])
		init_particles[0].position = np.array([0,0,0])
		init_particles[0].print_particle()
		print '\n'
		init_particles[1].mass = 5
		init_particles[1].momenta = np.array([0,30,0])
		init_particles[1].position = np.array([1,0,0])
		init_particles[1].print_particle()
		print '\n'
		init_particles[2].mass = 25
		init_particles[2].momenta = np.array([0,10,0])
		init_particles[2].position = np.array([4,0,0])
		init_particles[2].print_particle()

	#Initializes 4 equal mass objects, orbiting a common center of mass
	elif setup == 'test2':
		N = 4
		init_particles = initialize(N,T,L)
		init_particles[0].mass = 700
		init_particles[0].momenta = np.array([0,-8,0])
		init_particles[0].position = np.array([5,0,0])
		init_particles[0].print_particle()
		print '\n'
		init_particles[1].mass = 700
		init_particles[1].momenta = np.array([0,8,0])
		init_particles[1].position = np.array([-5,0,0])
		init_particles[1].print_particle()
		print '\n'
		init_particles[2].mass = 700
		init_particles[2].momenta = np.array([-8,0,0])
		init_particles[2].position = np.array([0,-5,0])
		init_particles[2].print_particle()
		print '\n'
		init_particles[3].mass = 700
		init_particles[3].momenta = np.array([8,0,0])
		init_particles[3].position = np.array([0,5,0])
		init_particles[3].print_particle()
	
	#Initializes a simple orbit: Central mass orbited by 1 mass mass
	elif setup == 'test3':
		N = 2
		init_particles = initialize(N,T,L)
		init_particles[0].mass = 1000
		init_particles[0].momenta = np.array([0,0,0])
		init_particles[0].position = np.array([0,0,0])
		init_particles[0].print_particle()
		print '\n'
		init_particles[1].mass = 1
		init_particles[1].momenta = np.array([0,30,0])
		init_particles[1].position = np.array([1,0,0])
		init_particles[1].print_particle()

	#Initializes N number of particles randomly distributed in the z=0 plane
	elif setup == 'test4':
		N = 50
		init_particles = initialize(N,T,L)
		init_particles[0].mass = 10000
		init_particles[0].radius = p_size*(init_particles[0].mass/300)**(1/3.0)
		init_particles[0].momenta = init_particles[0].momenta/init_particles[0].mass
		init_particles[0].position = np.array([0,0,0])
		init_particles[0].print_particle()
		for p in range(1,N):
			init_particles[p].mass = 100
			init_particles[p].momenta[2] = 0
			init_particles[p].momenta *= 80
			init_particles[p].position[2] = 0
			init_particles[p].print_particle()

	#Initializes 4 equal mass objects which interact with one another
	elif setup == 'test5':
		N = 4
		init_particles = initialize(N,T,L)
		init_particles[0].mass = 500
		init_particles[0].momenta = np.array([0,-11,0])
		init_particles[0].position = np.array([4,0,0])

		init_particles[1].mass = 500
		init_particles[1].momenta = np.array([0,11,0])
		init_particles[1].position = np.array([5,0,0]) 

		init_particles[2].mass = 500
		init_particles[2].momenta = np.array([0,-11,0])
		init_particles[2].position = np.array([-5,0,0])

		init_particles[3].mass = 500
		init_particles[3].momenta = np.array([0,11,0])
		init_particles[3].position = np.array([-4,0,0])

	print "\nRunning setup =", setup,'\n'

	#Array with 3 elements: x,y,z coordinates. Each element contains a list.
	#Each list contains the x or y or z coordinates for a single particle through all timesteps
	data = np.zeros(3,dtype=list)
	x,y,z =np.zeros(max_frames+1,dtype=list),np.zeros(max_frames+1,dtype=list),np.zeros(max_frames+1,dtype=list)

	#Array with 2 elements: Kinetic & Potential Energy. Each element is itself an array.
	#Each array contains the K or P energy for the system for all timesteps.
	energy_array = np.zeros(2, dtype=list)

	#Holds the Kinetic/Potential data for all particles at current timestep
	ke_tot, pe_tot = np.zeros(max_frames+1), np.zeros(max_frames+1) 
	
	#List records the velocities of each timestep 
	#max_velocity holds the highest recorded value
	velocity_list, max_velocity = [], 0	

	#Holds all the 'new' particles created from collisions
	new_particles = []

	#Loop that runs through 'max_frames' number of timesteps
	count = 0
	while count < max_frames:

		#Takes the original initialized list and makes an offical copy of it called 'particles'
		particles = init_particles.tolist()
		
		if velocity_list != []:
			#Updates max_velocity
			if max(velocity_list) > max_velocity:
				max_velocity = max(velocity_list)
		#Clear velocity_list for next iteration 
		velocity_list = []

		#Print the number of particles in simulation every 1000 timesteps
		if count % 1000 == 0: 
			print "Number of Particles: ", len(particles)

		#Print the number of particles in simulation every 100 timesteps
		if count % 100 == 0:
			print "C: ", round((count/max_frames)*100,3), "%"
			
		#Holds the x,y,z positions of all particles at current timestep
		x_f, y_f, z_f = [], [], []

		ke, pe, = 0, 0

		#Checks for collided particles (p.collis = True) and removesthem from particle array
		for k in range(len(init_particles)):
			p = init_particles[k]
			if p.collis:
				particles.remove(p)

		#Adds all newly created particles to the particle array		
		for new in new_particles:
			particles.append(new)
		init_particles = np.array(particles)
		#Clears the newly created particles from 'new_particles' for the next timestep
		new_particles = []

		#This list temporarily stores the new position and 'momenta'=velocity till RK4 is complete
		temp_pm = []

		#This loop runs through each particle and determines their new position and velocity
		i = 0
		for p in particles:

			#Stores the particles x,y,z position
			x_f.append(p.position[0])
			y_f.append(p.position[1])
			z_f.append(p.position[2])	

			#Calculates the sum Kinetic energy
			ke += p.kinetic()

			#Calculates the sum of the Potential energy
			i+=1
			for p2 in particles[i:]:
				pe += p.potential(p2)

			#Checks if a two particles are close enough to register a collision
			for p2 in particles:
				if p != p2:
					a = collision(p,p2)
					if a != None:
						new_particles.append(a)

			#Calculates the position and velocity of each particle
			temp_position, temp_momenta = rk4(p)
			#Adds pos/vel to a temporary list until all particles has been addressed
			temp_pm.append([temp_position, temp_momenta])

		#Updates the position and velocity of every particle from temporary list
		for j in range(len(particles)):
			p = particles[j]
			p.position = temp_pm[j][0]
			p.momenta  = temp_pm[j][1]
			if p.collis == False:
				velocity_list.append(np.sqrt(np.sum(p.momenta**2.0)))

		#Records the total Kintetic & Potential Energy for every particle for every timestep
		ke_tot[count] = ke
		pe_tot[count] = pe

		#Records the x,y,z coordinates of every particle for every timestep			
		x[count] = x_f
		y[count] = y_f
		z[count] = z_f			
		
		#End of loop at 1 to count and repeat until count >= max_frames
		count += 1

	#Puts all the Energy data into energy array
	energy_array[0] = ke_tot
	energy_array[1] = pe_tot

	#Wrapes all position coordinates into data array
	data[0] = x
	data[1] = y
	data[2] = z

	print "Highest Velocity Achieved in Run: ", max_velocity
	if max_velocity > top_v:
		print "HEY YOU'VE RECORDED A NEW HIGH VELOCITY: ", max_velocity

	#Writes the positions (xyz) of all the particles for every timestep to a text file
	fo = open("data_lunar.txt", "wb")
	for a in data:
		fo.write('\tMarker\t')
		for b in a:
			fo.write(str(b)+'\t')
	fo.close()
		
	#Plots a graph of Kinetic/Potential/Total Energy of system over all timesteps
	ke_list = np.array(energy_array[0])
	pe_list = np.array(energy_array[1])
	tot_energy = ke_list + pe_list
	steps = range(0,int(max_frames+1))

	#Prints out the percent of energy conserved over the course of the simulation
	print round(tot_energy[-1]/tot_energy[0] * 100,2) ,'%'' Energy Conserved'
	
	#Plots the Kinetic/Potential/Total Energy of the simulation
	plt.xlabel("Timesteps")
	plt.ylabel('Energy')
	title = "Timestep is "+ str(h)
	plt.title(title)
	plt.plot(steps, ke_list, 'r', steps, pe_list,'g', steps, tot_energy, 'k')
	plt.show()

	#Runs an animation of the particle's positions
	fig = plt.figure()
	ax = p3.Axes3D(fig)	
	ax.set_axis_off()
	ax.axis('off') 
	numframes = int(max_frames/skip_factor)
	anim = animation.FuncAnimation(fig, update_position, numframes, blit=False, init_func=init)
	plt.show()


main()