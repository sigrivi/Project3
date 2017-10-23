import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import linalg
import sys
import time

pi = np.pi
days = 365.25 			# number of days in a year

class Planet:
	def __init__(self, mass, r, v, name ):
		self.mass = mass
		self.r=np.asarray(r, dtype=np.float)
		self.v=np.asarray(v,dtype=np.float)
		self.a_s = 4*pi**2/1**2 #radius = 1
		self.a = np.zeros(len(r))
		self.name = name
	def total_energy(self):
		
		self.energy = 0.5 *self.mass *np.linalg.norm(self.v)**2 - 4*pi**2 *self.mass / np.linalg.norm(self.r)

		return(self.energy)
	def angular_momentum(self):
		self.l = np.linalg.norm(np.cross(self.r,self.v*self.mass))
		return(self.l)




class Acceleration:
	def __init__(self, Planets):
		self.N_p=len(Planets) ## Number of planets in the Planets list

	def Update(self, Planets, beta = 3.):
		for k in range(3): ## 3 is the length of the r and v vector

			for i in range(self.N_p):
				P=Planets[i]
				P.a[k] = 0
				c = 2*365.35*6*6*.1 #speed of light in AU/yr
				for j in range(self.N_p):
					if j==i:
						if P.name == "Mercury":
							l = np.linalg.norm(np.cross(P.r,P.v)) #angular momentum
							P.a[k] = P.a[k] -(4*pi**2*P.r[k]/np.linalg.norm(P.r)**beta) * (1+3*(l/(np.linalg.norm(P.r)*c))**2) ### relativistic correction
						else: 
							P.a[k] = P.a[k] -4*pi**2*P.r[k]/np.linalg.norm(P.r)**beta
						
					else: 
						distance = P.r-Planets[j].r
						#print(j, distance, Planets[j].r)
						P.a[k] = P.a[k] - distance[k]*Planets[j].mass/np.linalg.norm(distance)**beta
		return(Planets)



class DiffEqSolver:
	def __init__(self, h):
		self.h = h

	def EulerCromer(self, Planets):
		h = self.h
		for P in Planets:
			for i in range(len(P.v)):
				P.v[i] = P.v[i]+h*P.a[i]
				P.r[i] = P.r[i]+h*P.v[i]
		return Planets

	def VelocityVerlet(self,Planets, beta =3.):
		h = self.h
		for P in Planets:
			for i in range(len(P.v)):
				P.v[i] = P.v[i]+0.5*h*P.a[i]
				P.r[i] = P.r[i]+h*P.v[i]
		accel = Acceleration(Planets)
		accel.Update(Planets , beta)
		for P in Planets:
			for i in range(len(P.v)):
				P.v[i] = P.v[i]+0.5*h*P.a[i]




# solves the coupled differential equations. Planets is a list of planets, method is either "EulerCromer" or "VelocityVerlet". h is the step length. 
def solve(Planets, method, h, t_final =1., beta = 3.):
	solver = DiffEqSolver(h)
	accel = Acceleration(Planets)
	N_p = len(Planets)		# number of planets 
	N = np.int(t_final/h) 	# number of grid points
	x = np.zeros((N_p, N))	# array for x positions
	y = np.zeros((N_p, N))	# array for y positions


	# Calulate an position: 
	for j in range(N_p):
		x[j,0] = Planets[j].r[0]
		y[j,0] = Planets[j].r[1]
		
		if method == "EulerCromer":
			for i in range(1,N):
				accel.Update(Planets)
				solver.EulerCromer(Planets)
				x[j,i] = Planets[j].r[0]
				y[j,i] = Planets[j].r[1]
				#print("EulerCromer",np.linalg.norm(Planets[j].r[0]),np.linalg.norm(Planets[j].v[0]),np.linalg.norm(np.cross(Planets[j].r,Planets[j].v)))
		if method == "VelocityVerlet":
			for i in range(1,N):
				accel.Update(Planets,beta)
				solver.VelocityVerlet(Planets,beta)
				x[j,i] = Planets[j].r[0]
				y[j,i] = Planets[j].r[1]
				#print(Planets[j].l)
				#print("VelocityVerlet",np.linalg.norm(Planets[j].r),np.linalg.norm(Planets[j].v),np.linalg.norm(np.cross(Planets[j].r,Planets[j].v)))
	#print(Planets[0].r, Planets[0].v)
	return(x,y)




# A list only containing Earth orbiting the sun. Initial position x=1, y=0, Intital speed : vy = 2*pi
OnePlanet = [Planet(3e-6, [1,0,0], [0,2*pi,0], "Earth")]


#P_Newtonian = [MercuryNewtonian]
#x_Newtonian, y_Newtonian = solve(P_Newtonian, "VelocityVerlet",h)
#plt.plot(x_Newtonian[0], y_Newtonian[0])
#P_Relativistic = [Mercury]
#x_Relatvistic, y_Relativistic = solve(P_Relativistic, "VelocityVerlet",h)
#plt.plot(x_Relatvistic[0], y_Relativistic[0])
#plt.show()



def ComapareSolvers(h, t_final = 1.):
	P_E = [Planet(3e-6, [1,0,0], [0,2*pi,0], "Earth")]
	P_V = [Planet(3e-6, [1,0,0], [0,2*pi,0], "Earth")]
	initial_en_E = P_E[0].total_energy() # initial energy P_E
	initial_en_V = P_V[0].total_energy() # initial energy P_V
	initial_l_E = P_E[0].angular_momentum() # initial angular momentum P_E
	initial_l_V = P_V[0].angular_momentum() # initial angular momentum P_E
	l = np.linalg.norm(np.cross(P_E[0].r, P_E[0].v*P_E[0].mass))
	#print("before", initial_l_E/P_V[0].mass)

	x_verlet,y_verlet = solve(P_V, "VelocityVerlet", h, t_final)	# solve with Velocity Verlet

	x_euler,y_euler = solve(P_E, "EulerCromer", h, t_final) 		# solve with Euler Cromer

	diff_en_E = (initial_en_E - P_E[0].total_energy())/initial_en_E 	#relative difference in total energy Euler Cormer
	diff_en_V = (initial_en_V -  P_V[0].total_energy())/initial_en_V	#relatevi difference in total energy Velocity Verlet
	diff_l_E = (initial_l_E - P_E[0].angular_momentum())/initial_l_E	#relative difference in angular momentum Euler Cromer
	diff_l_V = (initial_l_V - P_V[0].angular_momentum())/initial_l_V	#relative difference in angular momentum Velocity Verlet
	#print("after",P_E[0].angular_momentum()/P_V[0].mass)
	l = np.linalg.norm(np.cross(P_E[0].r, P_E[0].v*P_E[0].mass))
	#print(P_E[0].r, P_E[0].v)
	
	plt.plot(x_euler[0],y_euler[0],label="EulerCromer")
	plt.plot(x_verlet[0],y_verlet[0],label="VelocityVerlet")
	plt.title("Solutions of the Sun-Earh system, h=%.2f" %h)
	plt.xlabel("x/AU")
	plt.ylabel("y/AU")
	plt.axis('equal')
	plt.legend()
	plt.show()
	plt.savefig('VV_vs_EC2',dpi=225)
	return(diff_en_E, diff_en_V,diff_l_E, diff_l_V)

(a,b,c,d) = ComapareSolvers(1/20, 1.)

H=10
t_final = 1
diff_en_E = np.zeros(H)
diff_en_V= np.zeros(H)
diff_l_E= np.zeros(H)
diff_l_V= np.zeros(H)
axis = np.zeros(H)
for i in range(1,H):
	print(i)
	h = 1/(1*i)
	N = np.int(t_final/h)
	axis[i-1] = h
	#diff_en_E[i-1], diff_en_V[i-1],diff_l_E[i-1], diff_l_V[i-1] = ComapareSolvers(h)
#plt.plot(np.log10(axis), np.log10(diff_en_E), np.log10(axis), np.log10(diff_en_V))
#plt.figure()
#plt.plot(axis)
#plt.savefig("totalenergy2.png")
#plt.xlabel("h")
#plt.ylabel("relative difference in energy")
#plt.show()
#print(axis.shape)
#plt.plot(axis, diff_l_E, axis, diff_l_V)
#plt.xlabel("h")
#plt.ylabel("relative difference in angular_momentum")
#plt.show()


def plot_difference(H): #does not worK!!!

	diff_en_E = np.zeros(10)
	diff_en_V= np.zeros(10)
	diff_l_E= np.zeros(10)
	diff_l_V= np.zeros(10)
	axis = np.zeros(10)
	for i in range(1,10):
		hh = 1/(10*i)
		N = np.int(t_final/hh)
		axis[i-1] = hh
		diff_en_E[i-1], diff_en_V[i-1],diff_l_E[i-1], diff_l_V[i-1] = ComapareSolvers(hh)
	plt.plot(axis, diff_en_E, axis, diff_en_V)
	#plt.xlabel("h")
	#plt.ylabel("relative difference in energy")
	plt.show()
	#plt.plot(axis, diff_l_E, axis, diff_l_V)
	#plt.xlabel("h")
	#plt.ylabel("relative difference in angular_momentum")
	#plt.show()


#plot_difference(10)
#plt.show()
