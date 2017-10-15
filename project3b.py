import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import linalg
import sys
import time

pi = np.pi

class Planet:
	def __init__(self, mass, r, v, name ):
		self.mass = mass
		self.r=np.asarray(r, dtype=np.float)
		self.v=np.asarray(v,dtype=np.float)
		self.a_s = 4*np.pi**2/1**2 #radius = 1
		self.a = np.zeros(len(r))
		self.name = name
		#self.yx = 0
		#self.forces = np.zeros(len(r))	



class Acceleration:
	def __init__(self, Planets):
		self.number_of_planets=len(Planets)

	def Update(self, Planets):
		for k in range(3): ## 3 is the length of the r and v vector

			for i in range(self.number_of_planets):
				P=Planets[i]
				P.a[k] = 0
				for j in range(self.number_of_planets):
					if j==i:
						P.a[k] = P.a[k] -4*pi**2*P.r[k]/np.linalg.norm(P.r)**3
						
					else: 
						distance = P.r-Planets[j].r
						#print(j, distance, Planets[j].r)
						P.a[k] = P.a[k] - distance[k]*Planets[j].mass/np.linalg.norm(distance)**3
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

	def VelocityVerlet(self,Planets):
		h = self.h
		for P in Planets:
			for i in range(len(P.v)):
				P.v[i] = P.v[i]+0.5*h*P.a[i]
				P.r[i] = P.r[i]+h*P.v[i]
		accel = Acceleration(Planets)
		accel.Update(Planets)
		for P in Planets:
			for i in range(len(P.v)):
				P.v[i] = P.v[i]+0.5*h*P.a[i]

N = 5000
t_stop = 100
h = 1/t_stop
days=365.25 #number of days in a year

# Initial values of the planets from 2017-Oct-04:
Mercury = Planet(1.65e-7, [-0.3789210,2.56223e-02,0],[-7.32704e-03*days, -2.68813e-02*days,0], "Mercury")
Venus = Planet(2.45e-6, [-.483111, .534117,0],[-1.49738e-02*days, -1.37846e-02*days,0 ], "Venus")
Earth = Planet(3e-6,[.985157,.191737e-01,0],[-3.48636e-3*days,1.68385e-02*days,0],"Earth")
Mars = Planet(3.3e-7, [-1.49970, 0.724618,0], [-5.52182e-03*days, -1.14214e-02*days,0],"Mars")
Jupiter = Planet(9.5e-4,[-4.63556,-2.8425,0],[3.85608e-03*days,-6.075643e-03*days,0],"Jupiter") 
Saturn = Planet(2.75e-4,[-0.421224,-10.0462,0],[5.26820e-03*days,-2.52167e-04*days,0], "Saturn")
Uranus = Planet(4.4e-5,[17.8825, 8.76632,0],[-1.75983e-03*days,3.34820e-03*days],"Uranus")
Neptune = Planet(0.515e-4, [28.6020, -8.86057, 0],[9.08471e-04*days, 3.01748e-03*days,0],"Neptune")


def ComapareSolvers(h):
	P_E = [Planet(3e-6, [1,0,0], [0,2*pi,0], "Earth")]
	P_V = [Planet(3e-6, [1,0,0], [0,2*pi,0], "Earth")]
	accel_E = Acceleration(P_E)
	accel_V = Acceleration(P_V)
	solver = DiffEqSolver(h)
	#P_E = OnePlanet.copy()
	#P_V = OnePlanet.copy()

	x_euler = np.zeros(N)
	y_euler = np.zeros(N)
	x_verlet = np.zeros(N)
	y_verlet = np.zeros(N)
	x_euler[0] = P_E[0].r[0]
	y_euler[0] = P_E[0].r[1]
	x_verlet[0] = P_V[0].r[0]
	y_verlet[0] = P_V[0].r[1]

	for i in range(1,N):
		accel_E.Update(P_E)
		solver.EulerCromer(P_E)
		x_euler[i] = P_E[0].r[0]
		y_euler[i] = P_E[0].r[1]
	plt.plot(x_euler,y_euler,label="EulerCromer")

	for i in range(1,N):
		accel_V.Update(P_V)
		solver.VelocityVerlet(P_V)
		x_verlet[i] = P_V[0].r[0]
		y_verlet[i] = P_V[0].r[1]
	plt.plot(x_verlet,y_verlet,label="VelocityVerlet")
	plt.legend()
	plt.show()

ComapareSolvers(h)

def CompareSolvers_2(Planets,h):
	accel = Acceleration(Planets)
	solver = DiffEqSolver(h)
	P_V = Planets.copy()
	P_E = Planets.copy()
	n_planets = len(Planets)
	x_euler = np.zeros((n_planets, N)) ## array for x positions
	y_euler = np.zeros((n_planets, N)) ## array for y positions
	x_verlet = np.zeros((n_planets, N)) ## array for x positions
	y_verlet = np.zeros((n_planets, N)) ## array for y positions

	# Calulate and plot position with velocity verlet: 
	for j in range(n_planets):
		x_verlet[j,0] = P_V[j].r[0]
		y_verlet[j,0] = P_V[j].r[1]
		
		for i in range(1,N):
			accel.Update(P_V)
			solver.VelocityVerlet(P_V)
			x_verlet[j,i] = P_V[j].r[0]
			y_verlet[j,i] = P_V[j].r[1]
		plt.plot(x_verlet[j],y_verlet[j],label="VelocityVerlet")

	# Calulate and plot position with euler cromer: 
	for j in range(n_planets):
		x_euler[j,0] = P_E[j].r[0]
		y_euler[j,0] = P_E[j].r[1]
		
		for i in range(1,N):
			accel.Update(P_E)
			solver.EulerCromer(P_E)
			x_euler[j,i] = P_E[j].r[0]
			y_euler[j,i] = P_E[j].r[1]
		plt.plot(x_euler[j],y_euler[j],label="EulerCromer")
		plt.legend()
		plt.show()
#CompareSolvers_2(OnePlanet,h)

# A list of all the planets
Planets = [Earth, Jupiter, Mercury, Venus, Mars, Saturn, Uranus, Neptune]

def PlotPlanets(Planets,N,h): #solves the coupled equations and plots the position of the planets

	solver = DiffEqSolver(h)
	accel = Acceleration(Planets)
	n_planets = len(Planets)

	x = np.zeros((n_planets, N)) ## array for x positions
	y = np.zeros((n_planets, N)) ## array for y positions

	# Calulate an plot position: 
	for j in range(n_planets):
		x[j,0] = Planets[j].r[0]
		y[j,0] = Planets[j].r[1]
		
		for i in range(1,N):
			accel.Update(Planets)
			solver.VelocityVerlet(Planets)
			x[j,i] = Planets[j].r[0]
			y[j,i] = Planets[j].r[1]
			
		plt.plot(x[j],y[j],label=Planets[j].name)
		plt.legend()
		
	#plt.xrange([-15,15])
	#plt.yrange([-15,15])
	#plt.axes('equal')
	plt.show()
	#plt.savefig('solar_system',dpi=225)

#PlotPlanets(Planets,N,h)