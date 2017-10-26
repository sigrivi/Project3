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
		self.l = (self.r[0]*self.v[1]-self.r[1]*self.v[0])
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


	# Calulate position: 
	for j in range(N_p):
		x[j,0] = Planets[j].r[0]
		y[j,0] = Planets[j].r[1]
		
		if method == "EulerCromer":
			for i in range(1,N):
				accel.Update(Planets)
				solver.EulerCromer(Planets)
				x[j,i] = Planets[j].r[0]
				y[j,i] = Planets[j].r[1]
				#print("EulerCromer",np.linalg.norm(Planets[j].r),np.linalg.norm(Planets[j].v),np.linalg.norm(np.cross(Planets[j].r,Planets[j].v)))
				#print(Planets[j].angular_momentum())
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

