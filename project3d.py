import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import linalg
import sys
import time

pi = np.pi
days = 365.25 			# number of days in a year

from project3b import Planet 		# attributes: mass, r, v, a_s, a, name, total_energy(), angular_momentum()
from project3b import Acceleration	# attributes: N_p, Update(Planets, beta = 3.)
from project3b import DiffEqSolver 	# attributes: h, EulerCromer(Planets), VelcityVerlet(Planets)
from project3b import solve 		# input: (Planets, method, h, t_final=1., beta = 3.)


def InitialVelocity(velocities,h, t_final = 3.):  # Input: a vector of initial x-velocities in multiples of pi, h: step length, t_final
	for v in velocities:
		Earth = [Planet(3e-6, [1,0,0], [0,v*pi,0], "Earth")]
		print(Earth[0].total_energy())

		x,y = solve(Earth, "VelocityVerlet", h, t_final)
		x,y = x[0],y[0]
		
		plt.plot(x,y, label =  "v= %.2f pi" %v)
	plt.legend()
	plt.axis('equal')
	plt.xlabel("x/AU")
	plt.ylabel("y/AU")
	plt.title("The orbit of Earth with different inital velocities. t_final = %d yr " %t_final)
	#plt.savefig('%d year' %t_final, dpi=225)
	plt.show()

velocities = [2,2.5,2**1.5, 3, 3.5]
InitialVelocity(velocities, 1/100, 1)

def GravitationalConstant(power, h, t_final = 1.):
	for p in power:
		Earth = [Planet(3e-6, [1,0,0], [0,2*pi,0], "Earth")]
		x,y = solve(Earth, "VelocityVerlet", h, t_final,  p)
		x,y = x[0],y[0]
		plt.plot(x,y, label = "beta = %.1f" %(p-1))
	plt.legend()
	plt.axis('equal')
	plt.xlabel("x/AU")
	plt.ylabel("y/AU")
	plt.title("The orbit of Earth with different expressions for gravitational force")
	#plt.savefig("gravitationalForce")
	plt.show()
GravitationalConstant([3,3.5, 4], 1/100, 5)