
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

def plot(Planets, x, y, delta_t): # Input: list of planets, array of x-positions, array of y-positions, delta_t: distance between points in the plot
	N_plot = len(x[0])//delta_t
	N_p = len(Planets)
	X = np.zeros(N_plot)
	Y = np.zeros(N_plot)
	

	for j in range(N_p):
		for k in range(N_plot): 	#samples N_plot points from the solution vectors x[j] and y[j]
			X[k] = x[j,k*delta_t]
			Y[k] = y[j,k*delta_t]
	
			
		plt.plot(x[j],y[j],label=Planets[j].name)
	plt.title("The solar system. t_final = 20. h = 1/100" )
	plt.xlabel("x/AU")
	plt.ylabel("y/AU")
	plt.axis('equal')
	plt.legend()
	plt.savefig('solar_system20',dpi=225)
	plt.show()

# Initial values of the planets from 2017-Oct-04:
Mercury = Planet(1.65e-7, [-0.3789210,2.56223e-02,0],[-7.32704e-03*days, -2.68813e-02*days,0], "Mercury")
MercuryNewtonian = Planet(1.65e-7, [-0.3789210,2.56223e-02,0],[-7.32704e-03*days, -2.68813e-02*days,0], "MercuryNewtonian")
Venus = Planet(2.45e-6, [-.483111, .534117,0],[-1.49738e-02*days, -1.37846e-02*days,0 ], "Venus")
Earth = Planet(3e-6,[.985157,.191737e-01,0],[-3.48636e-3*days,1.68385e-02*days,0],"Earth")
Mars = Planet(3.3e-7, [-1.49970, 0.724618,0], [-5.52182e-03*days, -1.14214e-02*days,0],"Mars")
Jupiter = Planet(9.5e-4,[-4.63556,-2.8425,0],[3.85608e-03*days,-6.075643e-03*days,0],"Jupiter") 
Saturn = Planet(2.75e-4,[-0.421224,-10.0462,0],[5.26820e-03*days,-2.52167e-04*days,0], "Saturn")
Uranus = Planet(4.4e-5,[17.8825, 8.76632,0],[-1.75983e-03*days,3.34820e-03*days],"Uranus")
Neptune = Planet(0.515e-4, [28.6020, -8.86057, 0],[9.08471e-04*days, 3.01748e-03*days,0],"Neptune")


# A list of all the planets 
Planets = [Mercury,Venus, Earth,Mars,  Jupiter, Saturn, Uranus, Neptune]
innerPlanets = [Mercury,Venus, Earth,Mars]
NoMercury = [Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune]

print("Mercury", Mercury.total_energy(), "Earth", Earth.total_energy(), "Neptune", Neptune.total_energy())
time0 = time.time()

x,y = solve(Planets, "VelocityVerlet", 1/100, 20)
plot(Planets, x,y,10)

print("Mercury", Mercury.total_energy(), "Earth", Earth.total_energy(), "Neptune", Neptune.total_energy()) 
time = time.time()-time0
print(time)
