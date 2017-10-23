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
	
#	plt.plot(x_euler[0],y_euler[0],label="EulerCromer")
#	plt.plot(x_verlet[0],y_verlet[0],label="VelocityVerlet")
#	plt.title("Solutions of the Sun-Earh system, h=%.2f" %h)
#	plt.xlabel("x/AU")
#	plt.ylabel("y/AU")
#	plt.axis('equal')
#	plt.legend()
#	plt.show()
#	plt.savefig('VV_vs_EC2',dpi=225)
	return(diff_en_E, diff_en_V,diff_l_E, diff_l_V)

(a,b,c,d) = ComapareSolvers(1/20, 1.)


def plot_difference(H, t_final =1.): #does not worK!!!

	diff_en_E = np.zeros(H)
	diff_en_V= np.zeros(H)
	diff_l_E= np.zeros(H)
	diff_l_V= np.zeros(H)
	axis = np.zeros(H)
	for i in range(1,H):
		h = 1/(10*i)
		N = np.int(t_final/h)
		axis[i-1] = h
		diff_en_E[i-1], diff_en_V[i-1],diff_l_E[i-1], diff_l_V[i-1] = ComapareSolvers(h)
	plt.plot(np.log10(axis), np.log10(diff_en_E), np.log10(axis), np.log10(diff_en_V))
	#plt.xlabel("h")
	#plt.ylabel("relative difference in energy")
	plt.show()
	#plt.plot(axis, diff_l_E, axis, diff_l_V)
	#plt.xlabel("h")
	#plt.ylabel("relative difference in angular_momentum")
	#plt.show()


plot_difference(20)

