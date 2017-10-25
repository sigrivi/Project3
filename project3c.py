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
	time0 = time.time()
	

	x_verlet,y_verlet = solve(P_V, "VelocityVerlet", h, t_final)	# solve with Velocity Verlet
	time_verlet = time.time()-time0

	time0 = time.time()
	x_euler,y_euler = solve(P_E, "EulerCromer", h, t_final) 		# solve with Euler Cromer
	time_euler = time.time()-time0

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
	return(diff_en_E, diff_en_V,diff_l_E, diff_l_V,time_euler, time_verlet)

#(a,b,c,d,e,f) = ComapareSolvers(1/20, 1.)


def plot_difference(H, t_final =1.): #plots the energy difference for h=1/(10H) to h=1/H

	diff_en_E = np.zeros(H)
	diff_en_V= np.zeros(H)
	diff_l_E= np.zeros(H)
	diff_l_V= np.zeros(H)
	axis = np.zeros(H)
	time_E = 0
	time_V = 0
	for i in range(1,H):
		h = 1/(10*i)
		N = np.int(t_final/h)
		axis[i-1] = h
		diff_en_E[i-1], diff_en_V[i-1],diff_l_E[i-1], diff_l_V[i-1], time_E, time_V = ComapareSolvers(h)
		print(N,"%.5f" %time_E,"%.5f" %time_V, time_V/time_E)

	
	plt.plot(np.log10(axis), np.log10(diff_en_E), label = "EulerCromer")
	plt.plot(np.log10(axis), np.log10(diff_en_V), label = "VelocityVerlet")
	plt.title("Relative difference in total energy as a function of step length h")
	plt.xlabel("log(h)")
	plt.ylabel("log(delta E)")
	plt. legend()
	plt.savefig('energy_loss',dpi=225)
	plt.show()



plot_difference(20)

