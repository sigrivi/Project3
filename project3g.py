import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import linalg
import sys
import time

pi = np.pi
days = 365.25 			# number of days in a year
h=1/(360*60*60*8)
from project3b import Planet 		# attributes: mass, r, v, a_s, a, name, total_energy(), angular_momentum()
from project3b import Acceleration	# attributes: N_p, Update(Planets, beta = 3.)
from project3b import DiffEqSolver 	# attributes: h, EulerCromer(Planets), VelcityVerlet(Planets)
from project3b import solve 		# input: (Planets, method, h, t_final=1., beta = 3.)


# Initial values of the planets from 2017-Oct-04:
Mercury = Planet(1.65e-7, [-0.3789210,2.56223e-02,0],[-7.32704e-03*days, -2.68813e-02*days,0], "Mercury")
MercuryNewtonian = Planet(1.65e-7, [-0.3789210,2.56223e-02,0],[-7.32704e-03*days, -2.68813e-02*days,0], "MercuryNewtonian")

# A list of all the planets 
OnlyMercury = [Mercury]
OnlyMercury_Newtonian = [MercuryNewtonian]

x_Newtonian, y_Newtonian = solve(OnlyMercury_Newtonian, "VelocityVerlet",h, 100) #The newtonian solution of the Mercury-Sun system
#plt.plot(x_Newtonian[0], y_Newtonian[0])

x_Relatvistic, y_Relativistic = solve(OnlyMercury, "VelocityVerlet",h, 100) #The relativistic sloution of the Mercury-Sun system

for i in range(len(x_Newtonian[0])):#
	if (x_Newtonian[0,i]**2+y_Newtonian[0,i]**2 - 0.3075**2)<1.e-7: # search for the perihelion (perihelion distance : 0.3075 AU)
		print(np.arctan(y_Newtonian[0,i]/x_Newtonian[0,i]), np.arctan(y_Relativistic[0,i]/x_Relatvistic[0,i])) # calculate perihelion angle


plt.plot(x_Relatvistic[0], y_Relativistic[0])
plt.axis('equal')
plt.title("Mercury")
plt.xlabel("x/AU")
plt.ylabel("y/AU")
plt.show()
#x_N = np.zeros(250)
#y_N= np.zeros(250)
#x_Nf = np.zeros(250)
#y_Nf = np.zeros(250)
#x_R=np.zeros(250)
#y_R= np.zeros(250)
#x_Rf =np.zeros(250)
#y_Rf=np.zeros(250)

#for i in range(250):
#	k = 100000-250
#	x_N[i],y_N[i] = x_Newtonian[0,k+i],y_Newtonian[0,k+i]
#	x_Nf[i],y_Nf[i] = x_Newtonian[0,k+i],y_Newtonian[0,k+i]
#	x_R[i],y_R[i] = x_Relatvistic[0,k+i],y_Relativistic[0,k+i]
#	x_Rf[i],y_Rf[i] = x_Relatvistic[0,k+i],y_Relativistic[0,k+i]

#plt.plot(x_N,y_N, label = "Newtonian0")
#plt.plot(x_Nf,y_Nf, label = "NewtonianFinal")
#plt.plot(x_R,y_R, label = "Relativistic0")
#plt.plot(x_Rf, y_Rf, label = "Relativisticfinal")
#plt.legend()
#plt.show()
