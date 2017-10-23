import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import linalg
import sys
import time

pi = np.pi
days = 365.25 			# number of days in a year
h=1/100
from project3b import Planet 		# attributes: mass, r, v, a_s, a, name, total_energy(), angular_momentum()
from project3b import Acceleration	# attributes: N_p, Update(Planets, beta = 3.)
from project3b import DiffEqSolver 	# attributes: h, EulerCromer(Planets), VelcityVerlet(Planets)
from project3b import solve 		# input: (Planets, method, h, t_final=1., beta = 3.)

# Initial values of the planets from 2017-Oct-04:
Mercury = Planet(1.65e-7, [-0.3789210,2.56223e-02,0],[-7.32704e-03*days, -2.68813e-02*days,0], "Mercury")
MercuryNewtonian = Planet(1.65e-7, [-0.3789210,2.56223e-02,0],[-7.32704e-03*days, -2.68813e-02*days,0], "MercuryNewtonian")

P_Newtonian = [MercuryNewtonian]
x_Newtonian, y_Newtonian = solve(P_Newtonian, "VelocityVerlet",h)
plt.plot(x_Newtonian[0], y_Newtonian[0])
P_Relativistic = [Mercury]
x_Relatvistic, y_Relativistic = solve(P_Relativistic, "VelocityVerlet",h)
plt.plot(x_Relatvistic[0], y_Relativistic[0])
plt.show()