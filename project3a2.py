import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import linalg
import sys
import time


### Initial values:
pi=np.pi
r = 1 			# radius in AU
v = 2*pi/r 		# speed in AU/yr
a = 4*pi**2 	# acceleration in AU/yr**2

t = 1. 			# t is number of years
N = 10        # N is number of grid points
h = t/N 		# h is step length

## Vectors for x and y position: 
x = np.zeros(N+1)
y = np.zeros(N+1)
vx = np.zeros(N+1)
vy = np.zeros(N+1)
x[0] = r
y[0] = 0
vx[0] = 0
vy[0] = v
radius_dev = 0
ortho_dev = 0 
## Euler Cromer for calculating position and velocity: 
for i in range(N):
	vx[i+1] = vx[i]-h*a*x[i]/r
	vy[i+1] = vy[i]-h*a*y[i]/r
	x[i+1] = x[i]+h*vx[i+1]
	y[i+1] = y[i]+h*vy[i+1]
	radius_dev = max(np.abs(r-(x[i]**2+y[i]**2)),radius_dev)	# to check the radius
	ortho_dev = max(vx[i]*x[i] +vy[i]*y[i], ortho_dev ) 		# to check orthogonality
print(ortho_dev)
print(radius_dev)
print(x[0]-x[N], y[0]-y[N]) 					# to check end point deviation
print(vx[0]**2+vy[0]**2 -(vx[N]**2+ vy[N]**2))	# to check v**2 deviation
axis = np.linspace(0,r, num =N)
plt.plot(x, y)
plt.axis('equal')
plt.title("Sun-Earth system. The sun is in the origin. Earth in circular orbit.")
plt.xlabel("x/AU")
plt.ylabel("y/AU")
plt.show()
#plt.savefig('sun-earth',dpi=225)

#k = 2*t**2*np.pi**2/N**2
#for i in range(N-1):
#	x[i+1] = x[i] + (t*2*pi/N)*vx[i]/v - k*x[i]
#	y[i+1] = y[i] + (t*2*pi/N)*vy[i]/v - k*y[i]
#	vx[i+1] = vx[i] - k*(x[i+1] + x[i])
#	vy[i+1] = vy[i] - k*(y[i+1] + y[i])



#print(abs(x[-1]-x[0]))
#print(t**2/N)
