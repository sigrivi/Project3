import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import linalg
import sys
import time


### Initial values:
pi=np.pi
r = 1
v =  2*pi*r

a=4*pi**2
#v**2/r

t = 1.
N = 1000        # number of grid points

## Constants:
h = t/N # step length
T=60*60*24*365
k = 2*t**2*np.pi**2/N**2


## 

x = np.zeros(N)
y= np.zeros(N)
vx = np.zeros(N)
vy = np.zeros(N)
x[0] = r
y[0] = 0
vx[0] = 0
vy[0] = v

for i in range(N-1):
	
	vx[i+1] = vx[i]-h*a*x[i]/r**3
	vy[i+1] = vy[i]-h*a*y[i]/r**3
	x[i+1] = x[i]+h*vx[i+1]
	y[i+1] = y[i]+h*vy[i+1]

print(vx[1], vy[1])
print(x[1],y[1])
#for i in range(N-1):
#	x[i+1] = x[i] + (t*2*pi/N)*vx[i]/v - k*x[i]
#	y[i+1] = y[i] + (t*2*pi/N)*vy[i]/v - k*y[i]
#	vx[i+1] = vx[i] - k*(x[i+1] + x[i])
#	vy[i+1] = vy[i] - k*(y[i+1] + y[i])

#axis = np.linspace(0,r, num =N)
#plt.plot(x, y)
#plt.axis('equal')
#plt.show()

#print(abs(x[-1]-x[0]))
#print(t**2/N)
