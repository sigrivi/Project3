import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import linalg
import sys
import time


### Initial values:
a = 100
v = 7
r = 1

t_final = 5
N = 100        # number of grid points

## Constants:
h = t_final/N # step length
k = h**2*a/(2*r)

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
	x[i+1] = x[i] + h*vx[i] - k*x[i]
	y[i+1] = y[i] + h*vy[i] - k*y[i]
	vx[i+1] = vx[i] - k*(x[i+1] + x[i])
	vy[i+1] = vy[i] - k*(y[i+1] + y[i])

axis = np.linspace(0,r, num =N)
plt.plot(axis, x, axis, y)
plt.show()
