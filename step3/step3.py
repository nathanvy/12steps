#/usr/bin/python
#1D diffusion
#
# du/dt + nu d2u/dx2 = 0 
# forward time differencing, central 2nd order space differencing
import numpy
import matplotlib.pyplot
import time, sys

nx = 41
dx = 2./(nx-1)
nt = 20
sigma = 0.2
nu = 0.3
dt = sigma*(dx**2)/nu

u = numpy.ones(nx)
u[0.5/dx : 1/dx+1] = 2 #same ICs

un = numpy.ones(nx)

for n in range(nt):
    un = u.copy()
    for i in range(1,nx-1):
        u[i] = un[i] + nu*dt/dx**2 *(un[i+1] - 2*un[i] + un[i-1])

matplotlib.pyplot.plot(numpy.linspace(0,2,nx), u)
matplotlib.pyplot.savefig("step3")
