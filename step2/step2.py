#/usr/bin/python
#nonlinear convection
#
# du/dt + u du/dx = 0 
# same discretization scheme as step 1
import numpy
import matplotlib.pyplot
import time, sys

nx = 81
dx = 2./(nx-1)
nt = 20
dt = 0.01

u = numpy.ones(nx)
u[0.5/dx : 1/dx+1] = 2 #same ICs

un = numpy.ones(nx)

for n in range(nt):
    un = u.copy()
    for i in range(1,nx):
        u[i] = un[i] - un[i]*dt/dx*(un[i]-un[i-1])

matplotlib.pyplot.plot(numpy.linspace(0,2,nx), u)
matplotlib.pyplot.savefig("step2")
