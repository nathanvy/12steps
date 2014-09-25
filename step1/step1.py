#/usr/bin/python

#for the 1D wave equation, del(u)/del(t) + c*del(u)/del(x)=0, 
import numpy
import matplotlib.pyplot
import time, sys

#discretizing using forward difference scheme for the time derivative and backward differencing for the space derivative,
# on a 1D grid i = 0...N we have,
# (u_n+1_i - u_n_i)/deltaT + c * (u_n_i - u_n_i-i)/deltaX = 0
nx = 81 # num points
dx = 2./(nx-1) #point spacing

nt = 25 #num time steps
dt = 0.025 # delta t
c = 1 #wave speed

# ICs
# u_0 = 2 on 0.5 <= x <=1 and u=1 everywhere else

u = numpy.ones(nx) # fill array u with 1's
u[0.5/dx : 1/dx+1]=2

print u

#temporary array for the next time step
un = numpy.ones(nx)

for n in range (nt): #loop from n=0 to n=nt
    un = u.copy() #copy existing values from u

    #for i in range(nx):
    for i in range (1,nx):
        u[i] = un[i] - c*dt/dx*(un[i]-un[i-1])

matplotlib.pyplot.plot(numpy.linspace(0,2,nx), u)
matplotlib.pyplot.savefig('step1')
