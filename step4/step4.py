#/usr/bin/python
#Burgers' eqn
#in 1D, looks something like
# du/dt + u du/dx = v d2u/dx2 which is a combination of nonlinear convection and diffusion
#
#using the same discretization scheme:
# (u_n+1_i - u_n_i)/deltaT + u_n_i(u_n_i - u_n_i-1)/deltaX = nu(u_n_i+i - 2u_n_i + u_n_i-1 )/deltaX^2
# which is forward differencing for time, rearward differenicing for space, and central 2nd-order differencing for the diffusion term
# therefore since we know the current space field, the only unknown is 
# u_n+1_i = u_n_i - u_n_i(deltaT/deltaX)(u_n_i-u_n_i-1) + v(deltaT/deltaX^2)(u_n_i+1 - 2u_n_i + u_n_i-1)
#
#our IC for this IVP will be: u = -2v/phi dphi/dx + 4
# where phi = exp(-x^2 / 4v) + exp(-(x-2pi)^2 / 4v)
#
#and our (periodic) BC will be u(0) = u(2pi)
import numpy
import sympy #symbolics lib
import matplotlib.pyplot
import time, sys

from sympy.utilities.lambdify import lambdify

# define symbols for sympy
x, nu, t = sympy.symbols('x nu t')
phi = sympy.exp( -(x-4*t)**2/(4*nu*(t+1)) ) + sympy.exp( -(x-4*t-2*numpy.pi)**2/(4*nu*(t+1)))
phiprime = phi.diff(x)
print phiprime

#debug shit
u = -2*nu*(phiprime/phi)+4
print u

ufunc = lambdify((t, x, nu), u)
print ufunc(1,4,3)


nx = 101
dx = 2*numpy.pi/(nx-1)
nt = 100
nu = 0.07
dt = dx*nu

x = numpy.linspace(0, 2*numpy.pi, nx)
#u = numpy.empty(nx)
un = numpy.empty(nx)
t = 0

u = numpy.asarray([ufunc(t,x0,nu) for x0 in x])

print u


matplotlib.pyplot.figure(figsize=(11,7), dpi=100)
matplotlib.pyplot.plot(x,u, marker = 'o', lw=2)
matplotlib.pyplot.xlim([0,2*numpy.pi])
matplotlib.pyplot.ylim([0,10])
matplotlib.pyplot.savefig("step4")
