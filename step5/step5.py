#/usr/bin/python
#
#in 2d space, x_i = x0 + i*deltaX and y_i = y0 + i*deltaY
#now u_i,j = u(x_i, y_i)
#which makes e.g. delu/delx at i,j = (u_i+1,j - u_i,j)/deltaX + O(deltaX)
#
#for 2d linear convection:
# du/dt + cdu/dx + cdu/dy = 0
#
# discretized,
# (u_n+1_i,j - u_n_i,j)/deltaT + c(u_n_i,j - u_n_i-1,j)/deltaX + c(u_n_i,j - u_n_i,j-1)/deltaY = 0
# where again we know the present flow field so the only unknown is u_n+1_i,j
#
#ICs: u(x) = 2 for 0.5<= x <= 1, 1 for everywhere else
#BCs: u = 1 for x=0, x=2, and y=0, y=2

from mpl_toolkits.mplot3d import Axes3D ## for 3D projected plots

import numpy
import matplotlib.pyplot

#setup
nx = 81
ny = 81
nt = 100
c = 1
dx = 2.0/(nx-1)
dy = 2.0/(ny-1)
sigma = 0.2
dt = sigma*dx

x = numpy.linspace(0,2,nx)
y = numpy.linspace(0,2,ny)

u = numpy.ones((ny, nx)) # creates a 1 x n vector of 1's
un = numpy.ones((ny, nx))

# ICs
u[0.5/dy:1/dy+1, 0.5/dx:1/dx+1] = 2

# plot that shit
#fig = matplotlib.pyplot.figure(figsize=(11,7), dpi = 100)
#ax = fig.gca(projection='3d')
#X, Y = numpy.meshgrid(x,y)
#surf = ax.plot_surface(X,Y,u[:])
#matplotlib.pyplot.savefig("step5-IC")

#shitty nested-loop method
#for n in range(nt+1):
#    un = u.copy()
#    for i in range(1, len(u)):
#        for j in range(1, len(u)):
#            u[i,j] = un[i,j] - (c*dt/dx*(un[i,j] - un[i-1,j])) - (c*dt/dy*(un[i,j]-un[i,j-1]))
#            u[0,:] = 1
#            u[-1,:] = 1
#            u[:,0] = 1
#            u[:,-1] = 1

#or, faster array method
for n in range(nt+1):
    un[:] = u[:]
    u[1:,1:] = un[1:,1:] - (c*dt/dx*(un[1:,1:] - un[0:-1,1:])) - (c*dt/dy*(un[1:,1:] - un[1:,0:-1]))
    u[0,:] = 1
    u[-1,:] = 1
    u[:,0] = 1
    u[:,-1] = 1

fig = matplotlib.pyplot.figure(figsize=(11,7), dpi = 100)
ax = fig.gca(projection='3d')
X, Y = numpy.meshgrid(x,y)
surf2 = ax.plot_surface(X,Y,u[:])
matplotlib.pyplot.savefig("step5")
